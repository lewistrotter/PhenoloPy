# phenolop
'''
This script contains functions for calculating per-pixel phenology metrics (phenometrics)
on a timeseries of vegetation index values (e.g. NDVI) stored in a xarray DataArray. The
methodology is mostly based on the methodology underpinning TIMESAT 3.3. Some code has also
been adapted from the great work of Chad from DEAfrica.

The script contains the following primary functions:
    1. 
    2. calc_vege_index: calculate one of several vegetation indices;
    3. resample: resamples data to a defined time interval (e.g. bi-monthly);
    4. removal_outliers: detects and removes timeseries outliers;
    5. interpolate: fill in missing nan values;
    6. smooth: applies one of various smoothers to raw vegetation data;
    7. correct_upper_envelope: shift timeseries upwards closer to upper envelope;
    8. detect_seasons: count num of seasons (i.e. peaks) in timeseries;
    9. calc_phenometrics: generate phenological metrics on timeseries.

The script also contains the following secondary functions:
    1. conform_dea_band_names: 

Links:
TIMESAT 3.3: http://web.nateko.lu.se/timesat/timesat.asp

Contacts: 
Lewis Trotter: lewis.trotter@postgrad.curtin.edu.au
'''

# import required libraries
import os, sys
import xarray as xr
import numpy as np
import pandas as pd
import math
import dask
from datetime import datetime, timedelta
from scipy.stats import zscore
from scipy.signal import savgol_filter, find_peaks
from scipy.ndimage import gaussian_filter
from statsmodels.tsa.seasonal import STL as stl
from datacube.utils.geometry import assign_crs


def conform_dea_band_names(ds):
    """
    Takes an xarray dataset containing spectral bands from various different
    satellite digital earth australia (DEA) products and conforms band names.
    Only use if you are using satellite data from DEA.
    
    Parameters
    ----------
    ds: xarray Dataset
        A two-dimensional or multi-dimensional array containing spectral 
        bands that will be renamed.

    Returns
    -------
    ds : xarray Dataset
        The original xarray Dataset inputted into the function, with renamed 
        spectral bands.
    """
    
    # notify user
    print('Conforming satellite bands')
    
    # check if dataset type
    if type(ds) != xr.Dataset:
        raise TypeError('Not a dataset. Please provide a xarray dataset.')
    
    # create band rename mapping dictionary
    band_map_dict = {
        'nbart_blue': 'blue',
        'nbart_green': 'green',
        'nbart_red': 'red',
        'nbart_nir': 'nir',
        'nbart_swir_1': 'swir1',
        'nbart_swir_2': 'swir2',
        'nbar_blue': 'blue',
        'nbar_green': 'green',
        'nbar_red': 'red',
        'nbar_nir': 'nir',
        'nbar_swir_1': 'swir1',
        'nbar_swir_2': 'swir2',
        'nbart_nir_1': 'nir',
        'nbart_red_edge_1': 'red_edge_1', 
        'nbart_red_edge_2': 'red_edge_2',    
        'nbart_swir_2': 'swir1',
        'nbart_swir_3': 'swir2',
        'nbar_nir_1': 'nir',
        'nbar_red_edge_1': 'red_edge_1',   
        'nbar_red_edge_2': 'red_edge_2',   
        'nbar_swir_2': 'swir1',
        'nbar_swir_3': 'swir2'
    }
 
    # rename bands in dataset to use conformed naming conventions
    bands_to_rename = {
        a: b for a, b in band_map_dict.items() if a in ds.variables
    }
    
    # apply the rename
    ds = ds.rename(bands_to_rename)
    
    # notify user
    print('> Satellite band names conformed successfully.\n')
    
    return ds


def calc_vege_index(ds, index='ndvi', drop=True):
    """
    Takes an xarray dataset containing spectral bands, calculates one of
    a set of remote sensing indices, and adds the resulting array as a 
    new variable called 'veg_index' in the original dataset.  
    
    Parameters
    ----------
    ds: xarray Dataset
        A two-dimensional or multi-dimensional array containing the spectral 
        bands required to calculate the index.
    index: str
        A string giving the name of the index to calculate:
        ndvi: Normalised Difference Vegetation Index (Rouse 1973)
        evi:  Enhanced Vegetation Index (Huete 2002)
        mavi: Moisture-Adjusted Vegetation Index
    drop: bool
        Provides the option to drop the original input variables, thus saving 
        space. If drop = True, returns only the index and its values.

    Returns
    -------
    ds : xarray Dataset
        The original xarray Dataset inputted into the function, with a 
        new varible named 'veg_index' containing the remote sensing index 
        as a DataArray. If drop = True, the new variable/s as DataArrays in 
        the original Dataset. 
    """
    
    # notify user
    print('Generating vegetation index: {0}'.format(index))
    
    # check if dataset type
    if type(ds) != xr.Dataset:
        raise TypeError('Not a dataset. Please provide a xarray dataset.')   
        
    # check for valid bands
    valid_bands = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'red_edge_1', 'red_edge_2']
    if not set(ds.variables) & set(valid_bands):
        raise ValueError('No valid spectral band names (e.g. red) exist in dataset. Aborting.')
        
    # get input band variables in order to drop these if drop=True
    if drop:
        bands_to_drop = list(ds.data_vars)
        print('> Drop bands set to True. Dropping these bands: {0}'.format(bands_to_drop))
        
    # calculate vegetation index
    if index is None:
        raise ValueError('No remote sensing index provided. Please provide an index.')
    elif index == 'ndvi':
        ds['veg_index'] = (ds.nir - ds.red) / (ds.nir + ds.red)
    elif index == 'evi':
        ds['veg_index'] = (2.5 * (ds.nir - ds.red)) / (ds.nir + 6 * ds.red - 7.5 * ds.blue + 1)
    elif index == 'mavi':
        ds['veg_index'] = (ds.nir - ds.red) / (ds.nir + ds.red + ds.swir1)
    else:
        raise ValueError('Selected index is not a valid index.')
        
    # drop no-longer required bands
    if drop:
        ds = ds.drop(bands_to_drop)
        
    # notify user
    print('> Vegetation index calculated successfully.\n')
        
    return ds


def remove_outliers(ds, method='median', user_factor=2, z_pval=0.05):
    """
    Takes an xarray dataset containing vegetation index variable and removes outliers within 
    the timeseries on a per-pixel basis. The resulting dataset contains the timeseries 
    with outliers set to nan. Based on methodology of TIMESAT. Can work on datasets with or
    without existing nan values.
    
    Parameters
    ----------
    ds: xarray Dataset
        A two-dimensional or multi-dimensional array containing a vegetation 
        index variable (i.e. 'veg_index').
    method: str
        The outlier detection method to apply to the dataset. The median method detects 
        outliers by calculating if values in pixel timeseries deviate more than a maximum 
        deviation (cutoff) from the median in a moving window (half window width = number 
        of values per year / 7) and it is lower than the mean value of its immediate neighbors 
        minus the cutoff or it is larger than the highest value of its immediate neighbor plus 
        The cutoff is the standard deviation of the entire time-series times a factor given by 
        the user. The second method, zscore, is similar but uses zscore to detect whether outlier
        is signficicantly (i.e. p-value) outside the population.
    user_factor: float
        An value between 0 to 10 which is used to 'multiply' the threshold cutoff. A higher factor 
        value results in few outliers (i.e. only the biggest outliers). Default factor is 2.
    z_pval: float
        The p-value for zscore method. A more significant p-value (i.e. 0.01) results in fewer
        outliers, a less significant p-value (i.e 0.1) results in more. Default is 0.05.
    Returns
    -------
    ds : xarray Dataset
        The original xarray Dataset inputted into the function, with a all detected outliers in the
        veg_index variable set to nan.
    """
    
    # notify user
    print('Outlier removal method: {0} with a user factor of: {1}'.format(method, user_factor))
    
    # check if type is xr dataset
    if type(ds) != xr.Dataset:
        raise TypeError('Not a dataset. Please provide a xarray dataset.')
        
    # check if time dimension is in dataset
    if 'time' not in list(ds.dims):
        raise ValueError('Time dimension not in dataset. Please ensure dataset has a time dimension.')
    
    # check if dataset contains veg_index variable
    if 'veg_index' not in list(ds.data_vars):
        raise ValueError('Vegetation index (veg_index) not in dataset. Please generate veg_index first.')
                        
    # check if dataset is 2D or above
    if len(ds['veg_index'].shape) == 1:
        raise Exception('Remove outliers does not operate on 1D datasets. Ensure it has an x, y and time dimension.')
        
    # check if user factor provided
    if user_factor <= 0:
        raise TypeError('User factor is less than 0. Please provide a value of 0 or above.')
        
    # check if pval provided if method is zscore
    if method == 'zscore' and z_pval not in [0.1, 0.05, 0.01]:
        raise ValueError('Zscore selected but invalid pvalue provided. Ensure it is either 0.1, 0.05 or 0.01.')    
        
    # remove outliers based on user selected method
    if method in ['median', 'zscore']:
        
        # calc cutoff val per pixel i.e. stdv of pixel multiply by user-factor 
        cutoffs = ds.std('time') * user_factor

        # generate outlier mask via median or zscore method
        if method == 'median':

            # calc mask of existing nan values (nan = True) in orig ds
            ds_mask = xr.where(ds.isnull(), True, False)

            # calc win size via num of dates in dataset
            win_size = int(len(ds['time']) / 7)
            win_size = int(win_size / int(len(ds.resample(time='1Y'))))

            if win_size < 3:
                win_size = 3
                print('> Generated roll window size less than 3, setting to default (3).')
            elif win_size % 2 == 0:
                win_size = win_size + 1
                print('> Generated roll window size is an even number, added 1 to make it odd ({0}).'.format(win_size))
            else:
                print('> Generated roll window size is: {0}'.format(win_size))

            # calc rolling median for whole dataset
            ds_med = ds.rolling(time=win_size, center=True).median()

            # calc nan mask of start/end nans from roll, replace them with orig vals
            med_mask = xr.where(ds_med.isnull(), True, False)
            med_mask = xr.where(ds_mask != med_mask, True, False)
            ds_med = xr.where(med_mask, ds, ds_med)

            # calc abs diff between orig ds and med ds vals at each pixel
            ds_diffs = abs(ds - ds_med)

            # calc mask of outliers (outlier = True) where absolute diffs exceed cutoff
            outlier_mask = xr.where(ds_diffs > cutoffs, True, False)

        elif method == 'zscore':

            # generate critical val from user provided p-value
            if z_pval == 0.01:
                crit_val = 2.3263
            elif z_pval == 0.05:
                crit_val = 1.6449
            elif z_pval == 0.1:
                crit_val = 1.2816
            else:
                raise ValueError('Zscore p-value not supported. Please use 0.1, 0.05 or 0.01.')

            # calc zscore, ignore nans in timeseries vectors
            zscores = ds.apply(zscore, nan_policy='omit', axis=0)

            # calc mask of outliers (outlier = True) where zscore exceeds critical value
            outlier_mask = xr.where(abs(zscores) > crit_val, True, False)

        # shift values left and right one time index and combine, get mean and max for each window
        lefts, rights = ds.shift(time=1).where(outlier_mask), ds.shift(time=-1).where(outlier_mask)
        nbr_means = (lefts + rights) / 2  # TODO - handle nans i.e. (1, np.nan) / 2 = np.nan... no!
        nbr_maxs = xr.ufuncs.fmax(lefts, rights)

        # keep nan only if middle val < mean of neighbours - cutoff or middle val > max val + cutoffs
        outlier_mask = xr.where((ds.where(outlier_mask) < (nbr_means - cutoffs)) | 
                                (ds.where(outlier_mask) > (nbr_maxs + cutoffs)), True, False)

        # flag outliers as nan in original da
        ds = xr.where(outlier_mask, np.nan, ds)
        
    else:
        raise ValueError('Provided method not supported. Please use median or zscore.')
        
    # check if any nans exist in dataset after resample and tell user
    if bool(ds.isnull().any()):
        print('> Warning: dataset contains nan values. You may want to interpolate next.')

    # notify user
    print('> Outlier removal successful.\n')

    return ds


def resample(ds, interval='1M', reducer='median'):
    """
    Takes an xarray dataset containing vegetation index variable and resamples
    to a new temporal resolution. The available time intervals are 1W (weekly),
    2W (bi-monthly) and 1M (monthly) resample intervals. The resulting dataset
    contains the new resampled veg_index variable.
    
    Parameters
    ----------
    ds: xarray Dataset
        A two-dimensional or multi-dimensional array containing a vegetation 
        index variable (i.e. 'veg_index').
    interval: str
        The new temporal interval which to resample the dataset to. Available
        intervals include 1W (weekly), 2W (bi-monthly) and 1M (monthly).
    reducer: str
        The reducer function to apply when downsampling. Can choose median or
        mean as the reducer.

    Returns
    -------
    ds : xarray Dataset
        The original xarray Dataset inputted into the function, with a 
        newly resampled 'veg_index' variable.
    """
    
    # notify user
    print('Resampling dataset interval: {0} via reducer: {1}'.format(interval, reducer))
    
    # check if type is xr dataset
    if type(ds) != xr.Dataset:
        raise TypeError('Not a dataset. Please provide a xarray dataset.')
        
    # check if time dimension is in dataset
    if 'time' not in list(ds.dims):
        raise ValueError('Time dimension not in dataset. Please ensure dataset has a time dimension.')
    
    # check if dataset contains veg_index variable
    if 'veg_index' not in list(ds.data_vars):
        raise ValueError('Vegetation index (veg_index) not in dataset. Please generate veg_index first.')

    # check if dataset is 2D or above
    if len(ds['veg_index'].shape) == 1:
        raise Exception('Resample does not operate on 1D datasets. Ensure it has an x, y and time dimension.')
        
    # resample based on user selected interval and reducer
    if interval in ['1W', '2W', '1M']:
        if reducer == 'mean':
            ds = ds.resample(time=interval).mean('time')
        elif reducer == 'median':
            ds = ds.resample(time=interval).median('time')
        else:
            raise ValueError('Provided reducer not supported. Please use mean or median.')
    else:
        raise ValueError('Provided resample interval not supported. Please use 1W, 2W or 1M.')
        
    # check if any nans exist in dataset after resample and tell user
    if bool(ds.isnull().any()):
        print('> Warning: dataset contains nan values. You should interpolate nan values next.')
        
    # notify user
    print('> Resample successful.\n')
    
    return ds


def group(ds, group_by='month', reducer='median'):
    """
    Takes an xarray dataset containing a vegetation index variable, groups and 
    reduces values based on a specified temporal group. The available group 
    options are by month only. The resulting dataset contains the new grouped 
    veg_index variable as a single year of weeks or months.
    
    Parameters
    ----------
    ds: xarray Dataset
        A two-dimensional or multi-dimensional array containing a vegetation 
        index variable (i.e. 'veg_index').
    group_by: str
        The groups which to reduce the dataset to. Available intervals only
        include month at this stage.
    reducer: str
        The reducer function to apply when downsampling. Can choose median or
        mean as the reducer.

    Returns
    -------
    ds : xarray Dataset
        The original xarray Dataset inputted into the function, with a 
        newly resampled 'veg_index' variable.
    """
    
    # notify user
    print('Group dataset interval: {0} via reducer: {1}'.format(group_by, reducer))
    
    # check if type is xr dataset
    if type(ds) != xr.Dataset:
        raise TypeError('Not a dataset. Please provide a xarray dataset.')
        
    # check if time dimension is in dataset
    if 'time' not in list(ds.dims):
        raise ValueError('Time dimension not in dataset. Please ensure dataset has a time dimension.')
    
    # check if dataset contains veg_index variable
    if 'veg_index' not in list(ds.data_vars):
        raise ValueError('Vegetation index (veg_index) not in dataset. Please generate veg_index first.')

    # check if dataset is 2D or above
    if len(ds['veg_index'].shape) == 1:
        raise Exception('Resample does not operate on 1D datasets. Ensure it has an x, y and time dimension.')
        
    # get all years in dataset, choose middle year in array
    years = np.array([year for year in ds.groupby('time.year').groups])
    year = np.take(years, years.size // 2)
    
    # notify user
    print('> Selecting year: {0} to re-label times after groupby.'.format(year))
          
    # group based on user selected interval and reducer
    if group_by in ['week', 'month']:
        if reducer == 'mean':
            ds = ds.groupby('time' + '.' + group_by).mean('time')
        elif reducer == 'median':
            ds = ds.groupby('time' + '.' + group_by).median('time')
        else:
            raise ValueError('Provided reducer not supported. Please use mean or median.')
    else:
        raise ValueError('Provided group_by interval not supported. Please use month.')
        
    # correct time index following group by
    if 'month' in list(ds.dims):
        ds = ds.rename({'month': 'time'})
        times = [datetime(year, int(dt), 1) for dt in ds['time']]
        ds['time'] = [np.datetime64(dt) for dt in times]
        
    elif 'week' in list(ds.dims):
        ds = ds.rename({'week': 'time'})
        times = [datetime.strptime('{0} {1} {2}'.format(year, int(dt), 0), '%Y %W %w') for dt in ds['time']]
        ds['time'] = [np.datetime64(dt) for dt in times]
        
    else:
        raise ValueError('Group_by was not found in dataset dimension. Aborting.')
    
    # check if any nans exist in dataset after resample and tell user
    if bool(ds.isnull().any()):
        print('> Warning: dataset contains nan values. You should interpolate nan values next.')
        
    # notify user
    print('> Group successful.\n')
    
    return ds

    
def interpolate(ds, method='interpolate_na'):
    """
    Takes an xarray dataset containing vegetation index variable and interpolates
    (linearly) all existing nan values within the timeseries using one of several
    methods. The resulting dataset contains the timeseries minus nan values.
    
    Parameters
    ----------
    ds: xarray Dataset
        A two-dimensional or multi-dimensional array containing a vegetation 
        index variable (i.e. 'veg_index').
    method: str
        The interpolation method to apply to the dataset to fill in nan values.
        Two methods are available. First is the built in xarray interpolate_na method, 
        which is robust but slow. The second is a custom DEA method called fast_fill, 
        which speeds up the process.
        
    Returns
    -------
    ds : xarray Dataset
        The original xarray Dataset inputted into the function, with a 
        newly interpolated 'veg_index' variable.
    """
    
    # notify user
    print('Interpolating dataset using method: {0}.'.format(method))
    
    # check if type is xr dataset
    if type(ds) != xr.Dataset:
        raise TypeError('Not a dataset. Please provide a xarray dataset.')
        
    # check if time dimension is in dataset
    if 'time' not in list(ds.dims):
        raise ValueError('Time dimension not in dataset. Please ensure dataset has a time dimension.')
    
    # check if dataset contains veg_index variable
    if 'veg_index' not in list(ds.data_vars):
        raise ValueError('Vegetation index (veg_index) not in dataset. Please generate veg_index first.')
                        
    # check if dataset is 2D or above
    if len(ds['veg_index'].shape) == 1:
        raise Exception('Interpolate does not operate on 1D datasets. Ensure it has an x, y and time dimension.')
        
    # resample based on user selected interval and reducer
    if method in ['interpolate_na', 'fast_fill']:
        
        if method == 'interpolate_na':
            
            # use internal xarray linear interpolate method along time dim
            ds = ds.interpolate_na(dim='time', method='linear')
            
        elif method == 'fast_fill':
            
            # grab x, y, time, etc. and reshape
            x, y, time, attrs = ds['veg_index'].x, ds['veg_index'].y, ds['veg_index'].time, ds['veg_index'].attrs
            da = ds['veg_index'].transpose("y", "x", "time").values
            
            # create nan mask and get indexes
            mask = np.isnan(da)
            idx = np.where(~mask, np.arange(mask.shape[-1]), 0)
            #np.maximum.accumulate(idx, axis=-1, out=idx)
            
            # build new grid as template
            i, j = np.meshgrid(np.arange(idx.shape[0]), np.arange(idx.shape[1]), indexing="ij")
            dat = da[i[:, :, np.newaxis], j[:, :, np.newaxis], idx]

            # if nan detected, fill it and add to template
            if np.isnan(np.sum(dat[:, :, 0])):
                fill = np.nanmean(dat, axis=-1)
                for t in range(dat.shape[-1]):
                    mask = np.isnan(dat[:, :, t])
                    if mask.any():
                        dat[mask, t] = fill[mask]
                    else:
                        break

            #stack back into da template
            dat = xr.DataArray(dat, attrs=attrs, 
                               coords={"x": x, "y": y, "time": time}, 
                               dims=["y", "x", "time"]
                              )

            # convert back to dataset
            ds = dat.to_dataset(name='veg_index')
            
    else:
        raise ValueError('Provided method not supported. Please use interpolate_na or fast_fill')
        
    # check if any nans exist in dataset after resample and tell user
    if bool(ds.isnull().any()):
        print('> Warning: dataset still contains nan values. The first and/or last time slices may be empty.')
        
    # notify user
    print('> Interpolation successful.\n')
    
    return ds


def smooth(ds, method='savitsky', window_length=3, polyorder=1, sigma=1):  
    """
    Takes an xarray dataset containing vegetation index variable and smoothes timeseries
    timeseries on a per-pixel basis. The resulting dataset contains a smoother timeseries. 
    Based on methodology of TIMESAT. Recommended that no nan values present in dataset.
    
    Parameters
    ----------
    ds: xarray Dataset
        A two-dimensional or multi-dimensional array containing a vegetation 
        index variable (i.e. 'veg_index').
    method: str
        The smoothing algorithm to apply to the dataset. The savitsky method uses the robust
        savitsky-golay smooting technique, as per TIMESAT. Symmetrical gaussian applies a simple 
        symmetrical gaussian. Asymmetrical gaussian applies an asymmetrical gaussian, resulting in
        a flatter peak. Double logistic applies two seperate logistic functions to give a flatter 
        peak based on TIMESAT. Default is savitsky.
    window_length: int
        The length of the filter window (i.e., the number of coefficients). Value must 
        be a positive odd integer. The larger the window length, the smoother the dataset.
        Default value is 3 (as per TIMESAT).
    polyorder: int
        The order of the polynomial used to fit the samples. Must be a odd number (int) and
        less than window_length.
    sigma: int
        Standard deviation for Gaussian kernel. The standard deviations of the Gaussian filter 
        must be provided as a single number between 1-9.
        
    Returns
    -------
    ds : xarray Dataset
        The original xarray Dataset as input into the function, with smoothed data in the
        veg_index variable.
    """
    
    # notify user
    print('Smoothing method: {0} with window length: {1} and polyorder: {2}.'.format(method, window_length, polyorder))
    
    # check if type is xr dataset
    if type(ds) != xr.Dataset:
        raise TypeError('> Not a dataset. Please provide a xarray dataset.')
        
    # check if time dimension is in dataset
    if 'time' not in list(ds.dims):
        raise ValueError('> Time dimension not in dataset. Please ensure dataset has a time dimension.')
    
    # check if dataset contains veg_index variable
    if 'veg_index' not in list(ds.data_vars):
        raise ValueError('> Vegetation index (veg_index) not in dataset. Please generate veg_index first.')
                        
    # check if dataset is 2D or above
    if len(ds['veg_index'].shape) == 1:
        raise Exception('> Remove outliers does not operate on 1D datasets. Ensure it has an x, y and time dimension.')
        
    # check if window length provided
    if window_length <= 0 or not isinstance(window_length, int):
        raise TypeError('> Window_length is <= 0 and/or not an integer. Please provide a value of 0 or above.')
        
    # check if user factor provided
    if polyorder <= 0 or not isinstance(polyorder, int):
        raise TypeError('> Polyorder is <= 0 and/or not an integer. Please provide a value of 0 or above.')
        
    # check if polyorder less than window_length
    if polyorder > window_length:
        raise TypeError('> Polyorder is > than window_length. Must be less than window_length.')
        
    # check if sigma is between 1 and 9
    if sigma < 1 or sigma > 9:
        raise TypeError('> Sigma is < 1 or > 9. Must be between 1 - 9.')
        
    # perform smoothing based on user selected method     
    if method in ['savitsky', 'symm_gaussian', 'asymm_gaussian', 'double_logistic']:
        if method == 'savitsky':
            
            # create savitsky smoother func
            def smoother(da, window_length, polyorder):
                return da.apply(savgol_filter, window_length=window_length, polyorder=polyorder, axis=0)
            
            # create kwargs dict
            kwargs = {'window_length': window_length, 'polyorder': polyorder}

        elif method == 'symm_gaussian':
            
            # create gaussian smoother func
            def smoother(da, sigma):
                return da.apply(gaussian_filter, sigma=sigma)
            
            # create kwargs dict
            kwargs = {'sigma': sigma}

        elif method == 'asymm_gaussian':
            raise ValueError('> Asymmetrical gaussian not yet implemented.')
            
        elif method == 'double_logistic':
            raise ValueError('> Double logistic not yet implemented.')
                
        # create template and map func to dask chunks
        temp = xr.full_like(ds, fill_value=np.nan)
        ds = xr.map_blocks(smoother, ds, template=temp, kwargs=kwargs)
        
    else:
        raise ValueError('Provided method not supported. Please use savtisky.')
        
    # check if any nans exist in dataset after resample and tell user
    if bool(ds.isnull().any()):
        print('> Warning: dataset contains nan values. You may want to interpolate next.')

    # notify user
    print('> Smoothing successful.\n')

    return ds


def calc_num_seasons(ds):
    """
    Takes an xarray Dataset containing vege values and calculates the number of
    seasons for each timeseries pixel. The number of seasons provides a count of number 
    of peaks in each pixel timeseries.

    Parameters
    ----------
    ds: xarray Dataset
        A two-dimensional or multi-dimensional array containing an Dataset of veg_index 
        and time values. 

    Returns
    -------
    da_num_seasons : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        number of seasons value detected across the timeseries at each pixel.
    """
    
    # notify user
    print('Beginning calculation of number of seasons.')
    
    # check if type is xr dataset
    if type(ds) != xr.Dataset:
        raise TypeError('> Not a dataset. Please provide a xarray dataset.')
        
    # check if time dimension is in dataset
    if 'time' not in list(ds.dims):
        raise ValueError('> Time dimension not in dataset. Please ensure dataset has a time dimension.')
    
    # check if dataset contains veg_index variable
    if 'veg_index' not in list(ds.data_vars):
        raise ValueError('> Vegetation index (veg_index) not in dataset. Please generate veg_index first.')
                        
    # check if dataset is 2D or above
    if len(ds['veg_index'].shape) == 1:
        raise Exception('> Remove outliers does not operate on 1D datasets. Ensure it has an x, y and time dimension.')
        
    # get height (val at 80% threshold) and dist between peals (4 seasons)
    #da_height = (ds['veg_index'].max('time') - ds['veg_index'].min('time')) * 0.80
    #distance = math.ceil(len(ds['time']) / 4)

    # set up calc peaks functions
    def calc_peaks(vec_data):
    
        # get height (val at 75% threshold) and dist between peals (4 seasons)
        height = np.nanquantile(vec_data, q=0.75)
        distance = math.ceil(len(vec_data) / 4)

        # calc num peaks
        peaks = find_peaks(vec_data, height=height, distance=distance)

        if peaks:
            return len(peaks[0])
        else:
            return 0
        
    # calculate nos using calc_funcs func
    print('> Calculating number of seasons.')
    da_nos = xr.apply_ufunc(calc_peaks, ds['veg_index'], 
                            input_core_dims=[['time']],
                            vectorize=True, 
                            dask='parallelized', 
                            output_dtypes=[np.int16])
    
    # convert type
    da_nos = da_nos.astype('int16')
    
    # rename
    da_nos = da_nos.rename('num_seasons')
    
    # notify user
    print('> Success!\n')
        
    return da_nos
 

def create_ds_template(da):
    """
    NOTE: this func is not yet used - requires next, unreleased ver of xarray. Will enable when ready.
    Takes an xarray data array containing veg_index values and creates an empty
    template dataset without a time dimension. This dataset is used to hold
    various phenometrics generated by phenolopy.

    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an array of veg_index
        values. 

    Returns
    -------
    ds : xarray Dataset
        A new xarray Dataset with numerous empty variables representing varioys
        phenometrics generated by phenolopy.
    """
    
    # notify user
    print('Creating a dateset template to hold various phenometrics.')
    
    # check if data array type
    if type(da) != xr.DataArray:
        raise TypeError('Not a data array. Please provide a xarray data array.')
        
    # check if time dimension is in dataset
    if 'time' not in list(da.dims):
        raise ValueError('Time dimension not in data array. Please ensure dataset has a time dimension.')
        
    # check if time dimension is in dataset
    if int(da['time'].count()) == 0:
        raise ValueError('Time dimension has no time slices. Please ensure data array contains data.')
        
    # create a single da template and remove time dimension
    temp_da = xr.full_like(da.isel(time=0), fill_value=np.nan)
    temp_da = temp_da.drop('time')
    
    # add each phenometric as a variable in template dataset
    temp_ds = xr.Dataset({
        'da_pos_values': temp_da.astype('float32'),
        'da_pos_times': temp_da.astype('int16'),
        'da_vos_values': temp_da.astype('float32'),
        'da_vos_times': temp_da.astype('int16'),
        'da_bse_values': temp_da.astype('float32'),
        'da_mos_values': temp_da.astype('float32'),        
        'aos_values': temp_da.astype('float32'),
        'da_sos_values': temp_da.astype('float32'),
        'da_sos_times': temp_da.astype('int16'),
        'da_eos_values': temp_da.astype('float32'),
        'da_eos_times': temp_da.astype('int16'),        
        'da_los_values': temp_da.astype('int16'), 
        'da_roi_values': temp_da.astype('float32'),
        'da_rod_values': temp_da.astype('float32'),
        'da_lios_values': temp_da.astype('float32'),
        'da_sios_values': temp_da.astype('float32'),
        'da_liot_values': temp_da.astype('float32'),
        'da_siot_values': temp_da.astype('float32')
    })   
    
    # notify user
    print('> Created blank dateset template successfully.\n')
    
    return temp_ds


def extract_crs(da):
    """
    Takes an xarray Dataset pulled from opendatacube and extracts crs metadata
    if exists. Returns None if not found.
    
    Parameters
    ----------
    da: xarray Dataset
        A single- or multi-dimensional array containing (or not) crs metadata.

    Returns
    -------
    crs: str
        A crs object.
    """
    
    # notify user
    print('Beginning extraction of CRS metadata.')
    try:
        # notify user
        print('> Extracting CRS metadata.')
        
        # extract crs metadata
        crs = da.geobox.crs
        
        # notify user
        print('> Success!\n')
        
    except:
        # notify user
        print('> No CRS metadata found. Returning None.\n')
        crs = None       
    
    return crs


def add_crs(ds, crs):
    """
    Takes an xarray Dataset adds previously extracted crs metadata, if exists. 
    Returns None if not found.
    
    Parameters
    ----------
    ds: xarray Dataset
        A single- or multi-dimensional array with/without crs metadata.

    Returns
    -------
    ds: xarray Dataset
        A Dataset with a new crs.
    """
    
    # notify user
    print('Beginning addition of CRS metadata.')
    try:
        # notify user
        print('> Adding CRS metadata.')
        
        # assign crs via odc utils
        ds = assign_crs(ds, str(crs))
        
        # notify user
        print('> Success!\n')
        
    except:
        # notify user
        print('> Could not add CRS metadata to data. Aborting.\n')
        pass
        
    return ds
    

def get_pos(da):
    """
    Takes an xarray DataArray containing veg_index values and calculates the vegetation 
    value and time (day of year) at peak of season (pos) for each timeseries per-pixel. 
    The peak of season is the maximum value in the timeseries, per-pixel.
    
    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values. 

    Returns
    -------
    da_pos_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at the peak of season (pos).
    da_pos_times : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) value detected at the peak of season (pos).
    """
    
    # notify user
    print('Beginning calculation of peak of season (pos) values and times.')   

    # get pos values (max val in each pixel timeseries)
    print('> Calculating peak of season (pos) values.')
    da_pos_values = da.max('time')
        
    # get pos times (day of year) at max val in each pixel timeseries)
    print('> Calculating peak of season (pos) times.')
    i = da.argmax('time', skipna=True)
    da_pos_times = da['time.dayofyear'].isel(time=i, drop=True)
    
    # convert type
    da_pos_values = da_pos_values.astype('float32')
    da_pos_times = da_pos_times.astype('int16')
    
    # rename vars
    da_pos_values = da_pos_values.rename('pos_values')
    da_pos_times = da_pos_times.rename('pos_times')

    # notify user
    print('> Success!\n')
    
    return da_pos_values, da_pos_times


def get_mos(da, da_peak_times):
    """
    Takes an xarray DataArray containing veg_index values and calculates the vegetation 
    values (time not available) at middle of season (mos) for each timeseries per-pixel. 
    The middle of season is the mean vege value and time (day of year) in the timeseries
    at 80% to left and right of the peak of season (pos) per-pixel.
    
    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values.
    da_peak_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel must be
        the time (day of year) value calculated at peak of season (pos) prior.

    Returns
    -------
    da_mos_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at the peak of season (pos).
    """
    
    # notify user
    print('Beginning calculation of middle of season (mos) values (times not possible).')  

    # get left and right slopes values
    print('> Calculating middle of season (mos) values.')
    slope_l = da.where(da['time.dayofyear'] <= da_peak_times)
    slope_r = da.where(da['time.dayofyear'] >= da_peak_times)
        
    # getupper 80% values in positive slope on left and right
    slope_l_upper = slope_l.where(slope_l >= (slope_l.max('time') * 0.8))
    slope_r_upper = slope_r.where(slope_r >= (slope_r.max('time') * 0.8))

    # get means of slope left and right
    slope_l_means = slope_l_upper.mean('time')
    slope_r_means = slope_r_upper.mean('time')

    # combine left and right veg_index means
    da_mos_values = (slope_l_means + slope_r_means) / 2
    
    # convert type
    da_mos_values = da_mos_values.astype('float32')
    
    # rename vars
    da_mos_values = da_mos_values.rename('mos_values')

    # notify user
    print('> Success!\n')
    
    #return da_mos_values
    return da_mos_values
    

def get_vos(da):
    """
    Takes an xarray DataArray containing veg_index values and calculates the vegetation 
    value and time (day of year) at valley of season (vos) for each timeseries per-pixel. 
    The valley of season is the minimum value in the timeseries, per-pixel.
    
    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values. 

    Returns
    -------
    da_vos_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at the valley of season (vos).
    da_vos_times : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) value detected at the valley of season (vos).
    """
    
    # notify user
    print('Beginning calculation of valley of season (vos) values and times.')

    # get vos values (min val in each pixel timeseries)
    print('> Calculating valley of season (vos) values.')
    da_vos_values = da.min('time')
    
    # get vos times (day of year) at min val in each pixel timeseries)
    print('> Calculating valley of season (vos) times.')
    i = da.argmin('time', skipna=True)
    da_vos_times = da['time.dayofyear'].isel(time=i, drop=True)
    
    # convert type
    da_vos_values = da_vos_values.astype('float32')
    da_vos_times = da_vos_times.astype('int16')
    
    # rename vars
    da_vos_values = da_vos_values.rename('vos_values')
    da_vos_times = da_vos_times.rename('vos_times')

    # notify user
    print('> Success!\n')
    
    return da_vos_values, da_vos_times


def get_bse(da, da_peak_times):
    """
    Takes an xarray DataArray containing veg_index values and calculates the vegetation 
    value base (bse) for each timeseries per-pixel. The base is calculated as the mean 
    value of two minimum values; the min of the slope to the left of peak of season, and
    the min of the slope to the right of the peak of season. Users must provide an existing
    peak of season (pos) data array, which can either be the max of the timeseries, or the
    middle of season (mos) values.
    
    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional DataArray containing an array of veg_index
        values.
    da_peak_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel must be
        the time (day of year) value calculated at either at peak of season (pos) or middle 
        of season (mos) prior.

    Returns
    -------
    da_bse_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        base (bse) veg_index value detected at across the timeseries at each pixel.
    """
    
    # notify user
    print('Beginning calculation of base (bse) values (times not possible).')

    # get vos values (min val in each pixel timeseries)
    print('> Calculating base (bse) values.')
    
    # split timeseries into left and right slopes via provided peak/middle values
    slope_l = da.where(da['time.dayofyear'] <= da_peak_times).min('time')
    slope_r = da.where(da['time.dayofyear'] >= da_peak_times).min('time')
    
    # get per pixel mean of both left and right slope min values 
    da_bse_values = (slope_l + slope_r) / 2
    
    # convert type
    da_bse_values = da_bse_values.astype('float32')
    
    # rename
    da_bse_values = da_bse_values.rename('bse_values')

    # notify user
    print('> Success!\n')
    
    return da_bse_values


def get_aos(da_peak_values, da_base_values):
    """
    Takes two xarray DataArrays containing the highest vege values (pos or mos) and the
    lowest vege values (bse or vos) and calculates the amplitude of season (aos). 
    The amplitude is calculated as the highest values minus the lowest values per pixel.

    Parameters
    ----------
    da_peak_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at either the peak (pos) or middle (mos) of season.
    da_base_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at either the base (bse) or valley (vos) of season.

    Returns
    -------
    da_aos_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        amplitude of season (aos) value detected between the peak and base vege values
        across the timeseries at each pixel.
    """
    
    # notify user
    print('Beginning calculation of amplitude of season (aos) values (times not possible).')

    # get aos values (peak - base in each pixel timeseries)
    print('> Calculating amplitude of season (aos) values.')
    da_aos_values = da_peak_values - da_base_values
    
    # convert type
    da_aos_values = da_aos_values.astype('float32')
    
    # rename
    da_aos_values = da_aos_values.rename('aos_values')
    
    # notify user
    print('> Success!\n')
        
    return da_aos_values


def get_sos(da, da_peak_times, da_base_values, da_aos_values, method, factor, thresh_sides, abs_value):
    """
    Takes several xarray DataArrays containing the highest vege values and times (pos or mos), 
    the lowest vege values (bse or vos), and the amplitude (aos) values and calculates the 
    vegetation values and times at the start of season (sos). Several methods can be used to
    detect the start of season; most are based on TIMESAT 3.3 methodology.

    Parameters
    ----------
    da : xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values. 
    da_peak_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) value detected at either the peak (pos) or middle (mos) of 
        season.
    da_base_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at either the base (bse) or valley (vos) of season.
    da_aos_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        amplitude of season (aos) value detected between the peak and base vege values
        across the timeseries at each pixel.
    method: str
        A string indicating which start of season detection method to use. Default is
        same as TIMESAT: seasonal_amplitude. The available options include:
        1. first_of_slope: lowest vege value of slope is sos (i.e. first lowest value).
        2. median_of_slope: middle vege value of slope is sos (i.e. median value).
        3. seasonal_amplitude: uses a percentage of the amplitude from base to find sos.
        4. absolute_value: users defined absolute value in vege index units is used to find sos.
        5. relative_amplitude: robust mean peak and base, and a factor of that area, used to find sos.
        6. stl_trendline: robust but slow - uses seasonal decomp LOESS method to find trend line and sos.
    factor: float
        A float value between 0 and 1 which is used to increase or decrease the amplitude
        threshold for the seasonal_amplitude method. A factor closer to 0 results in start 
        of season nearer to min value, a factor closer to 1 results in start of season
        closer to peak of season.
    thresh_sides: str
        A string indicating whether the sos value threshold calculation should be the min 
        value of left slope (one_sided) only, or use the bse/vos value (two_sided) calculated
        earlier. Default is two_sided, as per TIMESAT 3.3. That said, one_sided is potentially
        more robust.
    abs_value: float
        For absolute_value method only. Defines the absolute value in units of the vege index to
        which sos is defined. The part of the vege slope that the absolute value hits will be the
        sos value and time.

    Returns
    -------
    da_sos_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at the start of season (sos).
    da_sos_times : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) value detected at the start of season (sos).
    """
    
    # notify user
    print('Beginning calculation of start of season (sos) values and times.')
    
    # check factor
    if factor < 0 or factor > 1:
        raise ValueError('Provided factor value is not between 0 and 1. Aborting.')
            
    # check thresh_sides
    if thresh_sides not in ['one_sided', 'two_sided']:
        raise ValueError('Provided thresh_sides value is not one_sided or two_sided. Aborting.')
                    
    if method == 'first_of_slope':
        
        # notify user
        print('> Calculating start of season (sos) values via method: first_of_slope.')
          
        # get left slopes values, calc differentials, subset to positive differentials
        slope_l = da.where(da['time.dayofyear'] <= da_peak_times)
        slope_l_diffs = slope_l.differentiate('time')
        slope_l_pos_diffs = xr.where(slope_l_diffs > 0, True, False)
                
        # select vege values where positive on left slope
        slope_l_pos = slope_l.where(slope_l_pos_diffs)
        
        # get median of vege on pos left slope, calc vege dists from median
        slope_l_med = slope_l_pos.median('time')
        dists_from_median = slope_l_pos - slope_l_med 
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_median.isnull().all('time')
        dists_from_median = xr.where(mask, 0.0, dists_from_median)
        
        # get time index where min dist from median (first on slope)
        i = dists_from_median.argmin('time', skipna=True)
        
        # get vege start of season values and times (day of year)
        da_sos_values = slope_l_pos.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating start of season (sos) times via method: first_of_slope.')
        
        # get vege start of season times (day of year)
        da_sos_times = slope_l_pos['time.dayofyear'].isel(time=i, drop=True)

    elif method == 'median_of_slope':
        
        # notify user
        print('> Calculating start of season (sos) values via method: median_of_slope.')
          
        # get left slopes values, calc differentials, subset to positive differentials
        slope_l = da.where(da['time.dayofyear'] <= da_peak_times)
        slope_l_diffs = slope_l.differentiate('time')
        slope_l_pos_diffs = xr.where(slope_l_diffs > 0, True, False)
                
        # select vege values where positive on left slope
        slope_l_pos = slope_l.where(slope_l_pos_diffs)
        
        # get median of vege on pos left slope, calc absolute vege dists from median
        slope_l_med = slope_l_pos.median('time')
        dists_from_median = slope_l_pos - slope_l_med
        dists_from_median = abs(dists_from_median)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_median.isnull().all('time')
        dists_from_median = xr.where(mask, 0.0, dists_from_median)
        
        # get time index where min absolute dist from median (median on slope)
        i = dists_from_median.argmin('time', skipna=True)
        
        # get vege start of season values and times (day of year)
        da_sos_values = slope_l_pos.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating start of season (sos) times via method: median_of_slope.')
        
        # get vege start of season times (day of year)
        da_sos_times = slope_l_pos['time.dayofyear'].isel(time=i, drop=True)
        
    elif method == 'seasonal_amplitude':
        
        # notify user
        print('> Calculating start of season (sos) values via method: seasonal_amplitude.')
        
        # get left slopes values, calc differentials, subset to positive differentials
        slope_l = da.where(da['time.dayofyear'] <= da_peak_times)
        slope_l_diffs = slope_l.differentiate('time')
        slope_l_pos_diffs = xr.where(slope_l_diffs > 0, True, False)
                
        # select vege values where positive on left slope
        slope_l_pos = slope_l.where(slope_l_pos_diffs)
                   
        # use just the left slope min val (one), or use the bse/vos calc earlier (two) for sos
        if thresh_sides == 'one_sided':
            da_sos_values = (da_aos_values * factor) + slope_l.min('time')
        elif thresh_sides == 'two_sided':
            da_sos_values = (da_aos_values * factor) + da_base_values
                        
        # calc distance of pos vege from calculated sos value
        dists_from_sos_values = abs(slope_l_pos - da_sos_values)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_sos_values.isnull().all('time')
        dists_from_sos_values = xr.where(mask, 0.0, dists_from_sos_values)
            
        # get time index where min absolute dist from sos
        i = dists_from_sos_values.argmin('time', skipna=True)
                 
        # get vege start of season values and times (day of year)
        da_sos_values = slope_l_pos.isel(time=i, drop=True)
                
        # notify user
        print('> Calculating start of season (sos) times via method: seasonal_amplitude.')
        
        # get vege start of season times (day of year)
        da_sos_times = slope_l_pos['time.dayofyear'].isel(time=i, drop=True)
    
    elif method == 'absolute_value':
        
        # notify user
        print('> Calculating start of season (sos) values via method: absolute_value.')
        
        # get left slopes values, calc differentials, subset to positive differentials
        slope_l = da.where(da['time.dayofyear'] <= da_peak_times)
        slope_l_diffs = slope_l.differentiate('time')
        slope_l_pos_diffs = xr.where(slope_l_diffs > 0, True, False)
        
        # select vege values where positive on left slope
        slope_l_pos = slope_l.where(slope_l_pos_diffs)
        
        # calc abs distance of positive slope from absolute value
        dists_from_abs_value = abs(slope_l_pos - abs_value)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_abs_value.isnull().all('time')
        dists_from_abs_value = xr.where(mask, 0.0, dists_from_abs_value)
        
        # get time index where min absolute dist from sos (absolute value)
        i = dists_from_abs_value.argmin('time', skipna=True)
        
        # get vege start of season values and times (day of year)
        da_sos_values = slope_l_pos.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating start of season (sos) times via method: absolute_value.')
        
        # get vege start of season times (day of year)
        da_sos_times = slope_l_pos['time.dayofyear'].isel(time=i, drop=True)
        
    elif method == 'relative_value':

        # notify user
        print('> Calculating start of season (sos) values via method: relative_value.')
        
        # get left slopes values, calc differentials, subset to positive differentials
        slope_l = da.where(da['time.dayofyear'] <= da_peak_times)
        slope_l_diffs = slope_l.differentiate('time')
        slope_l_pos_diffs = xr.where(slope_l_diffs > 0, True, False)
        
        # select vege values where positive on left slope
        slope_l_pos = slope_l.where(slope_l_pos_diffs)

        # get relative amplitude via robust max and base (10% cut off either side)
        relative_amplitude = da.quantile(dim='time', q=0.90) - da.quantile(dim='time', q=0.10)
        
        # get sos value with user factor and robust mean base
        da_sos_values = (relative_amplitude * factor) + da.quantile(dim='time', q=0.10)
        
        # drop annoying quantile attribute from sos, ignore errors
        da_sos_values = da_sos_values.drop('quantile', errors='ignore')
           
        # calc abs distance of positive slope from sos values
        dists_from_sos_values = abs(slope_l_pos - da_sos_values)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_sos_values.isnull().all('time')
        dists_from_sos_values = xr.where(mask, 0.0, dists_from_sos_values)

        # get time index where min absolute dist from sos
        i = dists_from_sos_values.argmin('time', skipna=True)
                
        # get vege start of season values and times (day of year)
        da_sos_values = slope_l_pos.isel(time=i, drop=True)

        # notify user
        print('> Calculating start of season (sos) times via method: relative_value.')
        print('> Warning: this can take a long time.')
        
        # get vege start of season times (day of year)
        da_sos_times = slope_l_pos['time.dayofyear'].isel(time=i, drop=True)
        
    elif method == 'stl_trend':
        
        # notify user
        print('> Calculating end of season (eos) values via method: stl_trend.')
        
        # check if num seasons for stl is odd, +1 if not
        num_periods = len(da['time'])
        if num_periods % 2 == 0:
            num_periods = num_periods + 1
            print('> Number of stl periods is even number, added 1 to make it odd.')
        
        # prepare stl params
        stl_params = {
            'period': num_periods,
            'seasonal': 7,
            'trend': None,
            'low_pass': None,
            'robust': False
        }
        
        # prepare stl func
        def func_stl(v, period, seasonal, trend, low_pass, robust):
            return stl(v, period=period, seasonal=seasonal, trend=trend, low_pass=low_pass, robust=robust).fit().trend
        
        # notify user
        print('> Performing seasonal decomposition via LOESS. Warning: this can take a long time.')
        da_stl = xr.apply_ufunc(func_stl, da, 
                                input_core_dims=[['time']], 
                                output_core_dims=[['time']], 
                                vectorize=True, 
                                dask='parallelized', 
                                output_dtypes=[np.float32],
                                kwargs=stl_params)
        
        # notify user
        print('> Calculating start of season (sos) values via method: stl_trend.')
        
        # get left slopes values, calc differentials, subset to positive differentials
        slope_l = da.where(da['time.dayofyear'] <= da_peak_times)
        slope_l_diffs = slope_l.differentiate('time')
        slope_l_pos_diffs = xr.where(slope_l_diffs > 0, True, False)
        
        # select vege values where positive on left slope
        slope_l_pos = slope_l.where(slope_l_pos_diffs)
        
        # get min value left known pos date (todo double season issue)
        stl_l = da_stl.where(da_stl['time.dayofyear'] <= da_peak_times)
        
        # calc abs distance of positive slope from stl values
        dists_from_stl_values = abs(slope_l_pos - stl_l)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_stl_values.isnull().all('time')
        dists_from_stl_values = xr.where(mask, 0.0, dists_from_stl_values)

        # get time index where min absolute dist from sos
        i = dists_from_stl_values.argmin('time', skipna=True)
                
        # get vege start of season values and times (day of year)
        da_sos_values = slope_l_pos.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating start of season (sos) times via method: stl_trend.')
        
        # get vege start of season times (day of year)
        da_sos_times = slope_l_pos['time.dayofyear'].isel(time=i, drop=True)
    
    else:
        raise ValueError('Provided method not supported. Aborting.')
        
    # replace any all nan slices with first val and time in each pixel
    #da_sos_values = da_sos_values.where(~mask, slope_l.isel(time=0))
    #da_sos_times = da_sos_times.where(~mask, slope_l['time.dayofyear'].isel(time=0))
    da_sos_values = da_sos_values.where(~mask, np.nan)
    da_sos_times = da_sos_times.where(~mask, np.nan)
    
    # convert type
    da_sos_values = da_sos_values.astype('float32')
    da_sos_times = da_sos_times.astype('int16')
    
    # rename
    da_sos_values = da_sos_values.rename('sos_values')
    da_sos_times = da_sos_times.rename('sos_times')
            
    # notify user
    print('> Success!\n')
    
    return da_sos_values, da_sos_times


def get_eos(da, da_peak_times, da_base_values, da_aos_values, method, factor, thresh_sides, abs_value):
    """
    Takes several xarray DataArrays containing the highest vege values and times (pos or mos), 
    the lowest vege values (bse or vos), and the amplitude (aos) values and calculates the 
    vegetation values and times at the end of season (eos). Several methods can be used to
    detect the start of season; most are based on TIMESAT 3.3 methodology.

    Parameters
    ----------
    da : xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values. 
    da_peak_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) value detected at either the peak (pos) or middle (mos) of 
        season.
    da_base_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at either the base (bse) or valley (vos) of season.
    da_aos_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        amplitude of season (aos) value detected between the peak and base vege values
        across the timeseries at each pixel.
    method: str
        A string indicating which start of season detection method to use. Default is
        same as TIMESAT: seasonal_amplitude. The available options include:
        1. first_of_slope: lowest vege value of slope is eos (i.e. first lowest value).
        2. median_of_slope: middle vege value of slope is eos (i.e. median value).
        3. seasonal_amplitude: uses a percentage of the amplitude from base to find eos.
        4. absolute_value: users defined absolute value in vege index units is used to find eos.
        5. relative_amplitude: robust mean peak and base, and a factor of that area, used to find eos.
        6. stl_trendline: robust but slow - uses seasonal decomp LOESS method to find trend line and eos.
    factor: float
        A float value between 0 and 1 which is used to increase or decrease the amplitude
        threshold for the seasonal_amplitude method. A factor closer to 0 results in end 
        of season nearer to min value, a factor closer to 1 results in end of season
        closer to peak of season.
    thresh_sides: str
        A string indicating whether the sos value threshold calculation should be the min 
        value of left slope (one_sided) only, or use the bse/vos value (two_sided) calculated
        earlier. Default is two_sided, as per TIMESAT 3.3. That said, one_sided is potentially
        more robust.
    abs_value: float
        For absolute_value method only. Defines the absolute value in units of the vege index to
        which sos is defined. The part of the vege slope that the absolute value hits will be the
        sos value and time.

    Returns
    -------
    da_eos_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at the end of season (eos).
    da_eos_times : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) value detected at the end of season (eos).
    """
    
    # notify user
    print('Beginning calculation of end of season (eos) values and times.')
    
    # check factor
    if factor < 0 or factor > 1:
        raise ValueError('Provided factor value is not between 0 and 1. Aborting.')
            
    # check thresh_sides
    if thresh_sides not in ['one_sided', 'two_sided']:
        raise ValueError('Provided thresh_sides value is not one_sided or two_sided. Aborting.')
                    
    if method == 'first_of_slope':
        
        # notify user
        print('> Calculating end of season (eos) values via method: first_of_slope.')
          
        # get right slopes values, calc differentials, subset to negative differentials
        slope_r = da.where(da['time.dayofyear'] >= da_peak_times)
        slope_r_diffs = slope_r.differentiate('time')
        slope_r_neg_diffs = xr.where(slope_r_diffs < 0, True, False)
                
        # select vege values where negative on right slope
        slope_r_neg = slope_r.where(slope_r_neg_diffs)
        
        # get median of vege on neg right slope, calc vege dists from median
        slope_r_med = slope_r_neg.median('time')
        dists_from_median = slope_r_neg - slope_r_med 
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_median.isnull().all('time')
        dists_from_median = xr.where(mask, 0.0, dists_from_median)
        
        # get time index where min dist from median (first on slope)
        i = dists_from_median.argmin('time', skipna=True)
        
        # get vege end of season values and times (day of year)
        da_eos_values = slope_r_neg.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating end of season (eos) times via method: first_of_slope.')
        
        # get vege start of season times (day of year)
        da_eos_times = slope_r_neg['time.dayofyear'].isel(time=i, drop=True)

    elif method == 'median_of_slope':
        
        # notify user
        print('> Calculating end of season (eos) values via method: median_of_slope.')
          
        # get right slopes values, calc differentials, subset to positive differentials
        slope_r = da.where(da['time.dayofyear'] >= da_peak_times)
        slope_r_diffs = slope_r.differentiate('time')
        slope_r_neg_diffs = xr.where(slope_r_diffs < 0, True, False)
                
        # select vege values where negative on right slope
        slope_r_neg = slope_r.where(slope_r_neg_diffs)
        
        # get median of vege on neg right slope, calc absolute vege dists from median
        slope_r_med = slope_r_neg.median('time')
        dists_from_median = slope_r_neg - slope_r_med
        dists_from_median = abs(dists_from_median)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_median.isnull().all('time')
        dists_from_median = xr.where(mask, 0.0, dists_from_median)
        
        # get time index where min absolute dist from median (median on slope)
        i = dists_from_median.argmin('time', skipna=True)
        
        # get vege start of season values and times (day of year)
        da_eos_values = slope_r_neg.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating end of season (eos) times via method: median_of_slope.')
        
        # get vege end of season times (day of year)
        da_eos_times = slope_r_neg['time.dayofyear'].isel(time=i, drop=True)
        
    elif method == 'seasonal_amplitude':
        
        # notify user
        print('> Calculating end of season (eos) values via method: seasonal_amplitude.')
        
        # get right slopes values, calc differentials, subset to negative differentials
        slope_r = da.where(da['time.dayofyear'] >= da_peak_times)
        slope_r_diffs = slope_r.differentiate('time')
        slope_r_neg_diffs = xr.where(slope_r_diffs < 0, True, False)
        
        # select vege values where negative on right slope
        slope_r_neg = slope_r.where(slope_r_neg_diffs)       
               
        # use just the right slope min val (one), or use the bse/vos calc earlier (two) for sos
        if thresh_sides == 'one_sided':
            da_eos_values = (da_aos_values * factor) + slope_r.min('time')
        elif thresh_sides == 'two_sided':
            da_eos_values = (da_aos_values * factor) + da_base_values
            
        # calc distance of neg vege from calculated eos value
        dists_from_eos_values = abs(slope_r_neg - da_eos_values)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_eos_values.isnull().all('time')
        dists_from_eos_values = xr.where(mask, 0.0, dists_from_eos_values)
        
        # get time index where min absolute dist from eos
        i = dists_from_eos_values.argmin('time', skipna=True)
                
        # get vege end of season values and times (day of year)
        da_eos_values = slope_r_neg.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating end of season (eos) times via method: seasonal_amplitude.')
        
        # get vege end of season times (day of year)
        da_eos_times = slope_r_neg['time.dayofyear'].isel(time=i, drop=True)
    
    elif method == 'absolute_value':
        
        # notify user
        print('> Calculating end of season (eos) values via method: absolute_value.')
        
        # get right slopes values, calc differentials, subset to negative differentials
        slope_r = da.where(da['time.dayofyear'] >= da_peak_times)
        slope_r_diffs = slope_r.differentiate('time')
        slope_r_neg_diffs = xr.where(slope_r_diffs < 0, True, False)
        
        # select vege values where negative on right slope
        slope_r_neg = slope_r.where(slope_r_neg_diffs)
        
        # calc abs distance of negative slope from absolute value
        dists_from_abs_value = abs(slope_r_neg - abs_value)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_abs_value.isnull().all('time')
        dists_from_abs_value = xr.where(mask, 0.0, dists_from_abs_value)
        
        # get time index where min absolute dist from eos (absolute value)
        i = dists_from_abs_value.argmin('time', skipna=True)
        
        # get vege end of season values and times (day of year)
        da_eos_values = slope_r_neg.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating end of season (eos) times via method: absolute_value.')
        
        # get vege end of season times (day of year)
        da_eos_times = slope_r_neg['time.dayofyear'].isel(time=i, drop=True)
        
    elif method == 'relative_value':

        # notify user
        print('> Calculating end of season (eos) values via method: relative_value.')
        
        # get right slopes values, calc differentials, subset to negative differentials
        slope_r = da.where(da['time.dayofyear'] >= da_peak_times)
        slope_r_diffs = slope_r.differentiate('time')
        slope_r_neg_diffs = xr.where(slope_r_diffs < 0, True, False)
        
        # select vege values where negative on right slope
        slope_r_neg = slope_r.where(slope_r_neg_diffs)

        # get relative amplitude via robust max and base (10% cut off either side)
        relative_amplitude = da.quantile(dim='time', q=0.90) - da.quantile(dim='time', q=0.10)
        
        # get eos value with user factor and robust mean base
        da_eos_values = (relative_amplitude * factor) + da.quantile(dim='time', q=0.10)
        
        # drop annoying quantile attribute from eos, ignore errors
        da_eos_values = da_eos_values.drop('quantile', errors='ignore')
           
        # calc abs distance of negative slope from eos values
        dists_from_eos_values = abs(slope_r_neg - da_eos_values)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_eos_values.isnull().all('time')
        dists_from_eos_values = xr.where(mask, 0.0, dists_from_eos_values)

        # get time index where min absolute dist from eos
        i = dists_from_eos_values.argmin('time', skipna=True)
                
        # get vege end of season values and times (day of year)
        da_eos_values = slope_r_neg.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating end of season (eos) times via method: relative_value.')
        print('> Warning: this can take a long time.')
        
        # get vege end of season times (day of year)
        da_eos_times = slope_r_neg['time.dayofyear'].isel(time=i, drop=True)
        
    elif method == 'stl_trend':
        
        # notify user
        print('> Calculating end of season (eos) values via method: stl_trend.')
        
        # check if num seasons for stl is odd, +1 if not
        num_periods = len(da['time'])
        if num_periods % 2 == 0:
            num_periods = num_periods + 1
            print('> Number of stl periods is even number, added 1 to make it odd.')
        
        # prepare stl params
        stl_params = {
            'period': num_periods,
            'seasonal': 7,
            'trend': None,
            'low_pass': None,
            'robust': False
        }
        
        # prepare stl func
        def func_stl(v, period, seasonal, trend, low_pass, robust):
            return stl(v, period=period, seasonal=seasonal, trend=trend, low_pass=low_pass, robust=robust).fit().trend
        
        # notify user
        print('> Performing seasonal decomposition via LOESS. Warning: this can take a long time.')
        da_stl = xr.apply_ufunc(func_stl, da, 
                                input_core_dims=[['time']], 
                                output_core_dims=[['time']], 
                                vectorize=True, 
                                dask='parallelized', 
                                output_dtypes=[np.float32],
                                kwargs=stl_params)
        
        # notify user
        print('> Calculating end of season (eos) values via method: stl_trend.')
        
        # get right slopes values, calc differentials, subset to negative differentials
        slope_r = da.where(da['time.dayofyear'] >= da_peak_times)
        slope_r_diffs = slope_r.differentiate('time')
        slope_r_neg_diffs = xr.where(slope_r_diffs < 0, True, False)
        
        # select vege values where negative on right slope
        slope_r_neg = slope_r.where(slope_r_neg_diffs)
        
        # get min value right known pos date
        stl_r = da_stl.where(da_stl['time.dayofyear'] >= da_peak_times)
        
        # calc abs distance of negative slope from stl values
        dists_from_stl_values = abs(slope_r_neg - stl_r)
        
        # make mask for all nan pixels and fill with 0.0 (needs to be float)
        mask = dists_from_stl_values.isnull().all('time')
        dists_from_stl_values = xr.where(mask, 0.0, dists_from_stl_values)

        # get time index where min absolute dist from eos
        i = dists_from_stl_values.argmin('time', skipna=True)
                
        # get vege end of season values and times (day of year)
        da_eos_values = slope_r_eos.isel(time=i, drop=True)
        
        # notify user
        print('> Calculating end of season (eos) times via method: stl_trend.')
        
        # get vege end of season times (day of year)
        da_eos_times = slope_r_eos['time.dayofyear'].isel(time=i, drop=True)
        
    else:
        raise ValueError('Provided method not supported. Aborting.')

    # replace any all nan slices with last val and time in each pixel
    #da_eos_values = da_eos_values.where(~mask, slope_r.isel(time=-1))
    #da_eos_times = da_eos_times.where(~mask, slope_r['time.dayofyear'].isel(time=-1))
    da_eos_values = da_eos_values.where(~mask, np.nan)
    da_eos_times = da_eos_times.where(~mask, np.nan)
    
    # convert type
    da_eos_values = da_eos_values.astype('float32')
    da_eos_times = da_eos_times.astype('int16')
    
    # rename
    da_eos_values = da_eos_values.rename('eos_values')
    da_eos_times = da_eos_times.rename('eos_times')    
        
    # notify user
    print('> Success!\n')
    
    return da_eos_values, da_eos_times    


def get_los(da, da_sos_times, da_eos_times):
    """
    Takes two xarray DataArrays containing the start of season (sos) times (day of year) 
    and end of season (eos) times (day of year) and calculates the length of season (los). 
    This is calculated as eos day of year minus sos day of year per pixel.

    Parameters
    ----------
    da : xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values. 
    da_sos_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) detected at start of season (sos).
    da_eos_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) detected at end of season (eos).

    Returns
    -------
    da_los_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        length of season (los) value detected between the sos and eos day of year values
        across the timeseries at each pixel. The values in los represents number of days.
    """
    
    # notify user
    print('Beginning calculation of length of season (los) values (times not possible).')

    # get los values (eos day of year - sos day of year)
    print('> Calculating length of season (los) values.')
    da_los_values = da_eos_times - da_sos_times
    
    # correct los if negative values exist
    if xr.where(da_los_values < 0, True, False).any():
        
        # get max time (day of year) and negatives into data arrays
        da_max_times = da['time.dayofyear'].isel(time=-1)
        da_neg_values = da_eos_times.where(da_los_values < 0) - da_sos_times.where(da_los_values < 0)
        
        # replace negative values with max time 
        da_los_values = xr.where(da_los_values >= 0, da_los_values, da_max_times + da_neg_values)
        
        # drop time dim if exists
        da_los_values = da_los_values.drop({'time'}, errors='ignore')

    # convert type
    da_los_values = da_los_values.astype('int16')
    
    # rename
    da_los_values = da_los_values.rename('los_values')
    
    # notify user
    print('> Success!\n')
        
    return da_los_values    


def get_roi(da_peak_values, da_peak_times, da_sos_values, da_sos_times):
    """
    Takes four xarray DataArrays containing the peak season values (either pos or 
    mos) and times (day of year), and start of season (sos) values and times (day 
    of year). The rate of increase (roi) is calculated as the ratio of the difference 
    in peak and sos for vege values and time (day of year) per pixel. 

    Parameters
    ----------
    da_peak_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        vege value detected at either the peak (pos) or middle (mos) of season.
    da_peak_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) value detected at either the peak (pos) or middle (mos) of 
        season.
    da_sos_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        vege value detected at start of season (sos).
    da_sos_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) detected at start of season (sos).

    Returns
    -------
    da_roi_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        rate of increase value detected between the sos and peak values/times across the 
        timeseries at each pixel. The values in roi represents rate of vege growth.
    """
    
    # notify user
    print('Beginning calculation of rate of increase (roi) values (times not possible).')   

    # get ratio between the difference in peak and sos values and times
    print('> Calculating rate of increase (roi) values.')
    da_roi_values = (da_peak_values - da_sos_values) / (da_peak_times - da_sos_times)    

    # convert type
    da_roi_values = da_roi_values.astype('float32')
    
    # rename vars
    da_roi_values = da_roi_values.rename('roi_values')

    # notify user
    print('> Success!\n')
    
    return da_roi_values


def get_rod(da_peak_values, da_peak_times, da_eos_values, da_eos_times):
    """
    Takes four xarray DataArrays containing the peak season values (either pos or 
    mos) and times (day of year), and end of season (eos) values and times (day 
    of year). The rate of decrease (rod) is calculated as the ratio of the difference 
    in peak and eos for vege values and time (day of year) per pixel. 

    Parameters
    ----------
    da_peak_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        vege value detected at either the peak (pos) or middle (mos) of season.
    da_peak_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) value detected at either the peak (pos) or middle (mos) of 
        season.
    da_eos_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        vege value detected at end of season (eos).
    da_eos_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) detected at end of season (eos).

    Returns
    -------
    da_roi_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        rate of decrease value detected between the eos and peak values/times across the 
        timeseries at each pixel. The values in rod represents rate of vege decline.
    """
    
    # notify user
    print('Beginning calculation of rate of decrease (rod) values (times not possible).')   

    # get abs ratio between the difference in peak and eos values and times
    print('> Calculating rate of decrease (rod) values.')
    da_rod_values = abs((da_eos_values - da_peak_values) / (da_eos_times - da_peak_times))
    
    # convert type
    da_rod_values = da_rod_values.astype('float32')
    
    # rename vars
    da_rod_values = da_rod_values.rename('rod_values')

    # notify user
    print('> Success!\n')
    
    return da_rod_values


def get_lios(da, da_sos_times, da_eos_times):
    """
    Takes three xarray DataArrays containing vege values and sos/eos times (day of year) to
    calculate the long integral of season (lios) for each timeseries pixel. The lios is
    considered to act as a surrogate of vegetation productivity during growing season. The 
    long integral is calculated via the traperzoidal rule.

    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values.
    da_sos_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) detected at start of season (sos).
    da_eos_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) detected at end of season (eos).

    Returns
    -------
    da_lios_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        long integral of season (lios) value detected across the timeseries at each pixel.
    """
    
    # notify user
    print('Beginning calculation of long integral of season (lios) values (times not possible).')   

    # get vals between sos and eos times, replace any outside vals with 0
    print('> Calculating long integral of season (lios) values.')
    da_lios_values = da.where((da['time.dayofyear'] >= da_sos_times) &
                              (da['time.dayofyear'] <= da_eos_times), 0)
    
    # calculate lios using trapz (note: more sophisticated than integrate)
    da_lios_values = xr.apply_ufunc(np.trapz, da_lios_values, 
                                    input_core_dims=[['time']],
                                    dask='parallelized', 
                                    output_dtypes=[np.float32],
                                    kwargs={'dx': 1})
    
    # convert type
    da_lios_values = da_lios_values.astype('float32')
    
    # rename
    da_lios_values = da_lios_values.rename('lios_values')
    
    # notify user
    print('> Success!\n')
        
    return da_lios_values


def get_sios(da, da_sos_times, da_eos_times, da_base_values):
    """
    Takes four xarray DataArrays containing vege values,  sos and eos times and base
    values (vos or bse) and calculates the short integral of season (sios) for each 
    timeseries pixel. The sios is considered to act as a surrogate of vegetation productivity 
    minus the understorey vegetation during growing season. The short integral is calculated 
    via integrating the array with the traperzoidal rule minus the trapezoidal of the base.

    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values.
    da_sos_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) detected at start of season (sos).
    da_eos_times: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        time (day of year) detected at end of season (eos).
    da_base_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at either the base (bse) or valley (vos) of season.

    Returns
    -------
    da_sios_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        short integral of season (sios) value detected across the timeseries at each pixel.
    """
    
    # notify user
    print('Beginning calculation of short integral of season (sios) values (times not possible).')   

    # get veg vals between sos and eos times, replace any outside vals with 0
    print('> Calculating short integral of season (sios) values.')
    da_sios_values = da.where((da['time.dayofyear'] >= da_sos_times) &
                              (da['time.dayofyear'] <= da_eos_times), 0)
    
    # calculate sios using trapz (note: more sophisticated than integrate)
    da_sios_values = xr.apply_ufunc(np.trapz, da_sios_values, 
                                    input_core_dims=[['time']],
                                    dask='parallelized', 
                                    output_dtypes=[np.float32],
                                    kwargs={'dx': 1})
    
    # combine 2d base vals with 3d da, projecting const base val to pixel timeseries (i.e. a rectangle)
    da_sios_bse_values = da_base_values.combine_first(da)
    
    # get base vals between sos and eos times, replace any outside vals with 0
    da_sios_bse_values = da_sios_bse_values.where((da_sios_bse_values['time.dayofyear'] >= da_sos_times) &
                                                  (da_sios_bse_values['time.dayofyear'] <= da_eos_times), 0)
    
    # calculate trapz of base (note: more sophisticated than integrate)
    da_sios_bse_values = xr.apply_ufunc(np.trapz, da_sios_bse_values, 
                                    input_core_dims=[['time']],
                                    dask='parallelized', 
                                    output_dtypes=[np.float32],
                                    kwargs={'dx': 1})
    
    # remove base trapz from sios values
    da_sios_values = da_sios_values - da_sios_bse_values
    
    # convert type
    da_sios_values = da_sios_values.astype('float32')
    
    # rename
    da_sios_values = da_sios_values.rename('sios_values')
    
    # notify user
    print('> Success!\n')
        
    return da_sios_values


def get_liot(da):
    """
    Takes an xarray DataArray containing vege values and calculates the long integral of
    total (liot) for each timeseries pixel. The liot is considered to act as a surrogate of 
    vegetation productivty. The long integral is calculated via the traperzoidal rule.

    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values. 

    Returns
    -------
    da_liot_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        long integral of total (liot) value detected across the timeseries at each pixel.
    """
    
    # notify user
    print('Beginning calculation of long integral of total (liot) values (times not possible).')   

    # calculate liot using trapz (note: more sophisticated than integrate)
    print('> Calculating long integral of total (liot) values.')
    da_liot_values = xr.apply_ufunc(np.trapz, da, 
                                    input_core_dims=[['time']],
                                    dask='parallelized', 
                                    output_dtypes=[np.float32],
                                    kwargs={'dx': 1})
    
    # convert type
    da_liot_values = da_liot_values.astype('float32')
    
    # rename
    da_liot_values = da_liot_values.rename('liot_values')
    
    # notify user
    print('> Success!\n')
        
    return da_liot_values


def get_siot(da, da_base_values):
    """
    Takes an xarray DataArray containing vege values and calculates the short integral of
    total (siot) for each timeseries pixel. The siot is considered to act as a surrogate of 
    vegetation productivity minus the understorey vegetation. The short integral is calculated 
    via integrating the array with the traperzoidal rule minus the trapezoidal of the base.

    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values.
    da_base_values: xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        veg_index value detected at either the base (bse) or valley (vos) of season.

    Returns
    -------
    da_siot_values : xarray DataArray
        An xarray DataArray type with an x and y dimension (no time). Each pixel is the 
        short integral of total (siot) value detected across the timeseries at each pixel.
    """
    
    # notify user
    print('Beginning calculation of short integral of total (siot) values (times not possible).')   

    # calculate siot using trapz (note: more sophisticated than integrate)
    print('> Calculating short integral of total (siot) values.')
    da_siot_values = xr.apply_ufunc(np.trapz, da, 
                                    input_core_dims=[['time']],
                                    dask='parallelized', 
                                    output_dtypes=[np.float32],
                                    kwargs={'dx': 1})
    
    # combine 2d base vals with 3d da, projecting const base val to pixel timeseries (i.e. a rectangle)
    da_siot_bse_values = da_base_values.combine_first(da)
    
    # calculate trapz of base (note: more sophisticated than integrate)
    da_siot_bse_values = xr.apply_ufunc(np.trapz, da_siot_bse_values, 
                                    input_core_dims=[['time']],
                                    dask='parallelized', 
                                    output_dtypes=[np.float32],
                                    kwargs={'dx': 1})
    
    # remove base trapz from siot values
    da_siot_values = da_siot_values - da_siot_bse_values
    
    # convert type
    da_siot_values = da_siot_values.astype('float32')
    
    # rename
    da_siot_values = da_siot_values.rename('siot_values')
    
    # notify user
    print('> Success!\n')
        
    return da_siot_values
    

def calc_phenometrics(da, peak_metric='pos', base_metric='bse', method='first_of_slope', factor=0.5, thresh_sides='two_sided', abs_value=0):
    """
    Takes an xarray DataArray containing veg_index values and calculates numerous phenometrics
    (metrics that measure various aspects of plant phenoly or lifecycle). See the information
    within each metric method (e.g. get_pos) for a description of each.

    Parameters
    ----------
    da: xarray DataArray
        A two-dimensional or multi-dimensional array containing an DataArray of veg_index 
        and time values. 
    peak_metric: str
        Sets the highest value for each pixel timeseries to use for calculations that rely on 
        the highest value. Can be either pos (peak of season), which is the single highest
        value in the whole timeseries per pixel, or mos (middle of season), which is the mean 
        of the highest vales in top 90 percentile. Default is pos.
    base_metric: str
        Sets the lowest value for each pixel timeseries to use for calculations that rely on 
        the lowest value. Can be either vos (valley of season), which is the lowest possible
        value in the whole timeseries per pixel, or bse (base), which is the mean of the min
        value on the left and right slopes of the time series. Default is bse.
    method: str
        Sets the method used to determine pos (peak of season) and eos (end of season). Can be
        first_of_slope, median_of_slope, seasonal_amplitude, absolute_value, relative_amplitude,
        or stl_trend. See the get_pos and/or get_eos methods for more information.
    factor: float (>=0 and <=1)
        A float value between 0 and 1 which is used to increase or decrease the amplitude
        threshold for the get_pos and get_eos seasonal_amplitude method. A factor closer 
        to 0 results in start of season nearer to min value, a factor closer to 1 results in 
        start of season closer to peak of season.
    thresh_sides: str
        A string indicating whether the sos value threshold calculation should be the min 
        value of left slope (one_sided) only, or use the bse/vos value (two_sided) calculated
        earlier. Default is two_sided, as per TIMESAT 3.3. That said, one_sided is potentially
        more robust.
    abs_value: float
        For absolute_value method only. Defines the absolute value in units of the vege index to
        which sos is defined. The part of the vege slope that the absolute value hits will be the
        sos value and time.        
    
    Returns
    -------
    ds_phenos : xarray Dataset
        An xarray Dataset type with an x and y dimension (no time). Contains numerous
        variables representing the various phenometrics.
    """
    
    # notify user
    print('Initialising calculation of phenometrics.\n')
    
    # check if dask - not yet supported
    if dask.is_dask_collection(da):
        raise TypeError('Dask arrays not yet supported. Please compute first.')
    
    # check if dataset type
    if type(da) != xr.DataArray:
        raise TypeError('> Not a data array. Please provide a xarray data array.')
        
    # check if max metric parameters supported
    if peak_metric not in ['pos', 'mos']:
        raise ValueError('> The peak_metric parameter must be either pos or mos.')

    # check if min metric parameters supported
    if base_metric not in ['bse', 'vos']:
        raise ValueError('> The base_metric parameter must be either bse or vos.')
    
    # create template dataset to hold phenometrics
    # NOTE: disabled until next version of xarray - dask not yet enabled on all required funcs used by phenolopy
    #template = create_ds_template(da)
    
    # get crs info before work
    crs = extract_crs(da=da)
    
    # take a mask of all-nan slices for clean up at end and set all-nan to 0s
    da_all_nan_mask = da.isnull().all('time')
    da = da.where(~da_all_nan_mask, 0.0)
    
    # notify user
    print('Beginning calculation of phenometrics. This can take awhile - please wait.\n')
    
    # calc peak of season (pos) values and times
    da_pos_values, da_pos_times = get_pos(da=da)
    
    # calc valley of season (vos) values and times
    da_vos_values, da_vos_times = get_vos(da=da)
       
    # calc middle of season (mos) value (time not possible)
    da_mos_values = get_mos(da=da, da_peak_times=da_pos_times)
    
    # calc base (bse) values (time not possible).
    da_bse_values = get_bse(da=da, da_peak_times=da_pos_times)

    # calc amplitude of season (aos) values (time not possible). takes peak and base arrays
    if peak_metric == 'pos' and base_metric == 'bse':
        da_aos_values = get_aos(da_peak_values=da_pos_values, da_base_values=da_bse_values)
    elif peak_metric == 'pos' and base_metric == 'vos':
        da_aos_values = get_aos(da_peak_values=da_pos_values, da_base_values=da_vos_values)
    elif peak_metric == 'mos' and base_metric == 'bse':
        da_aos_values = get_aos(da_peak_values=da_mos_values, da_base_values=da_bse_values)
    elif peak_metric == 'mos' and base_metric == 'vos':
        da_aos_values = get_aos(da_peak_values=da_mos_values, da_base_values=da_vos_values)

    # calc start of season (sos) values and times. takes peak, base metrics and factor
    if base_metric == 'bse':
        da_sos_values, da_sos_times = get_sos(da=da, da_peak_times=da_pos_times, da_base_values=da_bse_values,
                                              da_aos_values=da_aos_values, method=method, factor=factor,
                                              thresh_sides='two_sided', abs_value=abs_value)
    elif base_metric == 'vos':
        da_sos_values, da_sos_times = get_sos(da=da, da_peak_times=da_pos_times, da_base_values=da_vos_values,
                                              da_aos_values=da_aos_values, method=method, factor=factor,
                                              thresh_sides='two_sided', abs_value=abs_value)   
    
    # calc end of season (eos) values and times. takes peak, base metrics and factor
    if base_metric == 'bse':
        da_eos_values, da_eos_times = get_eos(da=da, da_peak_times=da_pos_times, da_base_values=da_bse_values,
                                              da_aos_values=da_aos_values, method=method, factor=factor,
                                              thresh_sides='two_sided', abs_value=abs_value)
    elif base_metric == 'vos':
        da_eos_values, da_eos_times = get_eos(da=da, da_peak_times=da_pos_times, da_base_values=da_vos_values,
                                              da_aos_values=da_aos_values, method=method, factor=factor,
                                              thresh_sides='two_sided', abs_value=abs_value)
        
    
    # calc length of season (los) values (time not possible). takes sos and eos
    da_los_values = get_los(da=da, da_sos_times=da_sos_times, da_eos_times=da_eos_times)
    
    # calc rate of icnrease (roi) values (time not possible). takes peak array (pos)
    da_roi_values = get_roi(da_peak_values=da_pos_values, da_peak_times=da_pos_times,
                            da_sos_values=da_sos_values, da_sos_times=da_sos_times)
        
    # calc rate of decrease (rod) values (time not possible). takes peak array (pos)
    da_rod_values = get_rod(da_peak_values=da_pos_values, da_peak_times=da_pos_times,
                            da_eos_values=da_eos_values, da_eos_times=da_eos_times)
    
    # calc long integral of season (lios) values (time not possible)
    da_lios_values = get_lios(da=da, da_sos_times=da_sos_times, da_eos_times=da_eos_times)
    
    # calc short integral of season (sios) values (time not possible)
    if base_metric == 'bse':
        da_sios_values = get_sios(da=da, da_sos_times=da_sos_times, 
                                  da_eos_times=da_eos_times,
                                  da_base_values=da_bse_values)
    elif base_metric == 'vos':
        da_sios_values = get_sios(da=da, da_sos_times=da_sos_times, 
                                  da_eos_times=da_eos_times,
                                  da_base_values=da_vos_values)   
    
    # calc long integral of total (liot) values (time not possible)
    da_liot_values = get_liot(da=da)
    
    # calc short integral of total (siot) values (time not possible)
    if base_metric == 'bse':
        da_siot_values = get_siot(da=da, da_base_values=da_bse_values)
    elif base_metric == 'vos':
        da_siot_values = get_siot(da=da, da_base_values=da_vos_values)
        
    # create data array list
    da_list = [
        da_pos_values, 
        da_pos_times,
        da_mos_values, 
        da_vos_values, 
        da_vos_times,
        da_bse_values,
        da_aos_values,
        da_sos_values, 
        da_sos_times,
        da_eos_values, 
        da_eos_times,
        da_los_values,
        da_roi_values,
        da_rod_values,
        da_lios_values,
        da_sios_values,
        da_liot_values,
        da_siot_values
    ]
  
    # combine data arrays into one dataset
    ds_phenos = xr.merge(da_list)
    
    # set original all nan pixels back to nan
    ds_phenos = ds_phenos.where(~da_all_nan_mask)
    
    # add crs metadata back onto dataset
    ds_phenos = add_crs(ds=ds_phenos, crs=crs)
    
    # notify user
    print('Phenometrics calculated successfully!')
    
    return ds_phenos


def load_test_dataset(data_path='./data/'):
    tifs = []
    das = []

    for root, dirs, files in os.walk(data_path):
        for file in files:
            if file.endswith('.tif'):
                tifs.append(os.path.join(root, file))

    # sort year-month
    tifs.sort()
    
    for tif in tifs:
        tif_year = int(tif[-12:-8])
        tif_doy = int(tif[-7:-4])

        tif_date = datetime(tif_year, 1, 1) + timedelta(tif_doy - 1)
        tif_date = np.datetime64(tif_date, 'D')

        tif_da = xr.open_rasterio(tif)
        tif_da = tif_da.rename({'band': 'time'})
        tif_da['time'] = tif_da['time'].astype(np.datetime64)
        tif_da['time'] = [tif_date]

        das.append(tif_da)

    ds = xr.concat(das, dim='time')
    ds = ds.to_dataset(name='veg_index')
    
    return ds


def load_tifs(path):
    print('hey')