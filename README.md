# Phenolopy

## Introduction
Phenolopy (phenology + python) is a python-based library for analysing satellite timeseries data. Phenolopy has been designed to investigate the seasonality of satellite timeseries data and their relationship with dynamic vegetation properties such as phenology and temporal growth patterns. The temporal domain holds important information about short- and long-term changes within vegetation lifecycles.

Phenolopy is heavily based on the concept and methodology behind the TIMESAT software (see http://web.nateko.lu.se/timesat/timesat.asp). Unlike TIMESAT, however, Phenolopy has been designed to be a light-weight, easy to use python library that can be applied to any satellite imagery stack (e.g. MODIS, Landsat and Sentinel satellite sensors). Phenolopy also provides enhanced tools to resample timeseries data to different temporal intervals (e.g. weekly, bi-monthly or monthly spacing).

Phenolopy can be applied to derive numerous phenometrics from satellite imagery, as presented in Figure 1 below. These phenometrics can be derived for any type of land phenomena, including crop vegetation, native vegetation, and even moisture data. That said, the key target is vegetation.

## Phenometrics
Phenolopy can generate more than a dozen different metrics describing various aspects of vegetation phenology. These are explained below.
![alt text](https://github.com/lewistrotter/Phenolopy/blob/main/documentation/images/pheno_explain.png?raw=true)

The codes presented on the figure above translate to:

Code | Name | Description | Method | Value | Time
--- | --- | --- | --- | --- | ---
<sub>POS</sub> | <sub>Peak of Season</sub> | <sub>Highest vegetation value and time of season.</sub> | <sub>Maximum value in a timeseries.</sub> | :heavy_check_mark: |  :heavy_check_mark:

<sub>MOS</sub> | <sub>Middle of Season</sub> | <sub>Mean vegetation value and time of values in top 80% of season.</sub> | <sub>Mean value and time where the left and right slope edges have increased and decreased to the 80% level of the season, respectively.</sub> | :heavy_check_mark:  | :heavy_check_mark:
  
<sub>VOS | Valley of Season | Lowest vegetation value and time of season. | Minimum value in a timeseries. | X | X</sub>
<sub>BSE | Base | Mean of the lowest vegetation values in season. | Mean value of the lowest vegetation values to the left and right of Peak of Season. | X |  </sub>
<sub>SOS | Start of Season | Vegetation value and time at the start of season. | Six methods available: 1) seasonal amplitude; 2) absolute amplitude; 3) Relative amplitude; 4) LOESS STL Trend line; 5) First value of positive slope; and 6) Median value of positive slope. | X | X </sub>
<sub>EOS | End of season | Vegetation value and time at the end of season. | Six methods available: 1) seasonal amplitude; 2) absolute amplitude; 3) Relative amplitude; 4) LOESS STL Trend line; 5) First value of negative slope; and 6) Median value of negative slope. | X | X </sub>
<sub>LOS | Length of Season | Length of time (number of days) between the start and end of season. | The day of year at SOS minus EOS. |  | </sub>
<sub>ROI | Rate of Increase | The rate of vegetation "green up" at the beginning of season. | Calculated as the ratio of the difference between the left 20% and 80% levels and the corresponding time difference. |  |  </sub>
<sub>ROD | Rate of Decrease | The rate of vegetation "green down" at the end of season. | Calculated as the ratio of the difference between the right 20% and 80% levels and the corresponding time difference. |  |  </sub>
<sub>AOS | Amplitude of Season | The amplitude of vegetation values for season. | The difference between the maximum value and the VOS/BSE value. |  | 
<sub>SIOS | Short Integral of Season | Represents the seasonally active vegetation and provides a larger value for herbaceous vegetation cover and smaller value for evergreen vegetation cover. | Calculated using the trapezoidal rule on the total vegetation values between season start and end minus the VOS/BSE level value. | X |  </sub>
<sub>LIOS | Long Integral of Season | Represents the total productivity of vegetation when in season. | Calculated using the trapezoidal rule between the total vegetation values between season start and end. | X | </sub>
<sub>SIOT | Short Integral of Total | Represents total vegetation productivity throughout the season, and provides a larger value for herbaceous vegetation cover and smaller value for evergreen vegetation cover. | Calculated using the trapezoidal rule on the total vegetation values minus the VOS/BSE level value. | X | </sub>
<sub>LIOT | Long Integral of Total | Represents the total productivity of vegetation throughout the season. | Calculated using the trapezoidal rule between the total vegetation values between season start and end. | X | </sub>
<sub>NOS | Number of Seasons | Total number of seasons (i.e. prominent graph peaks) in timerseries. | Peaks detected using scipy find_peaks and any peaks are over 3 months apart. | | </sub>


## Key Technologies
- Python
- Xarray and Numpy
- OpenDataCube (ODC)
- Scipy and Statsmodel

## Demonstration
Todo.

## Sources
Todo.
