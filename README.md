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
POS | Peak of Season | Highest vegetation value and time of season. | Maximum value in a timeseries. | X | X
MOS | Middle of Season | Mean vegetation value and time of values in top 80% of season. | Mean value and time where the left and right slope edges have increased and decreased to the 80% level of the season, respectively. | X | X
VOS | Valley of Season | Lowest vegetation value and time of season. | Minimum value in a timeseries. | X | X
BSE | Base | Mean of the lowest vegetation values in season. | Mean value of the lowest vegetation values to the left and right of Peak of Season. | X | 

SOS | Start of Season | Vegetation value and time at the start of season. | Six methods available: 1) seasonal amplitude; 2) absolute amplitude; 3) Relative amplitude; 4) LOESS STL Trend line; 5) First value of positive slope; and 6) Median value of positive slope. | X | X
EOS | End of season | Vegetation value and time at the end of season. | Six methods available: 1) seasonal amplitude; 2) absolute amplitude; 3) Relative amplitude; 4) LOESS STL Trend line; 5) First value of negative slope; and 6) Median value of negative slope. | X | X
LOS | Length of season | TBA | TBA | TBA
ROI | Rate of increase | TBA | TBA | TBA
ROD | Rate of decrease | TBA | TBA | TBA
AOS | Amplitude of season | TBA | TBA | TBA
SIOS | Short integral of season | TBA | TBA | TBA
LIOS | Long integral of season | TBA | TBA | TBA
SIOT | Short integral of total | TBA | TBA | TBA
LIOT | Long integral of total | TBA | TBA | TBA
NOS | Number of Seasons | Total number of seasons (i.e. prominent graph peaks) in timerseries. | Peaks detected using scipy find_peaks and any peaks are over 3 months apart. | | 


## Key Technologies
- Python
- Xarray and Numpy
- OpenDataCube (ODC)
- Scipy and Statsmodel

## Demonstration
Todo.

## Sources
Todo.
