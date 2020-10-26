# Phenolopy
Phenolopy (phenology + python) is a python-based library for analysing satellite timeseries data. Phenolopy has been designed to investigate the seasonality of satellite timeseries data and their relationship with dynamic vegetation properties such as phenology and temporal growth patterns. The temporal domain holds important information about short- and long-term changes within vegetation lifecycles.

Phenolopy is heavily based on the concept and methodology behind the TIMESAT software (see http://web.nateko.lu.se/timesat/timesat.asp). Unlike TIMESAT, however, Phenolopy has been designed to be a light-weight, easy to use python library that can be applied to any satellite imagery stack (e.g. MODIS, Landsat and Sentinel satellite sensors). Phenolopy also provides enhanced tools to resample timeseries data to different temporal intervals (e.g. weekly, bi-monthly or monthly spacing).

Phenolopy can be applied to derive numerous phenometrics from satellite imagery, as presented in Figure 1 below. These phenometrics can be derived for any type of land phenomena, including crop vegetation, native vegetation, and even moisture data. That said, the key target is vegetation.

![alt text](https://github.com/lewistrotter/Phenolopy/blob/main/documentation/images/pheno_explain.png?raw=true)

The codes presented on the figure above translate to:

Code | Description | Method
--- | --- | ---
POS | Peak of season | TBA
MOS | Middle of season | TBA
NOS | Number of seasons | TBA
VOS | Valley of season | TBA
BSE | Base | TBA
SOS | Start of season | TBA
EOS | End of season | TBA
LOS | Length of season | TBA
ROI | Rate of increase | TBA
ROD | Rate of decrease | TBA
AOS | Amplitude of season | TBA
SIOS | Short integral of season | TBA
LIOS | Long integral of season | TBA
