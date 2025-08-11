# GRESO CO2 Growth Rate Analysis

**Global Rate of change of CO2 from Space-based Observations**

A Python package for calculating annual CO2 growth rates from OCO-2 satellite observations using the MIPV11 format and comparing with NOAA reference data.

## Overview

This project analyzes OCO-2 satellite data to compute global annual CO2 growth rates and validates results against NOAA ground-based measurements. The analysis uses 10-second resolution OCO-2 data to calculate spatially-weighted global averages and apply seasonal detrending to isolate long-term trends.


## Quick Start

### Prerequisites

- Python 3.7+
- Required packages: `numpy`, `pandas`, `xarray`, `matplotlib`, `scipy`, `python-dateutil`

### Installation

```bash
# Clone or download the project files
cd /path/to/project/

# Install dependencies (if not already available)
pip install numpy pandas xarray matplotlib scipy python-dateutil
```

### Running the Analysis

```bash
# Basic analysis with default settings (MIPV11 format)
python greso_analysis.py

# Analyze specific data type (OG, LNLG, or all)
python greso_analysis.py --data-type LNLG
```

### File Format Support

The analysis supports two OCO-2 data formats:

1. **MIPV11 Format (Default)**: Recommended format with improved calibration
   - Uses `xco2_2019_scale` variable
   - NetCDF4 format (`.nc4` files)
   - Download from: https://gml.noaa.gov/ccgg/OCO2_v11mip/download.php
   - **Direct file download**: https://gml.noaa.gov/aftp/user/andy/OCO-2/OCO2_b11.2_10sec_GOOD_r2.nc4

2. **Legacy Standard Format**: Original format for backward compatibility
   - Uses `xco2` variable
   - HDF5 format (`.h5` files)

```python
# Example: Using MIPV11 format (default)
from greso_analysis import AnalysisConfig, CO2GrowthRateAnalyzer

config = AnalysisConfig()
config.data_file = "/path/to/OCO2_b11.2_10sec_GOOD_r3.nc4"
# file_format is "mipv11" by default

analyzer = CO2GrowthRateAnalyzer(config)
results = analyzer.process_oco_data(["all"])

# Example: Using legacy standard format
from greso_analysis import create_standard_config

config = create_standard_config("/path/to/OCO2_b11.2-lite_10sec_GOOD_140906-250430.h5")
analyzer = CO2GrowthRateAnalyzer(config)
results = analyzer.process_oco_data(["all"])
```

## Data Requirements

### Input Data

**OCO-2 Satellite Data (MIPV11 Format - Default)**:
- **Recommended**: MIPV11 format files (e.g., `OCO2_b11.2_10sec_GOOD_r3.nc4`)
- **Download**: https://gml.noaa.gov/ccgg/OCO2_v11mip/download.php
- Format: NetCDF4 with 10-second resolution OCO-2 observations
- Variable: Uses `xco2_2019_scale` for improved calibration

**Legacy Standard Format (Optional)**:
- File: `OCO2_b11.2-lite_10sec_GOOD_140906-250430.h5`
- Format: HDF5 with 10-second resolution OCO-2 observations
- Variable: Uses `xco2`

**NOAA Reference Data**:
- Source: `ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_gr_gl.txt`
- Global CO2 growth rates from ground-based measurements
- Automatically downloaded during analysis

### Data Structure

The OCO-2 data should contain:
- `sounding_id`: Unique identifier with embedded timestamp
- `xco2`: Column-averaged CO2 concentration (ppm)
- `xco2_uncertainty`: Measurement uncertainty
- `latitude`, `longitude`: Geographic coordinates
- `data_type`: Quality flags (6=OG, <3=LNLG)
- Additional variables for filtering and analysis

## Usage Guide

### Configuration

Modify analysis parameters in the `AnalysisConfig` dataclass:

```python
config = AnalysisConfig(
    data_file="/path/to/your_data.nc4", # Path to your OCO-2 data file
    lat_range=[-50, 50],                   # Latitude bounds for analysis
    uncertainty_threshold=2.0,              # Maximum allowed uncertainty
    output_filename="annual_growth_rate.png" # Name for the output plot
)
```

### Main Components

1. **OCODataProcessor**: Handles data loading and preprocessing
2. **CO2GrowthRateAnalyzer**: Performs growth rate calculations
3. **VisualizationManager**: Creates plots and visualizations

### Output Files

- `annual_growth_rate.png`: Plot showing annual growth rates
- Console output: Tabular display of annual growth rates
- Log file: Detailed processing information

## Methodology

### Data Processing Pipeline

1. **Quality Filtering**: Apply uncertainty and data type filters
2. **Spatial Filtering**: Restrict to tropical/subtropical latitudes ([-50°, 50°])
3. **Temporal Processing**: Convert sounding IDs to datetime objects
4. **Spatial Weighting**: Area-weighted averaging by latitude
5. **Seasonal Detrending**: Remove seasonal cycles using harmonic analysis
6. **Growth Rate Calculation**: Compute annual growth rates from detrended data

### Scientific Approach

- **Spatial Coverage**: Global analysis with tropical focus
- **Temporal Resolution**: 10-second observations aggregated to monthly means
- **Detrending Method**: Harmonic analysis to remove seasonal cycles
- **Growth Rate Calculation**: Annual differences in deseasonalized CO2
- **Validation**: Direct comparison with NOAA global mean growth rates

## Results

### Visualization

The analysis generates a plot comparing OCO-2 derived annual CO2 growth rates with NOAA reference data:

![Annual CO2 Growth Rates](annual_growth_rate.png)

*Figure: Annual CO2 growth rates from OCO-2 satellite observations (blue) compared to NOAA ground-based measurements (red). The excellent agreement validates the satellite-based approach.*

### Example Output

```
============================================================
ANNUAL GROWTH RATES - ALL DATA
============================================================
Year       | Growth Rate (ppm/yr) | NOAA Growth Rate (ppm/yr)
------------------------------------------------------------
2015       |             2.9885 | 2.9500
2016       |             2.7095 | 2.8400
2017       |             2.0508 | 2.1400
2018       |             2.4630 | 2.3900
2019       |             2.7487 | 2.5000
2020       |             2.3091 | 2.3300
2021       |             2.2961 | 2.3700
2022       |             1.9195 | 2.2900
2023       |             3.0562 | 2.7300
2024       |             3.2199 | 3.7300
------------------------------------------------------------
MEAN       |             2.5046 | 2.5044
============================================================
```

### Validation

- **OCO-2 Mean Growth Rate**: 2.5050 ppm/yr
- **NOAA Mean Growth Rate**: 2.5044 ppm/yr
- **Difference**: 0.0006 ppm/yr (excellent agreement)

## File Structure

```
greso_analysis.py    # Main analysis script
utils.py            # Utility functions and classes
README.md           # This documentation file
annual_growth_rate.png  # Generated plot
```

## Development

### Code Architecture

The codebase follows object-oriented design principles:

- **AnalysisConfig**: Configuration management
- **OCODataProcessor**: Data loading and preprocessing
- **CO2GrowthRateAnalyzer**: Core analysis logic
- **VisualizationManager**: Plotting and visualization

### Key Functions

- `get_10Sec_data_baker()`: Load and filter OCO-2 MIPV11 data
- `CO2GR`: Main analysis class with binning and regional analysis
- `deseasonalize()`: Remove seasonal cycles
- `monthlyGrowthRate()`: Calculate monthly growth rates
- `areaWeighted()`: Spatial averaging with latitude weights

## Troubleshooting

### Common Issues

1. **File Not Found**: Ensure OCO-2 data file is in correct location
2. **Memory Issues**: Reduce dataset size or increase system memory
3. **Network Issues**: NOAA data download may fail - retry or use cached data

### Error Messages

- "No valid observations": Check data file path and quality filters
- "Interpolation failed": Reduce spatial/temporal filtering
- "Memory error": Process data in smaller chunks




### Testing

- Validate results against known benchmarks
- Test with different data types (OG, LNLG, all)
- Verify spatial and temporal filtering
- Check edge cases (empty datasets, missing data)

## License

This project is provided as-is for scientific research purposes. Please cite appropriately when using this code or results.

## References

**Primary Publication**:
Pandey, S., et al. (2024). Toward low-latency estimation of atmospheric CO₂ growth rates using satellite observations: Evaluating sampling errors of satellite and in situ observing approaches. *AGU Advances*, 5, e2023AV001145. https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023AV001145

## Acknowledgments

- OCO-2 Science Team for satellite data
- NOAA Global Monitoring Laboratory for reference data
- NASA/JPL for mission support and data access

## Authors

**Primary Developer**: Sudhanshu Pandey  
**Affiliation**: NASA Jet Propulsion Laboratory / California Institute of Technology  
**Contact**: sudhanshu.pandey@jpl.nasa.gov  

**Project**: GRESO (Global Rate Estimates from Satellite Observations)  
**Development Date**: August 2025  


Please cite the primary publication above when using this code or results in publications:

**Citation Format**:
```
Pandey, S., Miller, J.B., Basu, S., Liu, J., Weir, B., Byrne, B., Chevallier, F., Bowman, K.W., Liu, Z., Deng, F. and O’dell, C.W., 2024. Toward low‐latency estimation of atmospheric CO2 growth rates using satellite observations: Evaluating sampling errors of satellite and in situ observing approaches. *AGU Advances*, 5, e2023AV001145. https://doi.org/10.1029/2023AV001145
```
# stand_alone_GRESO_v2
