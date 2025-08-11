#!/usr/bin/env python3
"""
GRESO (Global Rate of change of CO2 from Space-based Observations) Analysis

This script analyzes OCO-2 satellite data to calculate annual CO2 growth rates
and compares them with NOAA Marine Boundary Layer reference data.

"""
import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from datetime import datetime as dt

import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import utils as ph
from utils import OCODataProcessor

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class AnalysisConfig:
    """Configuration parameters for GRESO analysis."""
    
    # Data file path (MIPV11 format by default)
    data_file: str = "/Users/pandeysu/Desktop/data_raw/OCO3_data/OCO2_10s_data/OCO2_b11.2_10sec_GOOD_r3.nc4"
    
    # File format: "standard" or "mipv11" (Baker format)
    file_format: str = "mipv11"  # "mipv11" (default) or "standard"
    
    # Analysis parameters
    lat_range: List[float] = None
    uncertainty_threshold: float = 2.0
    bin_days: int = 16
    start_time: dt = dt(2014, 10, 20)
    growth_rate_start_year: int = 2015
    growth_rate_end_year: int = 2025
    
    # Plotting parameters
    figure_size: Tuple[int, int] = (7, 4)
    figure_dpi: int = 300
    output_filename: str = "annual_growth_rate.png"
    
    # NOAA data URL
    noaa_url: str = "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_gr_gl.txt"
    
    def __post_init__(self):
        if self.lat_range is None:
            self.lat_range = [-50, 50]


def create_standard_config(data_file_path: str) -> AnalysisConfig:
    """Create configuration for legacy standard format analysis."""
    config = AnalysisConfig()
    config.data_file = data_file_path
    config.file_format = "standard"
    config.output_filename = "annual_growth_rate_standard.png"
    return config



class CO2GrowthRateAnalyzer:
    """Analyzes CO2 growth rates from OCO-2 data."""
    
    def __init__(self, config: AnalysisConfig):
        self.config = config
        self.data_processor = OCODataProcessor(config)
    
    def process_oco_data(self, data_types: List[str] = ["all"]) -> Tuple[Dict, ph.CO2GR]:
        """Process OCO-2 data and calculate regional time series."""
        oco_series = {}
        cgr = None
        
        for data_type in data_types:
            # Load and preprocess data
            raw_data = self.data_processor.load_oco_data(data_type)
            
            # Prepare data for CO2GR class
            processed_data = raw_data.copy()
            processed_data["dxco2"] = processed_data["xco2_uncertainty"]
            processed_data["model"] = processed_data["xco2"]
            processed_data["lon"] = processed_data["longitude"]
            processed_data["lat"] = processed_data["latitude"]
            
            # Initialize CO2GR processor
            cgr = ph.CO2GR()
            cgr.region = "globe"
            cgr.model = "obs"
            cgr.dat = processed_data.copy()
            
            # Bin data and calculate regional series
            cgr.binOCOdata(deldays=self.config.bin_days, 
                          st_time=self.config.start_time, 
                          aggmethod="median", 
                          giveCount=False)
            cgr.getRegionOCOSeris()
            
            # Store results
            oco_series[data_type] = cgr.obs_co2.copy()
            oco_series["time"] = cgr.obs_time.copy()
        
        return oco_series, cgr
    
    def calculate_growth_rates(self, time_series: np.ndarray, co2_data: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Calculate annual growth rates from CO2 time series."""
        # Remove NaN values
        valid_mask = ~np.isnan(co2_data)
        
        # Deseasonalize data
        deseason_time, deseason_co2 = ph.deseasonalize(time_series, co2_data, numharm=4)
        
        # Apply smoothing
        smoothed_co2 = ph.boxcar_smooth(deseason_co2, 3)
        
        # Calculate growth rates
        annual_years, annual_rates, monthly_years, monthly_rates = ph.giveGrowthRates(
            deseason_time, smoothed_co2, 
            self.config.growth_rate_start_year, 
            self.config.growth_rate_end_year
        )
        
        return annual_years, annual_rates
    
    def load_noaa_reference_data(self) -> pd.DataFrame:
        """Load NOAA Marine Boundary Layer reference data."""
        logger.info("Loading NOAA reference data")
        
        df = pd.read_csv(self.config.noaa_url, delim_whitespace=True, comment='#', 
                        header=None, names=["year", "growth", "std"])
        
        # Filter for analysis period
        valid_years = df["year"] > 2014
        return df[valid_years]


class VisualizationManager:
    """Handles plotting and visualization of results."""
    
    def __init__(self, config: AnalysisConfig):
        self.config = config
        
        # Define styling
        self.markers = {"all": "s", "LNLG": "s", "OG": "^"}
        self.colors = {"all": "black", "LNLG": "red", "OG": "green"}
    
    def create_growth_rate_plot(self, oco_results: Dict, noaa_data: pd.DataFrame, 
                               data_types: List[str] = ["all"]) -> plt.Figure:
        """Create comprehensive growth rate comparison plot."""
        fig, ax = plt.subplots(figsize=self.config.figure_size, dpi=self.config.figure_dpi)
        
        annual_growth_rates = {}
        
        for data_type in data_types:
            # Calculate growth rates
            analyzer = CO2GrowthRateAnalyzer(self.config)
            annual_years, annual_rates = analyzer.calculate_growth_rates(
                oco_results["time"], oco_results[data_type]
            )
            
            # Plot OCO-2 results
            ax.plot(annual_years - 0.5, annual_rates, 
                   marker=self.markers[data_type], 
                   color=self.colors[data_type], 
                   label=f"OCO-2 GRESO", 
                   linewidth=0.5, alpha=0.7, zorder=10)
            
            # Store results and calculate statistics
            annual_growth_rates[data_type] = annual_rates
            mean_rate = np.mean(annual_rates[:-1])
            
            # Print annual growth rates in a table
            print("\n" + "="*60)
            print(f"ANNUAL GROWTH RATES - {data_type.upper()} DATA")
            print("="*60)
            print(f"{'Year':<10} | {'OCO-2 Growth Rate (ppm/yr)':<18} | {'NOAA Growth Rate (ppm/yr)'}")
            print("-" * 60)
            
            # Print each year's growth rate
            for i, (year, rate) in enumerate(zip(annual_years, annual_rates)):
                noaa_row = noaa_data[noaa_data["year"] == int(year)]
                noaa_rate = noaa_row["growth"].values[0] if not noaa_row.empty else "N/A"
                print(f"{int(year):<10} | {rate:>18.2f} | {noaa_rate if isinstance(noaa_rate, str) else noaa_rate:.2f}")
            
            # Print mean values
            noaa_mean = np.mean(noaa_data["growth"][:-1])
            print("-" * 60)
            print(f"{'MEAN':<10} | {mean_rate:>18.2f} | {noaa_mean:.2f}")
            print("="*60 + "\n")
            
            logger.info(f"{data_type} mean growth rate: {mean_rate:.4f} ppm/yr")
            
        
        # Plot NOAA reference data
        ax.errorbar(noaa_data["year"], noaa_data["growth"], yerr=noaa_data["std"], 
                   fmt='x-', color='red', ecolor='gray', capsize=5, 
                   linewidth=0.5, label="NOAA MBL")
        
        noaa_mean = np.mean(noaa_data["growth"][:-1])
        logger.info(f"NOAA mean growth rate: {noaa_mean:.4f} ppm/yr")
        print(f"NOAA mean: {noaa_mean:.4f}")
        
        # Formatting
        ax.set_ylim(0, 4)
        ax.set_xlabel("Year")
        ax.set_ylabel("Annual Growth Rate (ppm/yr)")
        ax.set_title("Annual CO2 Growth Rate from OCO-2 and NOAA MBL")
        ax.set_ylim(1.5, 4)
        ax.set_xticks(range(2015, 2025, 1))
        ax.legend()
        ax.grid(linestyle=":", linewidth=0.5)
        
        plt.tight_layout()
        return fig, annual_growth_rates


def main():
    """Main analysis workflow."""
    logger.info("Starting GRESO CO2 growth rate analysis")
    
    # Initialize configuration
    config = AnalysisConfig()
    
    # Initialize analyzer
    analyzer = CO2GrowthRateAnalyzer(config)
    
    # Process OCO-2 data
    oco_results, cgr = analyzer.process_oco_data(data_types=["all"])
    
    # Load NOAA reference data
    noaa_data = analyzer.load_noaa_reference_data()
    
    # Create visualization
    viz_manager = VisualizationManager(config)
    fig, growth_rates = viz_manager.create_growth_rate_plot(oco_results, noaa_data)
    
    # Save results
    script_dir = Path(__file__).parent
    output_path = script_dir / config.output_filename
    fig.savefig(output_path, dpi=config.figure_dpi)
    logger.info(f"Plot saved to: {output_path}")
    
    plt.close(fig)
    
    logger.info("Analysis completed successfully")
    return growth_rates


if __name__ == "__main__":
    results = main()
