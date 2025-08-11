"""
Utility functions for CO2 Growth Rate Analysis Package
"""

import numpy as np
from datetime import datetime, timedelta
from typing import Union, List
from scipy import stats, optimize
from scipy.interpolate import interp1d
from dateutil.relativedelta import relativedelta
import pandas as pd
import matplotlib.pyplot as plt
from numpy import *


def datetime2year(dt_array: Union[np.ndarray, List]) -> Union[float, np.ndarray]:
    """
    Converts datetime objects to decimal years.
    
    Args:
        dt_array: Single datetime or array of datetime objects
        
    Returns:
        Decimal year(s)
    """
    # Handle single datetime object
    if isinstance(dt_array, datetime):
        year = dt_array.year
        start_of_year = datetime(year, 1, 1)
        if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
            days_in_year = 366
        else:
            days_in_year = 365
        day_of_year = (dt_array - start_of_year).days + (dt_array - start_of_year).seconds / 86400.0
        return year + day_of_year / days_in_year
    
    # Handle array of datetime objects
    dt_array = np.asarray(dt_array)
    if dt_array.size == 0:
        return np.array([])
    
    # Convert each datetime to decimal year
    result = []
    for dt_obj in dt_array.flat:
        if isinstance(dt_obj, datetime):
            year = dt_obj.year
            start_of_year = datetime(year, 1, 1)
            if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
                days_in_year = 366
            else:
                days_in_year = 365
            day_of_year = (dt_obj - start_of_year).days + (dt_obj - start_of_year).seconds / 86400.0
            result.append(year + day_of_year / days_in_year)
        else:
            # Handle numpy datetime64 or other formats
            dt_obj = pd.to_datetime(dt_obj).to_pydatetime()
            year = dt_obj.year
            start_of_year = datetime(year, 1, 1)
            if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
                days_in_year = 366
            else:
                days_in_year = 365
            day_of_year = (dt_obj - start_of_year).days + (dt_obj - start_of_year).seconds / 86400.0
            result.append(year + day_of_year / days_in_year)
    
    return np.array(result).reshape(dt_array.shape)


def year2datetime(year_float: Union[float, np.ndarray]) -> Union[datetime, np.ndarray]:
    """
    Converts decimal year(s) to datetime object(s).
    
    Args:
        year_float: Decimal year or array of decimal years
        
    Returns:
        Datetime object(s)
    """
    if np.isscalar(year_float):
        year = int(year_float)
        fraction = year_float - year
        start_of_year = datetime(year, 1, 1)
        if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
            days_in_year = 366
        else:
            days_in_year = 365
        return start_of_year + timedelta(days=fraction * days_in_year)
    else:
        return np.array([year2datetime(y) for y in year_float])


def give_grid_cell_area(lons_edges: np.ndarray, lats_edges: np.ndarray) -> np.ndarray:
    """
    Calculates the area of each cell in a regular lat/lon grid.
    
    Args:
        lons_edges: Longitude edges of grid cells
        lats_edges: Latitude edges of grid cells
        
    Returns:
        2D array of cell areas in km²
    """
    R_earth = 6371.0  # Earth radius in km
    
    # Convert to radians
    lons_rad = np.deg2rad(lons_edges)
    lats_rad = np.deg2rad(lats_edges)
    
    # Calculate differences
    dlon = np.diff(lons_rad)
    dlat = np.diff(lats_rad)
    
    # Create meshgrids
    dlon_grid, dlat_grid = np.meshgrid(dlon, dlat)
    _, lats_center_grid = np.meshgrid(
        (lons_edges[:-1] + lons_edges[1:]) / 2,
        (lats_edges[:-1] + lats_edges[1:]) / 2
    )
    
    # Calculate area
    area = R_earth**2 * dlon_grid * np.sin(np.deg2rad(lats_center_grid + dlat_grid/2)) - \
           R_earth**2 * dlon_grid * np.sin(np.deg2rad(lats_center_grid - dlat_grid/2))
    
    return np.abs(area)


def area_weighted_mean(glon_centers: np.ndarray, glat_centers: np.ndarray, 
                      data_grid: np.ndarray) -> float:
    """
    Calculates the area-weighted mean of a 2D data grid.
    
    Args:
        glon_centers: Longitude centers of grid cells
        glat_centers: Latitude centers of grid cells  
        data_grid: 2D data array
        
    Returns:
        Area-weighted mean value
    """
    # Create edges from centers
    dlon = glon_centers[1] - glon_centers[0] if len(glon_centers) > 1 else 1.0
    dlat = glat_centers[1] - glat_centers[0] if len(glat_centers) > 1 else 1.0
    
    lons_edges = np.concatenate([glon_centers - dlon/2, [glon_centers[-1] + dlon/2]])
    lats_edges = np.concatenate([glat_centers - dlat/2, [glat_centers[-1] + dlat/2]])
    
    # Calculate cell areas
    areas = give_grid_cell_area(lons_edges, lats_edges)
    
    # Mask invalid data
    valid_mask = ~np.isnan(data_grid)
    
    if not np.any(valid_mask):
        return np.nan
    
    # Calculate weighted mean
    weighted_sum = np.sum(data_grid[valid_mask] * areas[valid_mask])
    total_area = np.sum(areas[valid_mask])
    
    return weighted_sum / total_area if total_area > 0 else np.nan


def give_month_edges(decimal_years: np.ndarray) -> np.ndarray:
    """
    Given decimal years, return an array of decimal years at the start of each month.
    
    Args:
        decimal_years: Array of decimal years
        
    Returns:
        Array of month edge decimal years
    """
    if len(decimal_years) == 0:
        return np.array([])
    
    start_year = int(np.floor(np.min(decimal_years)))
    end_year = int(np.ceil(np.max(decimal_years))) + 1
    
    month_edges = []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            month_edges.append(year + (month - 1) / 12.0)
    
    return np.array(month_edges)


def boxcar_smooth(obs, jx=1):
    """
    Apply a boxcar (moving average) smoothing to the input array.
    
    Args:
        obs: Input array to be smoothed
        jx: Half-width of the smoothing window
        
    Returns:
        Smoothed array
    """
    nobs = obs.copy()
    for ii in range(jx, len(obs)-jx):
        nobs[ii] = np.mean(obs[ii-jx: ii+jx+1])
    return nobs

def getNOAAgrowthRates():
    url = "ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_gr_gl.txt"
    df = pd.read_csv(url, delim_whitespace=True, comment='#', header=None, names=["year", "grwoth", "std"])
    vd= df["year"] > 2014;    df= df[vd]
    #plt.plot(df["year"], df["grwoth"], marker='x', linestyle='-', label="NOAA MBL", color = "black", lw = 0.5)
    errorbar(df["year"], df["grwoth"], yerr=df["std"], fmt='x-', color='black', ecolor='gray', capsize=5, lw = 0.5, label="NOAA MBL")
    print ("NOAA mean", mean(df["grwoth"][:-1]))


# Functions from functions_OCO_GR.py needed for the main script

def error_mean(x):
    return np.std(x) / np.sqrt(len(x))

def meann(x):
    return np.nanmean(x)

def harmonics(params, x, numpoly, numharm):
    pi2 = 2*np.pi*x
    if numharm > 0:
        # create an array s of correct size by explicitly evaluating first harmonic
        s = params[numpoly]*np.sin(pi2) + params[numpoly+1]*np.cos(pi2)
        # do additional harmonics (nharm > 1)
        for i in range(1,numharm):
            ix = 2*i + numpoly    # index into params for harmonic coefficients
            s += params[ix]*np.sin((i+1)*pi2) + params[ix+1]*np.cos((i+1)*pi2)

        n = numpoly + 2*numharm
        if n < len(params):
            s = (1+params[n]*x) * s        # amplitude gain factor

        return s
    else:
        return 0

def fitFunc(params, x, numpoly, numharm):
    p = np.polyval(params[numpoly-1::-1], x)
    s = harmonics(params, x, numpoly, numharm)
    return p+s

def errfunc(p, x, y, numpoly, numharm):
    return y - fitFunc(p, x, numpoly, numharm)

def deseasonalize(obs_time, obs, numpoly=4, numharm=4, verbose=False):
    vd = np.where(~np.isnan(obs))
    pm = [1.0] * (numpoly + 2*numharm)
    params, covar, info, mesg, ier = optimize.leastsq(errfunc, pm, full_output=1, args=(obs_time[vd], obs[vd], numpoly, numharm))
    if verbose: 
        print("fit quality", info['qtf'][0])
    xi = obs_time[vd]
    yi = obs[vd]
    fi = fitFunc(params, xi, numpoly, numharm)
    return xi, np.polyval(params[numpoly-1::-1], xi) + yi - fi

def areaWeighted(glon, glat, aa):
    if len(np.where(~np.isnan(aa))[0]) == 0: 
        return np.nan
    
    # Handle 1D case (single latitude band)
    if aa.ndim == 1:
        # Simple area-weighted mean for 1D case
        area_weights = np.cos(np.deg2rad(glat))
        valid_mask = ~np.isnan(aa)
        if not np.any(valid_mask):
            return np.nan
        weighted_sum = np.sum(aa[valid_mask] * area_weights[valid_mask])
        total_weight = np.sum(area_weights[valid_mask])
        return weighted_sum / total_weight if total_weight > 0 else np.nan
    
    # Handle 2D case
    vvd = np.isnan(aa).all(axis=-1)
    bcc = np.zeros((aa.shape[0]))
    bcc[:] = np.nan
    bcc[~vvd] = np.nanmean(aa[~vvd, :], axis=1)
    aa = bcc
    
    # Create area weights based on latitude
    area_weights = np.cos(np.deg2rad(glat))
    valid_mask = ~np.isnan(aa)
    if not np.any(valid_mask):
        return np.nan
    
    weighted_sum = np.sum(aa[valid_mask] * area_weights[valid_mask])
    total_weight = np.sum(area_weights[valid_mask])
    return weighted_sum / total_weight if total_weight > 0 else np.nan

def giverMonthEdges(xxi):
    start_year = int(np.floor(np.min(xxi)))
    end_year = int(np.ceil(np.max(xxi))) + 1
    month_edges = []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            month_edges.append(year + (month - 1) / 12.0)
    return np.array(month_edges)

def monthlyGrowthRate(obs_time, mgrarea, verbose=False):
    vd = np.where(~np.isnan(mgrarea))
    deobs_time, deobs = obs_time[vd], mgrarea[vd]
    
    if len(deobs_time) < 2:
        return np.array([]), np.array([])
    
    f2 = interp1d(deobs_time, deobs, bounds_error=False, fill_value='extrapolate')
    obs_gr = []
    xx = []
    monthexx = giverMonthEdges(obs_time[vd])
    
    # Filter month edges to be within the data range
    min_time, max_time = np.min(deobs_time), np.max(deobs_time)
    valid_months = (monthexx >= min_time) & (monthexx <= max_time)
    monthexx = monthexx[valid_months]
    
    for ii, mt in enumerate(monthexx[:-1]):
        if monthexx[ii+1] <= max_time and monthexx[ii] >= min_time:
            obs_gr.append(f2(monthexx[ii+1]) - f2(monthexx[ii]))
            xx.append((monthexx[ii+1] + monthexx[ii])/2)
    return np.array(xx), np.array(obs_gr)

def giveAnnual(y):
    aaa = []
    for ii in np.arange(0, len(y)-12, 12):
        aaa.append(np.sum(y[ii:ii+12]))
    aaa.append(np.sum(y[len(y)-12:len(y)]))
    return np.array(aaa)

def giveGrowthRates(oxc, oyc, st_tim=2015, ed_tim=2021.0, verbose=False):
    oxc, oyc = deseasonalize(oxc, oyc)
    oyc = boxcar_smooth(oyc, 1)  # Simple smoothing instead of wavelet
    oxg, oyg = monthlyGrowthRate(oxc, oyc)
    vvd = np.logical_and(oxg > st_tim, oxg < ed_tim)
    oxg, oyg = oxg[vvd], oyg[vvd]
    annuy = giveAnnual(oyg)
    annux = giveAnnual(oxg) / 12
    return annux, annuy, oxg, oyg

class CO2GR:
    def __init__(self):
        self.inv_runs = ["IS", "LNLG", "LNLGIS", "OG", "LNLGOGIS"]
        self.regions = ['globe', "NET", "TRO", 'SET', "NH", "SH"]
        self.models = ["TM5-4DVAR", "OU", "Baker", "COLA", "AMES", "UT", "WOMBAT", "CAMS"]
        self.pgc2ppm = 2.124
        self.region = "globe"
        self.inv_run = "IS"
        self.flux_name = "TM5-4DVAR"
        self.model_name = "TM5-4DVAR"
        self.ttim_once = None
        self.miss_month = None
        self.interp_kind = 'linear'
        self.binres = 5
        self.rand_err = None

    def binOCOdata(self, deldays=16, st_time=None, end_time=None, giveCount=False, doMonth=False, verbose=True, aggmethod="mean"):
        dat = self.dat
        tim1 = st_time
        
        if self.binres == 20:
            glat = np.arange(-80, 90, 5)
            glon = np.arange(-170, 180, 5)
        elif self.binres == 10:
            glat = np.arange(-85, 90, 10)
            glon = np.arange(-175, 180, 10)
        elif self.binres == 5:
            glat = np.arange(-87.5, 90, 5)
            glon = np.arange(-177.5, 180, 5)
        
        mid_time = []
        if end_time is None:
            end_time = dat["tim"][-1]
            
        while tim1 < end_time:
            if doMonth:
                tim2 = tim1 + relativedelta(months=1)
            else:
                tim2 = tim1 + relativedelta(days=deldays)
            
            mid_time.append(datetime2year(tim1 + timedelta(days=(tim2-tim1).days/2)))
            vd = np.where((tim1 <= dat["tim"]) & (dat["tim"] < tim2))
            x, y, z = dat["lon"][vd], dat["lat"][vd], dat["model"][vd]
            
            bi_val = stats.binned_statistic_2d(x, y, z, bins=[len(glon), len(glat)], 
                                             range=[[-180, 180], [-90, 90]], statistic=aggmethod)
            temp = bi_val[0].transpose()[np.newaxis]
            
            errtemp = np.zeros_like(temp)
            ccount = np.zeros_like(temp)
            
            if giveCount:
                bi_val = stats.binned_statistic_2d(x, y, dat["dxco2"][vd], bins=[len(glon), len(glat)], 
                                                 range=[[-180, 180], [-90, 90]], statistic=error_mean)
                errtemp = bi_val[0].transpose()[np.newaxis]
                bi_val = stats.binned_statistic_2d(x, y, dat["dxco2"][vd], bins=[len(glon), len(glat)], 
                                                 range=[[-180, 180], [-90, 90]], statistic="count")
                ccount = bi_val[0].transpose()[np.newaxis]
            
            if tim1 == st_time:
                wcube = temp
                stderr = errtemp
                count = ccount
            else:
                wcube = np.append(wcube, temp, axis=0)
                stderr = np.append(stderr, errtemp, axis=0)
                count = np.append(count, ccount, axis=0)
            
            tim1 = tim2
        
        mid_time = np.array(mid_time)
        self.gridobs = {'lon': glon, 'lat': glat, 'model_xco2': wcube, 'dxco2': stderr, "time": mid_time, "count": count}

    def getRegionOCOSeris(self):
        region = self.region
        
        if self.binres == 20:
            latvd = np.arange(9)
        elif self.binres == 10:
            latvd = np.arange(18)
            if region == "NH": latvd = np.arange(9, 18)
            if region == "SH": latvd = np.arange(9)
            if region == "NET": latvd = np.arange(12, 18)
            if region == "SET": latvd = np.arange(6)
            if region == "TRO": latvd = np.arange(6, 12)
            if region == "NTRO": latvd = np.arange(9, 12)
            if region == "STRO": latvd = np.arange(6, 9)
        elif self.binres == 5:
            latvd = np.arange(36)
            if region == "NH": latvd = np.arange(9*2, 18*2)
            if region == "SH": latvd = np.arange(9*2)
            if region == "NET": latvd = np.arange(23, 36)
            if region == "SET": latvd = np.arange(13)
            if region == "TRO": latvd = np.arange(6*2, 12*2)
            if region == "NTRO": latvd = np.arange(18, 23)
            if region == "STRO": latvd = np.arange(13, 18)
            if region == "NH20+": latvd = np.arange(int(110/5), int(180/5))
            if region == "SH20-": latvd = np.arange(0, int(70/5))
            if region == "TRO20": latvd = np.arange(int(70/5), int(110/5))
        
        self.latvd = latvd
        self.miss_month = 2017.6
        
        # observed growth rate calculation
        glon, glat, wcube, errcube = self.gridobs["lon"], self.gridobs["lat"], self.gridobs["model_xco2"], self.gridobs['dxco2']
        wcube[errcube > 2] = np.nan
        mid_time = self.gridobs["time"]
        
        mgrarea = []
        for im in np.arange(wcube.shape[0]):
            if len(np.where(~np.isnan(wcube[im][self.latvd]))[0]) == 0:
                mgrarea.append(np.nan)
                continue
            ssss = areaWeighted(glon, glat[self.latvd], wcube[im][self.latvd])
            mgrarea.append(ssss)
        
        self.obs_time = np.array(mid_time)
        self.obs_co2 = np.array(mgrarea)


def get_10Sec_data_baker(file_path_10s, data_type="all", lat_range=[-50, 50]):
    """
    Load OCO-2 data from MIPV11 format files (Baker format).
    
    Parameters:
    -----------
    file_path_10s : str
        Path to the OCO-2 data file in MIPV11 format
    data_type : str
        Type of data to load ("OG", "LNLG", or "all")
    lat_range : list
        Latitude range for filtering [min_lat, max_lat]
    
    Returns:
    --------
    dict
        Dictionary containing filtered OCO-2 data
    """
    addat = {}
    
    if True:
        import xarray as xr
        fid = xr.open_dataset(file_path_10s)
        
        # Apply data type filters
        if data_type == "OG":
            vd1 = np.logical_and(fid["xco2_uncertainty"] < 2, fid["data_type"] == 6)
        elif data_type == "LNLG":
            vd1 = np.logical_and(fid["xco2_uncertainty"] < 2, fid["data_type"] < 3)
        else:  # "all"
            vd1 = np.logical_and(fid["xco2_uncertainty"] < 2, fid["data_type"] != 3)
        
        # Apply latitude filter
        vd2 = np.logical_and(fid.latitude > lat_range[0], fid.latitude < lat_range[1])
        vd = np.logical_and(vd1, vd2)
        
        # Extract variables (note: using xco2_2019_scale instead of xco2)
        reab = {}
        variables = ["sounding_id", "xco2_2019_scale", "xco2_apriori", "xco2_uncertainty", 
                    "latitude", "operation_mode", "longitude", "psurf_apriori", 
                    "xco2_averaging_kernel", "pressure_weight"]
        
        for kk in variables:
            reab[kk] = fid[kk].values[vd]
        
        # Convert sounding_id to datetime using the slower method (as in original)
        from datetime import datetime as dt
        reab['tim'] = np.array([dt.strptime(str(tt).split('.')[0], "%Y%m%d%H%M%S%f") 
                               for tt in reab["sounding_id"]])
        
        # Convert datetime to decimal year
        reab["tim_yf"] = datetime2year(reab["tim"])
        
        # Map xco2_2019_scale to xco2 for compatibility
        reab["xco2"] = reab["xco2_2019_scale"]
    
    return reab


class OCODataProcessor:
    """Handles loading and preprocessing of OCO-2 satellite data."""
    
    def __init__(self, config):
        """Initialize with AnalysisConfig object."""
        self.config = config
        
    def _extract_datetime_components(self, sounding_ids: np.ndarray) -> dict:
        """Extract datetime components from sounding_id using optimized integer arithmetic."""
        sounding_ids = sounding_ids.astype(np.int64)
        
        # Extract components using integer division and modulo
        years = sounding_ids // 10000000000
        remainder = sounding_ids % 10000000000
        months = remainder // 100000000
        remainder = remainder % 100000000
        days = remainder // 1000000
        remainder = remainder % 1000000
        hours = remainder // 10000
        remainder = remainder % 10000
        minutes = remainder // 100
        seconds = remainder % 100
        
        return {
            'years': years, 'months': months, 'days': days,
            'hours': hours, 'minutes': minutes, 'seconds': seconds
        }
    
    def _convert_to_decimal_years(self, datetime_components: dict) -> np.ndarray:
        """Convert datetime components to decimal years efficiently."""
        years = datetime_components['years']
        months = datetime_components['months']
        days = datetime_components['days']
        hours = datetime_components['hours']
        minutes = datetime_components['minutes']
        seconds = datetime_components['seconds']
        
        # Day of year calculation
        days_before_month = np.array([0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334])
        is_leap = ((years % 4 == 0) & (years % 100 != 0)) | (years % 400 == 0)
        
        day_of_year = days_before_month[months-1].astype(float) + days.astype(float)
        day_of_year += (is_leap & (months > 2)).astype(float)  # Add leap day
        
        # Add fractional day
        day_of_year += (hours.astype(float) + minutes.astype(float)/60.0 + 
                       seconds.astype(float)/3600.0) / 24.0
        
        # Calculate decimal year
        days_in_year = 365.0 + is_leap.astype(float)
        return years.astype(float) + (day_of_year - 1.0) / days_in_year
    
    def _create_datetime_objects(self, sounding_ids: np.ndarray) -> np.ndarray:
        """Create datetime objects efficiently using pandas."""
        import pandas as pd
        sounding_str = sounding_ids.astype(str).astype('U14')
        return pd.to_datetime(sounding_str, format='%Y%m%d%H%M%S', errors='coerce').to_pydatetime()
    
    def _apply_quality_filters(self, dataset, data_type: str) -> np.ndarray:
        """Apply quality filters based on data type and uncertainty."""
        # Data type filter
        if data_type == "OG":
            type_filter = np.logical_and(dataset["xco2_uncertainty"] < self.config.uncertainty_threshold, 
                                       dataset["data_type"] == 6)
        elif data_type == "LNLG":
            type_filter = np.logical_and(dataset["xco2_uncertainty"] < self.config.uncertainty_threshold, 
                                       dataset["data_type"] < 3)
        else:  # "all"
            type_filter = np.logical_and(dataset["xco2_uncertainty"] < self.config.uncertainty_threshold, 
                                       dataset["data_type"] != 3)
        
        # Latitude filter
        lat_filter = np.logical_and(dataset.latitude > self.config.lat_range[0], 
                                  dataset.latitude < self.config.lat_range[1])
        
        return np.logical_and(type_filter, lat_filter)
    
    def _load_standard_data(self, data_type: str) -> dict:
        """Load data using legacy standard format."""
        import xarray as xr
        dataset = xr.open_dataset(self.config.data_file)
        
        # Apply quality filters
        valid_mask = self._apply_quality_filters(dataset, data_type)
        
        # Extract filtered data
        data = {
            'sounding_id': dataset.sounding_id.values[valid_mask],
            'xco2': dataset.xco2.values[valid_mask],
            'xco2_apriori': dataset.xco2_apriori.values[valid_mask],
            'xco2_uncertainty': dataset.xco2_uncertainty.values[valid_mask],
            'latitude': dataset.latitude.values[valid_mask],
            'longitude': dataset.longitude.values[valid_mask],
            'psurf_apriori': dataset.psurf_apriori.values[valid_mask],
            'xco2_averaging_kernel': dataset.xco2_averaging_kernel.values[valid_mask],
            'pressure_weight': dataset.pressure_weight.values[valid_mask],
            'data_type': dataset.data_type.values[valid_mask]
        }
        
        # Convert sounding_id to datetime
        data['tim'] = self._create_datetime_objects(data['sounding_id'])
        data['tim_yf'] = self._convert_to_decimal_years(
            self._extract_datetime_components(data['sounding_id']))
        
        return data
    
    def _load_mipv11_data(self, data_type: str) -> dict:
        """Load data using MIPV11 format (default)."""
        try:
            return get_10Sec_data_baker(self.config.data_file, data_type, self.config.lat_range)
        except Exception as e:
            import logging
            logger = logging.getLogger(__name__)
            logger.error(f"Failed to load MIPV11 format data: {e}")
            logger.info("OCO-2 MIPV11 files can be downloaded from:")
            logger.info("https://gml.noaa.gov/ccgg/OCO2_v11mip/download.php")
            raise
    
    def load_oco_data(self, data_type: str = "all") -> dict:
        """Load OCO-2 data based on configured format."""
        import logging
        logger = logging.getLogger(__name__)
        
        logger.info(f"Loading OCO-2 data for type: {data_type}")
        logger.info(f"Using file format: {self.config.file_format}")
        
        if self.config.file_format == "standard":
            data = self._load_standard_data(data_type)
        else:  # mipv11 format (default)
            data = self._load_mipv11_data(data_type)
            
        logger.info(f"Loaded {len(data['xco2'])} valid observations")
        return data