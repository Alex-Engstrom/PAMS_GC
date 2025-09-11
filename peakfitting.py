# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 16:06:35 2025

@author: aengstrom
"""
import chromclass as cc
from pybaselines import Baseline
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import math
from scipy.signal import savgol_filter, find_peaks
from scipy.optimize import curve_fit
from scipy.integrate import simpson
from scipy.stats import chi2

#%%
def baseline(data, half_window=20, threshold=0.003, min_length=5, min_fwhm = None, smooth_half_window=5):
    """
    Baseline correction using pybaselines fastchrom
    Returns numpy arrays with retention time and data
    """
    x = data[0]  # Retention times
    y = data[1]  # Signal values
    
    try:
        
        # Initialize baseline fitter
        baseline_fitter = Baseline(x)
        
        # Apply fastchrom method
        baseline_vals, params = baseline_fitter.fastchrom(
            y,
            half_window=half_window,
            threshold=threshold,
            min_length=min_length,
            min_fwhm = min_fwhm,
            smooth_half_window=smooth_half_window
        )
        
        # The mask indicates baseline regions (True = baseline, False = peaks)
        baseline_mask = params['mask']
        
        # Invert to get peak regions (True = peaks, False = baseline)
        peak_mask = ~baseline_mask
        
        # Baseline corrected data
        corrected_y = y - baseline_vals
        
        # Calculate noise standard deviation*3 as approximate LOD
        baseline_noise = np.std(corrected_y[baseline_mask])
        noise_threshold = baseline_noise * 3
        
        # Create combined numpy arrays with retention time and data
        original_chrom = np.vstack((x, y))  # Time (sec) + original signal
        baseline_arr = np.vstack((x, baseline_vals))  # Time + baseline
        corrected_arr = np.vstack((x, corrected_y))  # Time + corrected signal
        full_data_dict = {
            'original_chrom': original_chrom,
            'baseline_arr': baseline_arr,
            'corrected_arr': corrected_arr,
            'baseline_mask': baseline_mask,
            'peak_mask': peak_mask,
            'baseline_noise': baseline_noise,  # 1-sigma noise
            'noise_threshold': noise_threshold,  # 3-sigma threshold
            'params': params
        }
        return full_data_dict
        
    except Exception as e:
        print(f"fastchrom baseline failed: {e}")
        # Fallback approach - return arrays with zero baseline
        zero_baseline = np.zeros_like(data[1])
        
        full_data_dict = {
            'original_chrom': np.vstack((x, data[1])),
            'baseline_arr': np.vstack((x, zero_baseline)),
            'corrected_arr': np.vstack((x, data[1])),
            'baseline_mask': np.ones_like(data[1], dtype=bool),
            'peak_mask': np.zeros_like(data[1], dtype=bool),
            'baseline_noise': 0.001,
            'noise_threshold': 0.003,
            'params': {}
        }
        return full_data_dict
from typing import List, Tuple



def find_peak_regions(mask_1d: np.ndarray) -> List[Tuple[int, int]]:
    """
    Find start and end indices of continuous True regions in a 1D boolean array.
    
    Args:
        mask_1d: 1D boolean array
        
    Returns:
        List of tuples (start_index, end_index) for each continuous True region
    """
    if mask_1d.size == 0:
        return []
    # Find where the mask changes from False to True or True to False
    changes = np.where(mask_1d[:-1] != mask_1d[1:])[0] + 1
    # If the array starts with True, add start at 0
    starts = []
    ends = []
    
    if mask_1d[0]:
        starts.append(0)
    
    # Add all changes from False→True as starts
    starts.extend(changes[1::2] if mask_1d[0] else changes[::2])
    # Add all changes from True→False as ends
    ends.extend(changes[::2] if mask_1d[0] else changes[1::2])
    
    # If the array ends with True, add end at the last index
    if mask_1d[-1]:
        ends.append(len(mask_1d))
    return list(zip(starts, ends))

def extract_valid_peak_regions_to_dataframe(data_2d: np.ndarray, mask_1d: np.ndarray, limit = None) -> pd.DataFrame:
    """
 
    """

    regions = find_peak_regions(mask_1d)
    # Pre-allocate data structure
    data_dict = {}
    
    for i, (start, end) in enumerate(regions):
        region_data = data_2d[:, start:end]
        if abs(np.mean(region_data[1,:])) > limit:

            region_name = f"region {i} {data_2d[0, start]:.3f} to {data_2d[0, end]:.3f}"
            
            data_dict[region_name] = [region_data]
        else:
            continue
    
    # Single DataFrame creation
    return pd.DataFrame.from_dict(data_dict)
def gaussian(x, amplitude, center, sigma):
    """Gaussian function for peak fitting"""
    return amplitude * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))

def identify_peak_centers(peak_region_data, noise_threshold):
    """
    Identify approximate peak centers within a peak region using scipy's find_peaks.
    
    Args:
        peak_region_data: 2D numpy array with [time, signal] for the peak region
        noise_threshold: Minimum height threshold for peak detection
        
    Returns:
        List of approximate peak center indices and their properties
    """
    time = peak_region_data[0]
    signal = peak_region_data[1]
    # Find peaks using scipy's find_peaks
    try:
        peaks, properties = find_peaks(
            signal, 
            height=noise_threshold, 
            prominence=noise_threshold,
            width=5,
            distance=.5
        )
    except Exception as e:
        print(f"Peak finding failed: {e}")
        return []  # Return empty list on failure
    
    peak_centers = []
    for i, peak_idx in enumerate(peaks):
        width_overrange = False
        width_max, width_min = time[peak_idx] + properties['widths'][i]/2, time[peak_idx] - properties['widths'][i]/2
        if width_max >= time[-1] or width_min <= time[0]:
            width_overrange = True
            min_acceptable_width = min(time[peak_idx]-time[0], time[-1]-time[peak_idx])
        peak_centers.append({
            'index': peak_idx,
            'time': time[peak_idx],
            'height': signal[peak_idx],
            'prominence': properties['prominences'][i] if 'prominences' in properties else None,
            'width': properties['widths'][i] if 'widths' in properties and width_overrange == False else min_acceptable_width
        })
    
    return peak_centers

def fit_single_gaussian(peak_region_data, peak_center_info, baseline_noise):
    """
    Fit a Gaussian function to a single peak within a region.
    
    Args:
        peak_region_data: 2D numpy array with [time, signal] for the peak region
        peak_center_info: Dictionary with peak center information
        baseline_noise: Standard deviation of baseline noise for error estimation
        
    Returns:
        Dictionary with fitting results and statistics
    """
    time = peak_region_data[0]
    signal = peak_region_data[1]
    center_idx = peak_center_info['index']
    
    # Extract a window around the peak center for fitting
    # Use a window of ~4 sigma (estimated from peak width)
    window_size = int(peak_center_info.get('width', 10) * 8 ) if 'width' in peak_center_info else 80
    window_size = max(10, min(window_size, len(time)))  # Reasonable bounds
    
    start_idx = max(0, center_idx - window_size // 2)
    end_idx = min(len(time), center_idx + window_size // 2)
    
    fit_time = time[start_idx:end_idx]
    fit_signal = signal[start_idx:end_idx]
    
    # Initial parameter guesses
    amplitude_guess = peak_center_info['height']
    center_guess = peak_center_info['time']
    sigma_guess = (fit_time[-1] - fit_time[0]) / 6  # Rough estimate
    
    # Parameter bounds
    bounds = (
        [0, fit_time[0], 0.1],  # lower bounds: amplitude, center, sigma
        [amplitude_guess * 2, fit_time[-1], (fit_time[-1] - fit_time[0]) / 2]  # upper bounds
    )
    
    try:
        # Perform Gaussian fit
        popt, pcov = curve_fit(
            gaussian, 
            fit_time, 
            fit_signal, 
            p0=[amplitude_guess, center_guess, sigma_guess],
            bounds=bounds,
            maxfev=10000
        )
        
        # Calculate fitted curve and residuals
        fitted_signal = gaussian(fit_time, *popt)
        residuals = fit_signal - fitted_signal
        
        # Calculate statistics
        chi_squared = np.sum((residuals / baseline_noise) ** 2)
        dof = len(fit_time) - 3  # degrees of freedom
        reduced_chi2 = chi_squared / dof if dof > 0 else np.inf
        p_value = 1 - chi2.cdf(chi_squared, dof) if dof > 0 else 0
        
        # Calculate area under the peak
        area = simpson(fitted_signal, fit_time)
        
        # Calculate FWHM
        fwhm = 2.355 * popt[2]  # FWHM = 2.355 * sigma
        
        return {
            'success': True,
            'parameters': popt,  # [amplitude, center, sigma]
            'covariance': pcov,
            'fitted_curve': np.vstack((fit_time, fitted_signal)),
            'residuals': residuals,
            'area': area,
            'fwhm': fwhm,
            'chi_squared': chi_squared,
            'reduced_chi2': reduced_chi2,
            'p_value': p_value,
            'dof': dof,
            'fit_window': (start_idx, end_idx)
        }
        
    except Exception as e:
        print(f"Gaussian fit failed for peak at {peak_center_info['time']:.3f}s: {e}")
        return {
            'success': False,
            'error': str(e),
            'parameters': None
        }
def multi_gaussian(x, *params):
    """
    Sum of multiple Gaussian functions for fitting overlapping peaks.
    
    Args:
        x: Independent variable (time)
        *params: Alternating [amplitude, center, sigma] for each peak
        
    Returns:
        Sum of Gaussian functions
    """
    result = np.zeros_like(x)
    for i in range(0, len(params), 3):
        amplitude = params[i]
        center = params[i + 1]
        sigma = params[i + 2]
        result += gaussian(x, amplitude, center, sigma)
    return result

def fit_multi_gaussian(peak_region_data, peak_centers, baseline_noise, max_peaks=100):
    """
    Fit multiple overlapping Gaussian peaks within a region.
    
    Args:
        peak_region_data: 2D numpy array with [time, signal] for the peak region
        peak_centers: List of peak center information dictionaries
        baseline_noise: Standard deviation of baseline noise for error estimation
        max_peaks: Maximum number of peaks to attempt fitting
        
    Returns:
        Dictionary with multi-Gaussian fitting results including individual peak arrays
    """
    time = peak_region_data[0]
    signal = peak_region_data[1]
    
    # Limit the number of peaks to fit
    if len(peak_centers) > max_peaks:
        print(f"Warning: Limiting fitting to {max_peaks} peaks out of {len(peak_centers)} detected")
        peak_centers = sorted(peak_centers, key=lambda x: x['height'], reverse=True)[:max_peaks]
    
    # Prepare initial parameters and bounds
    initial_params = []
    lower_bounds = []
    upper_bounds = []
    for peak in peak_centers:
        # Initial guesses
        amplitude_guess = peak['height']
        center_guess = peak['time']
        
        # Estimate sigma from peak width if available, otherwise use reasonable guess
        if 'width' in peak and peak['width'] is not None:
            sigma_guess = peak['width'] / 2.355  # Convert FWHM to sigma
        else:
            sigma_guess = (time[-1] - time[0]) / (len(peak_centers) * 6)
        
        initial_params.extend([amplitude_guess, center_guess, sigma_guess])
        
        # Parameter bounds
        lower_bounds.extend([0, time[0], 0.1])
        upper_bounds.extend([
            amplitude_guess * 3, 
            time[-1], 
            (time[-1] - time[0]) / 2
        ])
    
    try:
        # Perform multi-Gaussian fit
        popt, pcov = curve_fit(
            multi_gaussian, 
            time, 
            signal, 
            p0=initial_params,
            bounds=(lower_bounds, upper_bounds),
            maxfev=10000
        )
        
        # Calculate fitted curve and residuals
        fitted_signal = multi_gaussian(time, *popt)
        residuals = signal - fitted_signal
        
        # Calculate statistics
        chi_squared = np.sum((residuals / baseline_noise) ** 2)
        dof = len(time) - len(popt)  # degrees of freedom
        reduced_chi2 = chi_squared / dof if dof > 0 else np.inf
        p_value = 1 - chi2.cdf(chi_squared, dof) if dof > 0 else 0
        
        # Extract individual peak parameters and calculate their properties
        individual_peaks = []
        individual_peak_arrays = []  # List to store numpy arrays for each peak
        total_area = 0
        
        for i in range(0, len(popt), 3):
            amplitude, center, sigma = popt[i:i+3]
            
            # Calculate individual peak signal
            individual_signal = gaussian(time, amplitude, center, sigma)
            
            # Create numpy array for this individual peak [time, signal]
            individual_array = np.vstack((time, individual_signal))
            individual_peak_arrays.append(individual_array)
            
            # Calculate individual peak area
            area = simpson(individual_signal, time)
            total_area += area
            
            # Calculate FWHM
            fwhm = 2.355 * sigma
            
            # Estimate parameter uncertainties from covariance matrix
            if pcov is not None and not np.isinf(pcov).any():
                param_errors = np.sqrt(np.diag(pcov))
                amplitude_err = param_errors[i] if i < len(param_errors) else np.nan
                center_err = param_errors[i+1] if i+1 < len(param_errors) else np.nan
                sigma_err = param_errors[i+2] if i+2 < len(param_errors) else np.nan
            else:
                amplitude_err = center_err = sigma_err = np.nan
            
            individual_peaks.append({
                'amplitude': amplitude,
                'center': center,
                'sigma': sigma,
                'fwhm': fwhm,
                'area': area,
                'area_percentage': (area / total_area) * 100 if total_area > 0 else 0,
                'amplitude_error': amplitude_err,
                'center_error': center_err,
                'sigma_error': sigma_err,
                'peak_index': i // 3 + 1  # Add peak index for reference
            })
        
        return {
            'success': True,
            'parameters': popt,
            'covariance': pcov,
            'fitted_curve': np.vstack((time, fitted_signal)),  # Combined fit
            'individual_peaks': individual_peaks,  # Peak parameters
            'individual_peak_arrays': individual_peak_arrays,  # Numpy arrays for each peak
            'total_area': total_area,
            'residuals': residuals,
            'chi_squared': chi_squared,
            'reduced_chi2': reduced_chi2,
            'p_value': p_value,
            'dof': dof,
            'num_peaks_fitted': len(peak_centers),
            'time_array': time,  # Original time array for reference
            'signal_array': signal  # Original signal array for reference
        }
        
    except Exception as e:
        print(f"Multi-Gaussian fit failed: {e}")
        return {
            'success': False,
            'error': str(e),
            'parameters': None,
            'individual_peaks': [],
            'individual_peak_arrays': []  # Empty list even on failure
        }
def analyze_chromatogram(chrom):
    base = baseline(chrom, half_window=30, threshold=0.002, min_length=2, min_fwhm = None, smooth_half_window=5)
    corrected_array = base['corrected_arr']
    peak_mask = base['peak_mask']
    baseline_mask = base['baseline_mask']
    threshold = base['noise_threshold']
    original_time = corrected_array[0, :]
    original_signal = corrected_array[1, :]
    continuous_regions_list = find_peak_regions(peak_mask.T)
    
    peak_info = pd.DataFrame()
    
    # Initialize the fitted peaks array with zeros
    fitted_peaks_all = np.zeros_like(original_signal)

    for start, stop in continuous_regions_list:
        # Extract the region from the ORIGINAL arrays (not masked)
        region_time = original_time[start:stop]
        region_signal = original_signal[start:stop]
        
        # Identify peak centers using the original data
        peak_centers = identify_peak_centers(np.vstack((region_time, region_signal)), 
                                            noise_threshold=threshold)
        
        # Fit multi-Gaussian using the original data
        gauss_fit = fit_multi_gaussian(np.vstack((region_time, region_signal)), 
                                      peak_centers, 
                                      baseline_noise=base['baseline_noise'])  # Use actual baseline noise
        
        if gauss_fit['success']:
            # Add peak area to peak_info dictionary
            for peak in gauss_fit['individual_peaks']:
                #round peak center to nearest multiple of 0.2 so that it aligns with the chromatogram collection interval
                retention_time = np.round(np.round(peak['center']/0.2)*0.2, decimals = 1)
                print(retention_time)
                area = peak['area']
                peak_info[retention_time] = [area]
            # Get the fitted curve for this region
            fitted_curve = gauss_fit['fitted_curve']
            
            # Store in the correct positions of the full array
            fitted_peaks_all[start:stop] = fitted_curve[1, :]
            
            # Plot individual region for debugging
            plt.figure(figsize=(10, 6))
            plt.plot(region_time, region_signal, 'b-', label='Original', linewidth=2)
            plt.plot(fitted_curve[0, :], fitted_curve[1, :], 'r-', label='Fitted', linewidth=2)
            
            # Plot individual peaks
            for i, peak_array in enumerate(gauss_fit['individual_peak_arrays']):
                plt.plot(peak_array[0, :], peak_array[1, :], '--', 
                        label=f'Peak {i+1}', alpha=0.7)
            
            plt.legend()
            plt.title(f'Region {start}-{stop}')
            plt.show()
    residuals = original_signal - fitted_peaks_all
    plt.figure(figsize=(100, 8))
    plt.plot(original_time, original_signal, 'b-', label='Original Signal', linewidth=1, zorder = 1)
    plt.scatter(original_time, fitted_peaks_all, s = 1, label='Gaussian Fit', color = 'red', zorder = 2)
    plt.xlabel('Retention Time')
    plt.ylabel('Signal Intensity')
    plt.legend()
    plt.title('Full Chromatogram Gaussian Peak Fitting Results')
    for center in peak_info.columns:
        plt.axvline(x = center, color = 'green', linestyle='--', linewidth=1)
    plt.show()
    
    plt.plot(original_time, residuals)
    plt.title('Residuals')
    return peak_info




if __name__ == "__main__":
    c1 = cc.Chromatogram(r"C:\Users\aengstrom\AutoGC\RB\validation\08\processed\ezchrom_outputs\cdf\rbsh01ndat-Back Signal.cdf",dataformat="cdf")
    chrom = c1.chromatogram
    peak_info = analyze_chromatogram(chrom)
    
                                       
