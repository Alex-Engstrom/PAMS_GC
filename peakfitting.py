# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 16:06:35 2025

@author: aengstrom
"""
import chromclass as cc
from pybaselines import Baseline
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.signal import savgol_filter, find_peaks
from scipy.optimize import curve_fit
from scipy.integrate import simpson
from scipy.stats import chi2

# Fundamental constants for pA to ion count conversion
ELEMENTARY_CHARGE = 1.60217662e-19  # Coulombs per electron
PICOAMPS_PER_AMP = 1e12             # pA per A

#%%
def baseline(data, half_window=20, threshold=0.003, min_length=5, smooth_half_window=5):
    """
    Baseline correction using pybaselines fastchrom
    Returns numpy arrays with retention time and data
    """
    x = data[0]  # Retention times
    y = data[1]  # Signal values
    
    try:
        # Convert x to minutes if it's in seconds
        if max(x) > 100:  # If values are large, assume seconds
            x_minutes = x / 60.0
        else:
            x_minutes = x
        
        # Initialize baseline fitter
        baseline_fitter = Baseline(x_minutes)
        
        # Apply fastchrom method
        baseline_vals, params = baseline_fitter.fastchrom(
            y,
            half_window=half_window,
            threshold=threshold,
            min_length=min_length,
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
        original_chrom = np.vstack((x_minutes, y))  # Time (min) + original signal
        baseline_arr = np.vstack((x_minutes, baseline_vals))  # Time + baseline
        corrected_arr = np.vstack((x_minutes, corrected_y))  # Time + corrected signal
        
        return {
            'original_chrom': original_chrom,
            'baseline_arr': baseline_arr,
            'corrected_arr': corrected_arr,
            'baseline_mask': baseline_mask,
            'peak_mask': peak_mask,
            'baseline_noise': baseline_noise,  # 1-sigma noise
            'noise_threshold': noise_threshold,  # 3-sigma threshold
            'params': params
        }
        
    except Exception as e:
        print(f"fastchrom baseline failed: {e}")
        # Fallback approach - return arrays with zero baseline
        x_minutes = data[0] / 60.0 if max(data[0]) > 100 else data[0]
        zero_baseline = np.zeros_like(data[1])
        
        return {
            'original_chrom': np.vstack((x_minutes, data[1])),
            'baseline_arr': np.vstack((x_minutes, zero_baseline)),
            'corrected_arr': np.vstack((x_minutes, data[1])),
            'baseline_mask': np.ones_like(data[1], dtype=bool),
            'peak_mask': np.zeros_like(data[1], dtype=bool),
            'baseline_noise': 0.001,
            'noise_threshold': 0.003,
            'params': {}
        }

def gaussian(x, amplitude, center, sigma):
    """Gaussian function for peak fitting"""
    return amplitude * np.exp(-((x - center) ** 2) / (2 * sigma ** 2))

def picoamps_to_ion_counts_per_second(pA_signal):
    """Convert picoamp signal to ion counts per second"""
    amps = pA_signal / PICOAMPS_PER_AMP
    electrons_per_second = amps / ELEMENTARY_CHARGE
    return electrons_per_second  # ions per second

def ion_counts_per_second_to_picoamps(ions_per_second):
    """Convert ion counts per second back to picoamps"""
    amps = ions_per_second * ELEMENTARY_CHARGE
    pA_signal = amps * PICOAMPS_PER_AMP
    return pA_signal

def calculate_uncertainty_pA_corrected(pA_signal, baseline_noise_pA=0.001, acquisition_time=0.2):
    """
    Calculate proper uncertainty for pA signals
    Includes both Poisson counting statistics and baseline noise
    """
    # Convert to ion counts per second
    ions_per_second = picoamps_to_ion_counts_per_second(pA_signal)
    baseline_ions_per_second = picoamps_to_ion_counts_per_second(baseline_noise_pA)
    
    # Convert to total ions in acquisition time
    total_ions = ions_per_second * acquisition_time
    baseline_total_ions = baseline_ions_per_second * acquisition_time
    
    # Combined uncertainty: sqrt(Poisson + baseline^2) for total counts
    uncertainty_total_ions = np.sqrt(np.maximum(total_ions, 0) + baseline_total_ions**2)
    
    # Convert back to uncertainty in pA (per second basis)
    uncertainty_ions_per_second = uncertainty_total_ions / acquisition_time
    uncertainty_pA = ion_counts_per_second_to_picoamps(uncertainty_ions_per_second)
    
    return uncertainty_pA

def calculate_chi_squared_corrected(observed_pA, expected_pA, uncertainties_pA, n_params=3):
    """
    Calculate chi-squared with proper uncertainty handling
    """
    # Calculate chi-squared
    residuals = observed_pA - expected_pA
    chi_sq = np.sum((residuals / uncertainties_pA) ** 2)
    
    # Degrees of freedom
    dof = len(observed_pA) - n_params
    reduced_chi_sq = chi_sq / dof if dof > 0 else np.nan
    
    # P-value
    p_value = 1 - chi2.cdf(chi_sq, dof) if dof > 0 else np.nan
    
    # RMS residuals
    rms_residuals = np.sqrt(np.mean(residuals ** 2))
    
    # Fit quality assessment
    if not np.isfinite(reduced_chi_sq):
        fit_quality = 'undefined'
    elif reduced_chi_sq < 1.5:
        fit_quality = 'excellent'
    elif reduced_chi_sq < 3.0:
        fit_quality = 'good'
    elif reduced_chi_sq < 5.0:
        fit_quality = 'fair'
    else:
        fit_quality = 'poor'
    
    return {
        'chi_squared': chi_sq,
        'reduced_chi_squared': reduced_chi_sq,
        'degrees_of_freedom': dof,
        'p_value': p_value,
        'rms_residuals': rms_residuals,
        'fit_quality': fit_quality
    }

def find_contiguous_regions(mask):
    """
    Find contiguous True regions in a boolean mask
    """
    regions = []
    in_region = False
    start_idx = 0
    
    for i, value in enumerate(mask):
        if value and not in_region:
            # Start of new region
            in_region = True
            start_idx = i
        elif not value and in_region:
            # End of current region
            in_region = False
            if i - start_idx > 1:  # Only add regions with at least 2 points
                regions.append((start_idx, i))
    
    # Handle case where region continues to end of array
    if in_region and len(mask) - start_idx > 1:
        regions.append((start_idx, len(mask)))
    
    return regions

def improve_peak_edge_detection(time, signal, center_idx, baseline_noise, min_points=10):
    """
    Improved peak edge detection that captures the full peak shape
    Returns integer indices for proper slicing
    """
    # Calculate first and second derivatives for better edge detection
    derivative = np.gradient(signal, time)
    second_derivative = np.gradient(derivative, time)
    
    # Find left edge
    left_edge = int(center_idx)  # Ensure integer
    baseline_level = np.median(signal[max(0, center_idx-50):max(0, center_idx-10)])
    threshold = baseline_level + 2 * baseline_noise
    
    # Move left until we hit baseline or inflection point
    for i in range(center_idx, max(0, center_idx-200), -1):
        i = int(i)  # Ensure integer
        if i <= 0:
            left_edge = 0
            break
        
        # Check multiple conditions for edge detection
        signal_condition = signal[i] <= threshold
        derivative_condition = abs(derivative[i]) < baseline_noise * 0.5
        second_derivative_condition = second_derivative[i] > 0  # Concave up
        
        if signal_condition and (derivative_condition or second_derivative_condition):
            left_edge = i
            # Move a bit further to ensure we capture the full rise
            left_edge = max(0, left_edge - min(5, left_edge))
            break
    
    # Find right edge
    right_edge = int(center_idx)  # Ensure integer
    for i in range(center_idx, min(len(signal), center_idx+200)):
        i = int(i)  # Ensure integer
        if i >= len(signal) - 1:
            right_edge = len(signal) - 1
            break
        
        # Check multiple conditions for edge detection
        signal_condition = signal[i] <= threshold
        derivative_condition = abs(derivative[i]) < baseline_noise * 0.5
        second_derivative_condition = second_derivative[i] > 0  # Concave up
        
        if signal_condition and (derivative_condition or second_derivative_condition):
            right_edge = i
            # Move a bit further to ensure we capture the full fall
            right_edge = min(len(signal) - 1, right_edge + 5)
            break
    
    # Ensure minimum number of points for fitting
    current_width = right_edge - left_edge
    if current_width < min_points:
        # Expand symmetrically
        expansion = min_points - current_width
        left_edge = max(0, left_edge - expansion // 2)
        right_edge = min(len(signal) - 1, right_edge + expansion // 2)
    
    # Final integer conversion
    left_edge = int(left_edge)
    right_edge = int(right_edge)
    
    return left_edge, right_edge


def fit_gaussian_to_full_peak(time, signal_pA, center_idx, baseline_noise_pA=0.001, 
                            acquisition_time=0.2, max_iterations=3):
    """
    Iterative Gaussian fitting that ensures the full peak is captured
    """
    best_fit = None
    best_chi_sq = float('inf')
    
    for iteration in range(max_iterations):
        # Adjust edge detection based on iteration
        expansion_factor = 1 + iteration * 0.5  # Expand search area each iteration
        
        left_edge, right_edge = improve_peak_edge_detection(
            time, signal_pA, center_idx, baseline_noise_pA,
            min_points=15 * expansion_factor
        )
        
        # Extract the peak region - ensure integer indices
        left_edge = int(left_edge)
        right_edge = int(right_edge)
        peak_indices = slice(left_edge, right_edge + 1)
        
        x_fit = time[peak_indices]
        y_fit_pA = signal_pA[peak_indices]
        
        # Calculate uncertainties
        uncertainties_pA = np.array([calculate_uncertainty_pA_corrected(s, baseline_noise_pA, acquisition_time) 
                                   for s in y_fit_pA])
        
        # Initial parameter guesses with better estimates
        amplitude_guess = np.max(y_fit_pA) - np.min(y_fit_pA)
        center_guess = time[center_idx] if left_edge <= center_idx <= right_edge else np.mean(x_fit)
        
        # Estimate sigma from peak width (more robust)
        half_max = (np.max(y_fit_pA) + np.min(y_fit_pA)) / 2
        above_half = np.where(y_fit_pA > half_max)[0]
        if len(above_half) > 1:
            fwhm = x_fit[above_half[-1]] - x_fit[above_half[0]]
            sigma_guess = fwhm / (2 * np.sqrt(2 * np.log(2)))
        else:
            sigma_guess = (x_fit[-1] - x_fit[0]) / 4
        
        # Wider bounds to accommodate full peaks
        bounds = (
            [0.1 * amplitude_guess, x_fit[0], 0.01 * sigma_guess],
            [10 * amplitude_guess, x_fit[-1], 10 * sigma_guess]
        )
        
        try:
            # Fit with more iterations and better method
            popt, pcov = curve_fit(
                gaussian, x_fit, y_fit_pA,
                p0=[amplitude_guess, center_guess, sigma_guess],
                sigma=uncertainties_pA,
                absolute_sigma=True,
                bounds=bounds,
                maxfev=5000,
                method='trf'  # Trust Region Reflective algorithm
            )
            
            # Calculate fitted curve
            y_fitted_pA = gaussian(x_fit, *popt)
            
            # Calculate area
            area_pA_min = simpson(y_fitted_pA, x_fit)
            area_ions = area_pA_min * 60 * picoamps_to_ion_counts_per_second(1)
            center_ion_count = picoamps_to_ion_counts_per_second(popt[0]) * acquisition_time
            
            # Calculate chi-squared
            chi_sq_stats = calculate_chi_squared_corrected(y_fit_pA, y_fitted_pA, uncertainties_pA)
            
            # Parameter uncertainties
            perr = np.sqrt(np.diag(pcov)) if pcov is not None else np.zeros(3)
            
            current_fit = {
                'amplitude_pA': popt[0],
                'center': popt[1],
                'sigma': popt[2],
                'fwhm': 2.355 * popt[2],
                'area_pA_min': area_pA_min,
                'amplitude_ions_per_acq': center_ion_count,
                'total_ions': area_ions,
                'amplitude_err_pA': perr[0],
                'center_err': perr[1],
                'sigma_err': perr[2],
                'chi_squared': chi_sq_stats['chi_squared'],
                'reduced_chi_squared': chi_sq_stats['reduced_chi_squared'],
                'degrees_of_freedom': chi_sq_stats['degrees_of_freedom'],
                'p_value': chi_sq_stats['p_value'],
                'rms_residuals_pA': chi_sq_stats['rms_residuals'],
                'fit_quality': chi_sq_stats['fit_quality'],
                'fitted_curve_pA': y_fitted_pA,
                'fit_x': x_fit,
                'actual_data_pA': y_fit_pA,
                'uncertainties_pA': uncertainties_pA,
                'left_edge': left_edge,
                'right_edge': right_edge
            }
            
            # Keep the best fit (lowest chi-squared)
            if current_fit['reduced_chi_squared'] < best_chi_sq:
                best_fit = current_fit
                best_chi_sq = current_fit['reduced_chi_squared']
                
        except Exception as e:
            print(f"Gaussian fit iteration {iteration+1} failed: {e}")
            continue
    
    return best_fit

def calculate_trapezoidal_area_corrected(time, signal_pA, left_edge, right_edge, acquisition_time=0.2):
    """
    Calculate trapezoidal area with proper unit conversion to total ions
    """
    if left_edge is None or right_edge is None or right_edge <= left_edge:
        return 0.0, 0.0
    
    start_idx = max(0, left_edge)
    end_idx = min(len(signal_pA), right_edge + 1)
    
    if end_idx <= start_idx:
        return 0.0, 0.0
    
    # Trapezoidal integration in pA·min
    area_pA_min = np.trapezoid(signal_pA[start_idx:end_idx], time[start_idx:end_idx])
    
    # Convert to total ions: pA·min → ions
    total_ions = area_pA_min * 60 * picoamps_to_ion_counts_per_second(1)
    
    return total_ions, area_pA_min

def refine_peak_center(signal, approximate_center, window=5):
    """
    Refine peak center position by finding local maximum in a window
    Returns integer index
    """
    start = max(0, int(approximate_center) - window)
    end = min(len(signal), int(approximate_center) + window + 1)
    
    window_signal = signal[start:end]
    local_max_idx = np.argmax(window_signal)
    
    return int(start + local_max_idx)

def find_peaks_scipy_strategy(time, signal_smooth, baseline_noise, 
                            min_snr, min_prominence, min_width,
                            acquisition_time, max_reduced_chi_sq):
    """
    Primary peak detection using scipy's find_peaks with adaptive parameters
    """
    peaks_info = []
    
    # Adaptive prominence calculation based on signal range
    signal_range = np.max(signal_smooth) - np.min(signal_smooth)
    adaptive_prominence = max(min_prominence * signal_range, baseline_noise * min_snr)
    
    # Try multiple width settings to catch different peak types
    width_settings = [min_width, min_width * 2, min_width * 3]
    
    for width_setting in width_settings:
        try:
            peak_indices, properties = find_peaks(
                signal_smooth,
                height=baseline_noise * min_snr,
                prominence=adaptive_prominence,
                width=width_setting,
                rel_height=0.5
            )
            
            for i, peak_idx in enumerate(peak_indices):
                peak_idx = int(peak_idx)  # Ensure integer
                # Skip if we already found this peak
                if any(abs(p['center_index'] - peak_idx) < 5 for p in peaks_info):
                    continue
                
                # Estimate peak edges using the improved function
                left_edge, right_edge = improve_peak_edge_detection(
                    time, signal_smooth, peak_idx, baseline_noise
                )
                
                # Fit Gaussian using the improved fitting function
                gaussian_fit = fit_gaussian_to_full_peak(
                    time, signal_smooth, peak_idx, baseline_noise, acquisition_time
                )
                
                if (gaussian_fit and gaussian_fit['reduced_chi_squared'] <= max_reduced_chi_sq and
                    gaussian_fit['amplitude_pA'] >= baseline_noise * min_snr):
                    
                    peak_info = create_peak_info(time, signal_smooth, peak_idx, left_edge, right_edge,
                                               gaussian_fit, baseline_noise, acquisition_time)
                    peaks_info.append(peak_info)
                    
        except Exception as e:
            print(f"Scipy peak finding failed for width {width_setting}: {e}")
            continue
    
    return peaks_info

def create_peak_info(time, signal, center_idx, left_edge, right_edge,
                   gaussian_fit, baseline_noise, acquisition_time):
    """
    Create a standardized peak information dictionary
    """
    trapezoidal_ions, trapezoidal_pA_min = calculate_trapezoidal_area_corrected(
        time, signal, left_edge, right_edge, acquisition_time
    )
    
    snr = gaussian_fit['amplitude_pA'] / baseline_noise
    
    return {
        'center_index': int(center_idx),
        'center_time': float(time[center_idx]),
        'center_signal_pA': float(signal[center_idx]),
        'center_signal_ions': float(picoamps_to_ion_counts_per_second(signal[center_idx]) * acquisition_time),
        'left_edge_index': int(left_edge),
        'left_edge_time': float(time[left_edge]),
        'right_edge_index': int(right_edge),
        'right_edge_time': float(time[right_edge]),
        'peak_width': float(time[right_edge] - time[left_edge]),
        'trapezoidal_area_ions': float(trapezoidal_ions),
        'trapezoidal_area_pA_min': float(trapezoidal_pA_min),
        'gaussian_fit': gaussian_fit,
        'fitted_amplitude_pA': gaussian_fit['amplitude_pA'],
        'fitted_center': gaussian_fit['center'],
        'fitted_sigma': gaussian_fit['sigma'],
        'fitted_fwhm': gaussian_fit['fwhm'],
        'gaussian_area_pA_min': gaussian_fit['area_pA_min'],
        'gaussian_area_ions': gaussian_fit['total_ions'],
        'chi_squared': gaussian_fit['chi_squared'],
        'reduced_chi_squared': gaussian_fit['reduced_chi_squared'],
        'degrees_of_freedom': gaussian_fit['degrees_of_freedom'],
        'p_value': gaussian_fit['p_value'],
        'rms_residuals_pA': gaussian_fit['rms_residuals_pA'],
        'fit_quality': gaussian_fit['fit_quality'],
        'fitted_amplitude_ions': gaussian_fit['amplitude_ions_per_acq'],
        'signal_to_noise': snr,
        'baseline_noise_pA': baseline_noise
    }

def find_peaks_derivative_strategy(time, signal_smooth, baseline_noise, existing_peaks,
                                 min_snr, acquisition_time, max_reduced_chi_sq):
    """
    Secondary peak detection using derivative analysis to find peaks
    that might be missed by scipy's algorithm
    """
    additional_peaks = []
    
    # Calculate first and second derivatives
    derivative = np.gradient(signal_smooth, time)
    second_derivative = np.gradient(derivative, time)
    
    # Find zero crossings in first derivative (potential peak centers)
    zero_crossings = np.where(np.diff(np.sign(derivative)))[0]
    
    for crossing in zero_crossings:
        crossing = int(crossing)  # Ensure integer
        # Skip if this is near an existing peak
        if any(abs(p['center_index'] - crossing) < 10 for p in existing_peaks):
            continue
        
        # Check if this is a maximum (negative second derivative)
        if crossing > 0 and crossing < len(second_derivative) - 1:
            if second_derivative[crossing] < 0:  # This indicates a maximum
                # Refine peak center around the zero crossing
                refined_center = refine_peak_center(signal_smooth, crossing, 5)
                
                # Estimate edges using the improved function
                left_edge, right_edge = improve_peak_edge_detection(
                    time, signal_smooth, refined_center, baseline_noise
                )
                
                # Check if this peak is significant
                peak_height = signal_smooth[refined_center]
                if peak_height >= baseline_noise * min_snr:
                    # Fit Gaussian using the improved function
                    gaussian_fit = fit_gaussian_to_full_peak(
                        time, signal_smooth, refined_center, baseline_noise, acquisition_time
                    )
                    
                    if (gaussian_fit and gaussian_fit['reduced_chi_squared'] <= max_reduced_chi_sq):
                        peak_info = create_peak_info(time, signal_smooth, refined_center, left_edge, right_edge,
                                                   gaussian_fit, baseline_noise, acquisition_time)
                        additional_peaks.append(peak_info)
    
    return additional_peaks

def remove_duplicate_peaks(peaks_info, rt_tolerance=0.01):
    """
    Remove duplicate peaks that are too close in retention time
    """
    unique_peaks = []
    seen_centers = set()
    
    for peak in sorted(peaks_info, key=lambda x: x['fitted_amplitude_pA'], reverse=True):
        # Check if this peak is too close to one we've already seen
        is_duplicate = False
        for seen_center in seen_centers:
            if abs(peak['fitted_center'] - seen_center) < rt_tolerance:
                is_duplicate = True
                break
        
        if not is_duplicate:
            unique_peaks.append(peak)
            seen_centers.add(peak['fitted_center'])
    
    return unique_peaks

def find_all_significant_peaks(time, signal_smooth, baseline_noise, min_snr=3.0, 
                             min_prominence=0.01, acquisition_time=0.2):
    """
    Comprehensive peak detection using multiple strategies to ensure all significant peaks are found
    """
    all_peak_indices = set()
    
    # Strategy 1: scipy find_peaks with multiple parameter sets
    peak_strategies = [
        # Standard parameters
        {'height': baseline_noise * min_snr, 'prominence': min_prominence, 'width': 3},
        # Broader peaks
        {'height': baseline_noise * min_snr, 'prominence': min_prominence * 0.5, 'width': 5},
        # Very broad peaks
        {'height': baseline_noise * min_snr, 'prominence': min_prominence * 0.3, 'width': 8},
        # High sensitivity
        {'height': baseline_noise * (min_snr * 0.7), 'prominence': min_prominence * 0.2, 'width': 2},
    ]
    
    for params in peak_strategies:
        try:
            peak_indices, _ = find_peaks(signal_smooth, **params)
            all_peak_indices.update(peak_indices)
        except:
            continue
    
    # Strategy 2: Find peaks using gradient analysis (for very broad peaks)
    gradient = np.gradient(signal_smooth, time)
    second_gradient = np.gradient(gradient, time)
    
    # Find potential peaks where gradient changes sign (zero crossings)
    zero_crossings = np.where(np.diff(np.sign(gradient)))[0]
    
    for crossing in zero_crossings:
        if crossing < 5 or crossing > len(signal_smooth) - 5:
            continue
        
        # Check if this is a maximum (negative second derivative)
        if second_gradient[crossing] < 0:
            # Refine the peak position
            window = min(10, crossing, len(signal_smooth) - crossing - 1)
            local_max_idx = np.argmax(signal_smooth[crossing-window:crossing+window+1])
            refined_peak = crossing - window + local_max_idx
            
            # Check if peak is significant
            if signal_smooth[refined_peak] > baseline_noise * min_snr:
                all_peak_indices.add(refined_peak)
    
    # Strategy 3: Manual peak finding for very large peaks that might be missed
    # Look for regions where signal is significantly above baseline
    above_threshold = signal_smooth > baseline_noise * min_snr * 2
    regions = find_contiguous_regions(above_threshold)
    
    for start, end in regions:
        if end - start < 3:  # Skip very small regions
            continue
        
        # Find local maximum in this region
        region_signal = signal_smooth[start:end]
        local_max_idx = np.argmax(region_signal)
        peak_idx = start + local_max_idx
        
        # Only add if not already found and significant
        if (peak_idx not in all_peak_indices and 
            signal_smooth[peak_idx] > baseline_noise * min_snr):
            all_peak_indices.add(peak_idx)
    
    # Convert to sorted list and remove duplicates that are too close
    peak_indices = sorted(all_peak_indices)
    filtered_peaks = []
    
    for i, peak_idx in enumerate(peak_indices):
        # Skip peaks that are too close to previous ones
        if i > 0 and abs(peak_idx - peak_indices[i-1]) < 5:
            continue
        filtered_peaks.append(peak_idx)
    
    return filtered_peaks

def enhanced_find_peaks_complete(results, smooth_window=21, smooth_order=3, 
                               min_snr=3.0, min_prominence=0.01, min_width=3,
                               acquisition_time=0.2, max_reduced_chi_sq=10.0):
    """
    Complete enhanced peak finding with improved edge detection and fitting
    Returns both the peaks info and all detected indices for debugging
    """
    time = results['corrected_arr'][0]
    signal = results['corrected_arr'][1]
    baseline_noise = results['baseline_noise']
    
    # Smooth the signal
    if len(signal) > smooth_window:
        signal_smooth = savgol_filter(signal, smooth_window, smooth_order)
    else:
        signal_smooth = signal
    
    # Find all significant peaks using comprehensive detection
    all_peak_indices = find_all_significant_peaks(
        time, signal_smooth, baseline_noise, min_snr, min_prominence, acquisition_time
    )
    
    print(f"Found {len(all_peak_indices)} potential peak centers")
    
    # Fit each detected peak
    peaks_info = []
    successful_fits = 0
    
    for peak_idx in all_peak_indices:
        peak_idx = int(peak_idx)
        
        # Estimate peak edges using the improved function
        left_edge, right_edge = improve_peak_edge_detection(
            time, signal_smooth, peak_idx, baseline_noise
        )
        
        # Fit Gaussian using the improved fitting function
        gaussian_fit = fit_gaussian_to_full_peak(
            time, signal_smooth, peak_idx, baseline_noise, acquisition_time
        )
        
        if (gaussian_fit and gaussian_fit['reduced_chi_squared'] <= max_reduced_chi_sq and
            gaussian_fit['amplitude_pA'] >= baseline_noise * min_snr):
            
            peak_info = create_peak_info(
                time, signal_smooth, peak_idx, left_edge, right_edge,
                gaussian_fit, baseline_noise, acquisition_time
            )
            peaks_info.append(peak_info)
            successful_fits += 1
        else:
            print(f"Peak at {time[peak_idx]:.3f} min failed fitting: "
                  f"χ²={gaussian_fit['reduced_chi_squared'] if gaussian_fit else 'N/A'}, "
                  f"SNR={signal_smooth[peak_idx]/baseline_noise:.1f}")
    
    # Remove duplicates and sort
    final_peaks = remove_duplicate_peaks(peaks_info)
    final_peaks.sort(key=lambda x: x['fitted_center'])
    
    print(f"Successful fits: {successful_fits}/{len(all_peak_indices)}")
    print(f"After filtering: {len(final_peaks)} valid peaks")
    
    return final_peaks, all_peak_indices

def manual_peak_inspection(results, min_snr=5.0):
    """
    Manual inspection function to help identify missed peaks
    """
    time = results['corrected_arr'][0]
    signal = results['corrected_arr'][1]
    baseline_noise = results['baseline_noise']
    
    # Find regions significantly above baseline
    threshold = baseline_noise * min_snr
    above_threshold = signal > threshold
    
    regions = find_contiguous_regions(above_threshold)
    
    print("Potential peak regions (manual inspection):")
    for i, (start, end) in enumerate(regions):
        if end - start < 5:  # Skip very small regions
            continue
        
        region_max = np.max(signal[start:end])
        region_center = start + np.argmax(signal[start:end])
        region_snr = region_max / baseline_noise
        
        print(f"Region {i+1}: RT={time[region_center]:.3f} min, "
              f"Max={region_max:.2f} pA, SNR={region_snr:.1f}, "
              f"Width={time[end]-time[start]:.3f} min")
    
    return regions

#%% Plotting functions

def plot_peak_detection_debug(results, all_peak_indices, peaks_info, figsize=(15, 10)):
    """
    Debug plot to show all detected peaks and which ones were kept
    """
    time = results['corrected_arr'][0]
    signal = results['corrected_arr'][1]
    signal_smooth = savgol_filter(signal, 21, 3) if len(signal) > 21 else signal
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    
    # Plot 1: Show all detected peaks
    ax1.plot(time, signal, 'b-', alpha=0.6, label='Signal', linewidth=1)
    ax1.plot(time, signal_smooth, 'g-', alpha=0.8, label='Smoothed', linewidth=1)
    
    # Mark all potential peaks
    for i, peak_idx in enumerate(all_peak_indices):
        ax1.plot(time[peak_idx], signal_smooth[peak_idx], 'ro', 
                markersize=6, alpha=0.7, label='Detected' if i == 0 else "")
    
    ax1.set_title('All Detected Peak Centers')
    ax1.set_ylabel('Signal (pA)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Show final fitted peaks
    ax2.plot(time, signal, 'b-', alpha=0.6, label='Signal', linewidth=1)
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(peaks_info)))
    for i, (peak, color) in enumerate(zip(peaks_info, colors)):
        fit = peak['gaussian_fit']
        ax2.plot(fit['fit_x'], fit['fitted_curve_pA'], color=color, linewidth=2, 
                label=f'Peak {i+1}' if i < 5 else "")
        ax2.axvline(x=peak['fitted_center'], color=color, linestyle='--', alpha=0.7)
    
    ax2.set_title('Final Fitted Peaks')
    ax2.set_xlabel('Retention Time (min)')
    ax2.set_ylabel('Signal (pA)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
def plot_baseline_results(results, title_suffix=""):
    """
    Plot the baseline correction results
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(28, 21))
    
    # Extract data from arrays
    time_min = results['original_chrom'][0]
    original_signal = results['original_chrom'][1]
    baseline_vals = results['baseline_arr'][1]
    corrected_signal = results['corrected_arr'][1]
    
    # Plot 1: Original data with baseline
    ax1.plot(time_min, original_signal, 
             label="Raw Chromatogram", color="blue", alpha=0.6)
    ax1.plot(time_min, baseline_vals, 
             label="Estimated Baseline", color="red", linewidth=2)
    ax1.set_xlabel("Retention Time (min)")
    ax1.set_ylabel("Signal (FID)")
    ax1.set_title(f"Baseline Correction {title_suffix}")
    ax1.legend()
    ax1.grid(True)
    
    # Plot 2: Baseline corrected data
    ax2.plot(time_min, corrected_signal, 
             label="Corrected Chromatogram", color="green", alpha=0.8)
    ax2.axhline(y=results['noise_threshold'], color='orange', linestyle='--',
                label=f'Noise Threshold (3σ = {results["noise_threshold"]:.6f})')
    ax2.axhline(y=results['baseline_noise'], color='red', linestyle=':',
                label=f'Baseline Noise (1σ = {results["baseline_noise"]:.6f})')
    ax2.set_xlabel("Retention Time (min)")
    ax2.set_ylabel("Corrected Signal")
    ax2.set_title("Baseline Corrected Signal with Noise Threshold")
    ax2.legend()
    ax2.grid(True)
    
    plt.tight_layout()
    plt.show()

def print_detailed_peak_info(peaks_info, peak_numbers=None):
    """
    Print detailed information about specific peaks with chi-squared statistics
    """
    if peak_numbers is None:
        peak_numbers = range(1, len(peaks_info) + 1)
    
    print("="*120)
    print("DETAILED PEAK INFORMATION WITH χ² STATISTICS")
    print("="*120)
    
    for peak_num in peak_numbers:
        if peak_num < 1 or peak_num > len(peaks_info):
            continue
            
        peak = peaks_info[peak_num - 1]
        fit = peak['gaussian_fit']
        
        print(f"\nPEAK {peak_num}:")
        print(f"  Retention Time:     {peak['fitted_center']:.6f} ± {fit['center_err']:.6f} min")
        print(f"  Amplitude (pA):     {peak['fitted_amplitude_pA']:.6f} ± {fit['amplitude_err_pA']:.6f}")
        print(f"  Amplitude (ions):   {peak['fitted_amplitude_ions']:.1f} ions/point")
        print(f"  Sigma (Width):      {peak['fitted_sigma']:.6f} ± {fit['sigma_err']:.6f} min")
        print(f"  FWHM:               {peak['fitted_fwhm']:.6f} min")
        print(f"  Gaussian Area (pA): {peak['gaussian_area_pA_min']:.6f} pA·min")
        print(f"  Gaussian Area (ions): {peak['gaussian_area_ions']:.1f} ions")
        print(f"  Trapezoidal Area:   {peak['trapezoidal_area_ions']:.1f} ions")
        print(f"  Chi-squared (χ²):   {peak['chi_squared']:.2f}")
        print(f"  Reduced χ² (χ²/ν):  {peak['reduced_chi_squared']:.2f}")
        print(f"  Degrees of Freedom: {peak['degrees_of_freedom']}")
        print(f"  P-value:            {peak['p_value']:.6f}")
        print(f"  RMS Residuals (pA): {peak['rms_residuals_pA']:.6f}")
        print(f"  Signal-to-Noise:    {peak['signal_to_noise']:.1f}")
        print(f"  Fit Quality:        {peak['fit_quality']}")
        print(f"  Baseline Noise:     {peak['baseline_noise_pA']:.6f} pA")
        print(f"  Edges:              {peak['left_edge_time']:.3f} - {peak['right_edge_time']:.3f} min")
        
        # Sanity check
        if peak['fitted_amplitude_ions'] > peak['gaussian_area_ions']:
            print(f"  WARNING: Center ions > Total ions! ({peak['fitted_amplitude_ions']:.1f} > {peak['gaussian_area_ions']:.1f})")

def plot_gaussian_fits(results, peaks_info, plot_individual_peaks=True, plot_summary=True, 
                      plot_comprehensive=True, acquisition_time=0.2, figsize=(15, 10)):
    """
    Plot Gaussian fits with residuals and statistical information
    
    Parameters:
    -----------
    results : dict
        Output from baseline() function
    peaks_info : list
        List of peak information dictionaries from Gaussian fitting
    plot_individual_peaks : bool
        Whether to plot individual peak fits
    plot_summary : bool
        Whether to plot summary statistics
    plot_comprehensive : bool
        Whether to plot comprehensive chromatogram view
    acquisition_time : float
        Acquisition time in seconds
    figsize : tuple
        Figure size
    """
    
    if not peaks_info:
        print("No peaks to plot")
        return
    
    if plot_individual_peaks:
        plot_individual_peak_fits(results, peaks_info, acquisition_time, figsize)
    
    if plot_summary:
        plot_summary_statistics(peaks_info, figsize)
        
    if plot_comprehensive:
        plot_comprehensive_chromatogram(results, peaks_info, figsize)

def plot_individual_peak_fits(results, peaks_info, acquisition_time=0.2, figsize=(15, 10)):
    """
    Plot individual peak fits with residuals and statistical information
    """
    time = results['corrected_arr'][0]
    signal = results['corrected_arr'][1]
    baseline_noise = results['baseline_noise']
    
    n_peaks = len(peaks_info)
    if n_peaks == 0:
        return
    
    n_cols = min(3, n_peaks)
    n_rows = (n_peaks + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows * 2, n_cols, figsize=(figsize[0], figsize[1] * 1.5))
    
    # Handle different subplot configurations
    if n_peaks == 1:
        axes = np.array([[axes[0]], [axes[1]]])
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    elif n_rows == 1:
        axes = axes.reshape(2, n_cols)
    
    fig.suptitle('Individual Peak Gaussian Fits with χ² Statistics', fontsize=16, fontweight='bold')
    
    colors = plt.cm.tab10(np.linspace(0, 1, n_peaks))
    
    for idx, (peak, color) in enumerate(zip(peaks_info, colors)):
        fit = peak['gaussian_fit']
        
        # Get row and column indices
        row_fit = 2 * (idx // n_cols)
        row_res = 2 * (idx // n_cols) + 1
        col = idx % n_cols
        
        # Get the axes
        if n_cols > 1:
            ax_fit = axes[row_fit, col]
            ax_res = axes[row_res, col]
        else:
            ax_fit = axes[row_fit][col]
            ax_res = axes[row_res][col]
        
        # Plot actual data and fit
        ax_fit.plot(fit['fit_x'], fit['actual_data_pA'], 'o', color=color, 
                   markersize=4, alpha=0.7, label='Actual Data')
        ax_fit.plot(fit['fit_x'], fit['fitted_curve_pA'], '-', color=color, 
                   linewidth=2, label='Gaussian Fit')
        
        # Add uncertainty bands if available
        if 'uncertainties_pA' in fit:
            ax_fit.fill_between(fit['fit_x'], 
                              fit['fitted_curve_pA'] - fit['uncertainties_pA'],
                              fit['fitted_curve_pA'] + fit['uncertainties_pA'],
                              color=color, alpha=0.2, label='Uncertainty')
        
        # Add title with statistical information
        title = (f'Peak {idx+1}: RT={peak["fitted_center"]:.3f} min\n'
                f'Ampl: {peak["fitted_amplitude_pA"]:.2f} pA, '
                f'S/N: {peak["signal_to_noise"]:.1f}\n'
                f'χ²/ν={peak["reduced_chi_squared"]:.2f}, '
                f'Quality: {peak["fit_quality"]}')
        ax_fit.set_title(title, fontsize=10)
        ax_fit.set_xlabel('Retention Time (min)')
        ax_fit.set_ylabel('Signal (pA)')
        ax_fit.legend(fontsize=8)
        ax_fit.grid(True, alpha=0.3)
        
        # Plot residuals
        residuals = fit['actual_data_pA'] - fit['fitted_curve_pA']
        ax_res.plot(fit['fit_x'], residuals, 'ko', markersize=3, alpha=0.7)
        ax_res.axhline(y=0, color='r', linestyle='-', alpha=0.7)
        
        # Add uncertainty band if available
        if 'uncertainties_pA' in fit:
            ax_res.fill_between(fit['fit_x'], -fit['uncertainties_pA'], fit['uncertainties_pA'],
                              color='gray', alpha=0.3, label='Uncertainty band')
        
        ax_res.set_title(f'Residuals (RMS: {peak["rms_residuals_pA"]:.4f} pA)')
        ax_res.set_xlabel('Retention Time (min)')
        ax_res.set_ylabel('Residual (pA)')
        ax_res.legend(fontsize=8)
        ax_res.grid(True, alpha=0.3)
    
    # Hide empty subplots
    total_subplots = n_rows * n_cols
    for idx in range(n_peaks, total_subplots):
        row_fit = 2 * (idx // n_cols)
        row_res = 2 * (idx // n_cols) + 1
        col = idx % n_cols
        
        if n_cols > 1:
            if row_fit < axes.shape[0] and col < axes.shape[1]:
                axes[row_fit, col].set_visible(False)
            if row_res < axes.shape[0] and col < axes.shape[1]:
                axes[row_res, col].set_visible(False)
        else:
            if row_fit < len(axes) and col < len(axes[row_fit]):
                axes[row_fit][col].set_visible(False)
            if row_res < len(axes) and col < len(axes[row_res]):
                axes[row_res][col].set_visible(False)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    plt.show()

def plot_summary_statistics(peaks_info, figsize=(12, 10)):
    """
    Plot summary statistics of the Gaussian fits
    """
    if not peaks_info:
        return
        
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle('Gaussian Fit Summary Statistics', fontsize=16, fontweight='bold')
    
    peak_numbers = range(1, len(peaks_info) + 1)
    
    # Plot 1: Reduced chi-squared values
    red_chi_sq = [p['reduced_chi_squared'] for p in peaks_info]
    bars1 = ax1.bar(peak_numbers, red_chi_sq, alpha=0.7, color='skyblue')
    ax1.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Ideal χ²/ν = 1')
    ax1.axhline(y=3, color='orange', linestyle=':', alpha=0.7, label='Acceptable χ²/ν = 3')
    ax1.set_xlabel('Peak Number')
    ax1.set_ylabel('Reduced Chi-squared (χ²/ν)')
    ax1.set_title('Fit Quality by Peak')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Add values on bars
    for i, (bar, val) in enumerate(zip(bars1, red_chi_sq)):
        ax1.text(bar.get_x() + bar.get_width()/2, val + max(red_chi_sq)*0.01, 
                f'{val:.2f}', ha='center', va='bottom', fontsize=8)
    
    # Plot 2: Signal-to-Noise ratio
    snr_values = [p['signal_to_noise'] for p in peaks_info]
    bars2 = ax2.bar(peak_numbers, snr_values, alpha=0.7, color='lightgreen')
    ax2.set_xlabel('Peak Number')
    ax2.set_ylabel('Signal-to-Noise Ratio')
    ax2.set_title('Peak Signal-to-Noise Ratios')
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Peak widths (FWHM)
    fwhm_values = [p['fitted_fwhm'] for p in peaks_info]
    bars3 = ax3.bar(peak_numbers, fwhm_values, alpha=0.7, color='lightcoral')
    ax3.set_xlabel('Peak Number')
    ax3.set_ylabel('FWHM (min)')
    ax3.set_title('Peak Widths (Full Width at Half Maximum)')
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Total ion counts
    ion_counts = [p['gaussian_area_ions'] for p in peaks_info]
    bars4 = ax4.bar(peak_numbers, ion_counts, alpha=0.7, color='gold')
    ax4.set_xlabel('Peak Number')
    ax4.set_ylabel('Total Ion Counts')
    ax4.set_title('Total Ion Counts per Peak')
    ax4.grid(True, alpha=0.3)
    
    # Format y-axis for ion counts to avoid scientific notation
    ax4.ticklabel_format(style='plain', axis='y')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.show()

def plot_comprehensive_chromatogram(results, peaks_info, figsize=(15, 10)):
    """
    Plot the complete chromatogram with all Gaussian fits overlaid
    """
    time = results['corrected_arr'][0]
    signal = results['corrected_arr'][1]
    baseline_noise = results['baseline_noise']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)
    
    # Plot 1: Complete chromatogram with fits
    ax1.plot(time, signal, 'b-', alpha=0.7, label='Corrected Signal', linewidth=1)
    
    # Create a combined fit curve
    combined_fit = np.zeros_like(signal)
    time_indices = {t: i for i, t in enumerate(time)}
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(peaks_info)))
    
    for idx, (peak, color) in enumerate(zip(peaks_info, colors)):
        fit = peak['gaussian_fit']
        
        # Add this Gaussian to the combined fit
        for t, y in zip(fit['fit_x'], fit['fitted_curve_pA']):
            if t in time_indices:
                combined_fit[time_indices[t]] += y
        
        # Plot individual peak fits
        ax1.plot(fit['fit_x'], fit['fitted_curve_pA'], color=color, linewidth=2, 
                alpha=0.8, label=f'Peak {idx+1} Fit')
        
        # Mark peak centers
        ax1.axvline(x=peak['fitted_center'], color=color, linestyle='--', alpha=0.7)
        
        # Add peak labels
        ax1.text(peak['fitted_center'], peak['fitted_amplitude_pA'] * 1.1, 
                f'{idx+1}', ha='center', va='bottom', fontsize=10, color=color,
                bbox=dict(boxstyle="circle,pad=0.3", facecolor=color, alpha=0.2))
    
    ax1.plot(time, combined_fit, 'r-', alpha=0.5, linewidth=2, label='Combined Fit')
    ax1.set_title('Chromatogram with Gaussian Fits')
    ax1.set_xlabel('Retention Time (min)')
    ax1.set_ylabel('Signal (pA)')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Residuals of combined fit
    residuals = signal - combined_fit
    ax2.plot(time, residuals, 'k-', alpha=0.7, linewidth=1, label='Residuals')
    ax2.axhline(y=0, color='r', linestyle='-', alpha=0.5)
    ax2.fill_between(time, -3*baseline_noise, 3*baseline_noise, 
                    color='gray', alpha=0.2, label='±3σ noise band')
    
    # Calculate and display overall fit statistics
    overall_rms = np.sqrt(np.mean(residuals**2))
    within_noise = np.sum(np.abs(residuals) < 3 * baseline_noise) / len(residuals) * 100
    
    ax2.set_title(f'Residuals of Combined Fit (RMS: {overall_rms:.4f} pA, {within_noise:.1f}% within 3σ)')
    ax2.set_xlabel('Retention Time (min)')
    ax2.set_ylabel('Residual (pA)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

def plot_fit_quality_vs_snr(peaks_info, figsize=(10, 8)):
    """
    Plot fit quality vs signal-to-noise ratio
    """
    if not peaks_info:
        return
        
    fig, ax = plt.subplots(figsize=figsize)
    
    snr_values = [p['signal_to_noise'] for p in peaks_info]
    red_chi_sq = [p['reduced_chi_squared'] for p in peaks_info]
    ion_counts = [p['gaussian_area_ions'] for p in peaks_info]
    
    # Color by ion count
    sc = ax.scatter(snr_values, red_chi_sq, c=ion_counts, cmap='viridis', 
                   s=100, alpha=0.7, edgecolors='black')
    
    ax.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Ideal χ²/ν = 1')
    ax.axhline(y=3, color='orange', linestyle=':', alpha=0.7, label='Acceptable χ²/ν = 3')
    
    ax.set_xlabel('Signal-to-Noise Ratio')
    ax.set_ylabel('Reduced Chi-squared (χ²/ν)')
    ax.set_title('Fit Quality vs Signal-to-Noise Ratio')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add colorbar
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Total Ion Counts')
    
    plt.tight_layout()
    plt.show()

def plot_peak_detection_overview(results, peaks_info, figsize=(15, 12)):
    """
    Plot an overview of peak detection showing original signal, baseline, and detected peaks
    """
    time = results['original_chrom'][0]
    original_signal = results['original_chrom'][1]
    baseline_vals = results['baseline_arr'][1]
    corrected_signal = results['corrected_arr'][1]
    peak_mask = results['peak_mask']
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize)
    
    # Plot 1: Original signal with baseline
    ax1.plot(time, original_signal, 'b-', alpha=0.7, label='Original Signal', linewidth=1)
    ax1.plot(time, baseline_vals, 'r-', alpha=0.8, label='Baseline', linewidth=2)
    ax1.set_title('Original Signal with Baseline')
    ax1.set_ylabel('Signal (pA)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Baseline-corrected signal with peak regions
    ax2.plot(time, corrected_signal, 'g-', alpha=0.7, label='Corrected Signal', linewidth=1)
    
    # Shade peak regions
    peak_regions = find_contiguous_regions(peak_mask)
    for start, end in peak_regions:
        if end - start > 1:  # Only shade regions with multiple points
            ax2.fill_between(time[start:end], 0, corrected_signal[start:end], 
                           color='orange', alpha=0.3, label='Peak Region' if start == peak_regions[0][0] else "")
    
    # Mark detected peak centers
    for i, peak in enumerate(peaks_info):
        ax2.axvline(x=peak['center_time'], color='red', linestyle='--', alpha=0.7)
        ax2.text(peak['center_time'], peak['center_signal_pA'] * 1.05, f'{i+1}', 
                ha='center', va='bottom', fontsize=10, color='red',
                bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.8))
    
    ax2.set_title('Baseline-Corrected Signal with Detected Peaks')
    ax2.set_ylabel('Corrected Signal (pA)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Peak mask and noise threshold
    ax3.plot(time, peak_mask.astype(float), 'k-', alpha=0.8, label='Peak Mask', linewidth=1)
    ax3.axhline(y=0.5, color='r', linestyle='--', alpha=0.7, label='Threshold')
    ax3.set_yticks([0, 1])
    ax3.set_yticklabels(['Baseline', 'Peak'])
    ax3.set_title('Peak Detection Mask')
    ax3.set_xlabel('Retention Time (min)')
    ax3.set_ylabel('Region Type')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
def plot_peak_fit_quality(results, peaks_info, figsize=(15, 12)):
    """
    Plot to evaluate the quality of peak fits and show if full peaks are captured
    """
    time = results['corrected_arr'][0]
    signal = results['corrected_arr'][1]
    
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize)
    
    # Plot 1: Full chromatogram with fits
    ax1.plot(time, signal, 'b-', alpha=0.6, label='Signal', linewidth=1)
    
    colors = plt.cm.tab10(np.linspace(0, 1, len(peaks_info)))
    
    for idx, (peak, color) in enumerate(zip(peaks_info, colors)):
        fit = peak['gaussian_fit']
        
        # Plot the fit
        ax1.plot(fit['fit_x'], fit['fitted_curve_pA'], color=color, linewidth=2, 
                label=f'Peak {idx+1} Fit')
        
        # Plot the actual data in the fit region
        ax1.plot(fit['fit_x'], fit['actual_data_pA'], 'o', color=color, 
                markersize=3, alpha=0.7)
        
        # Mark the edges
        ax1.axvline(x=peak['left_edge_time'], color=color, linestyle=':', alpha=0.5)
        ax1.axvline(x=peak['right_edge_time'], color=color, linestyle=':', alpha=0.5)
        
        # Add peak number
        ax1.text(peak['fitted_center'], peak['fitted_amplitude_pA'] * 1.1, 
                f'{idx+1}', ha='center', va='bottom', fontsize=10, color=color)
    
    ax1.set_title('Peak Fits with Edge Detection')
    ax1.set_ylabel('Signal (pA)')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Residuals for each peak
    for idx, (peak, color) in enumerate(zip(peaks_info, colors)):
        fit = peak['gaussian_fit']
        residuals = fit['actual_data_pA'] - fit['fitted_curve_pA']
        ax2.plot(fit['fit_x'], residuals, 'o-', color=color, markersize=3, 
                alpha=0.7, label=f'Peak {idx+1}')
    
    ax2.axhline(y=0, color='k', linestyle='-', alpha=0.5)
    ax2.set_title('Fit Residuals by Peak')
    ax2.set_xlabel('Retention Time (min)')
    ax2.set_ylabel('Residual (pA)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Fit quality metrics
    peak_numbers = range(1, len(peaks_info) + 1)
    chi_sq_values = [p['reduced_chi_squared'] for p in peaks_info]
    rms_values = [p['rms_residuals_pA'] for p in peaks_info]
    
    ax3.bar(peak_numbers, chi_sq_values, alpha=0.6, label='Reduced χ²', color='skyblue')
    ax3.bar(peak_numbers, rms_values, alpha=0.6, label='RMS Residuals (pA)', color='lightcoral')
    
    ax3.axhline(y=1, color='red', linestyle='--', label='Ideal χ²/ν = 1')
    ax3.set_xlabel('Peak Number')
    ax3.set_title('Fit Quality Metrics')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

def diagnose_peak_fitting(results, peaks_info, specific_peaks=None):
    """
    Diagnostic function to check if peaks are being properly fitted
    """
    if specific_peaks is None:
        specific_peaks = range(len(peaks_info))
    
    print("=== PEAK FITTING DIAGNOSIS ===")
    
    for peak_idx in specific_peaks:
        if peak_idx >= len(peaks_info):
            continue
            
        peak = peaks_info[peak_idx]
        fit = peak['gaussian_fit']
        
        print(f"\nPeak {peak_idx + 1}:")
        print(f"  Retention time: {peak['fitted_center']:.3f} min")
        print(f"  Fit range: {peak['left_edge_time']:.3f} - {peak['right_edge_time']:.3f} min")
        print(f"  Fit width: {peak['right_edge_time'] - peak['left_edge_time']:.3f} min")
        print(f"  Data points in fit: {len(fit['fit_x'])}")
        print(f"  Amplitude: {peak['fitted_amplitude_pA']:.2f} pA")
        print(f"  FWHM: {peak['fitted_fwhm']:.3f} min")
        print(f"  Reduced χ²: {peak['reduced_chi_squared']:.2f}")
        print(f"  RMS residuals: {peak['rms_residuals_pA']:.4f} pA")
        
        # Check if the fit captures the peak properly
        actual_max = np.max(fit['actual_data_pA'])
        fitted_max = np.max(fit['fitted_curve_pA'])
        ratio = fitted_max / actual_max if actual_max > 0 else 0
        
        if ratio < 0.9:
            print(f"  WARNING: Fit captures only {ratio:.1%} of peak height!")
        elif ratio > 1.1:
            print(f"  WARNING: Fit overestimates peak height by {(ratio-1):.1%}!")
        else:
            print(f"  Fit quality: Good ({ratio:.1%} of actual height)")


# Update the main processing to include debugging and manual inspection
if __name__ == "__main__":
    # Load and process chromatogram
    d1 = cc.Day(r"C:\AutoGCData\RB", "20250802")
    
    for i, c in enumerate(d1.chromatograms[1:2]):
        print(f"Processing chromatogram {i+1}")
        
        chrom_data = c.back.chromatogram
        if chrom_data is None:
            continue
            
        # Perform baseline correction
        results = baseline(chrom_data)
        
        # Plot baseline results first
        plot_baseline_results(results, f"Chromatogram {i+1} - Baseline")
        
        # Manual inspection to see what peaks should be there
        print(f"\n=== MANUAL PEAK INSPECTION ===")
        manual_regions = manual_peak_inspection(results, min_snr=5.0)
        print(f"Found {len(manual_regions)} potential peak regions manually")
        
        # Use the improved peak finding
        print(f"\n=== AUTOMATED PEAK DETECTION ===")
        peaks_info, all_peak_indices = enhanced_find_peaks_complete(
            results, 
            min_snr=3.0,
            min_prominence=0.005,
            acquisition_time=0.2
        )
        
        if peaks_info:
            print(f"Automated detection found {len(peaks_info)} peaks")
            
            # Create a modified version that returns all detected indices for debugging
            # For debugging purposes, let's also get all potential peaks
            time = results['corrected_arr'][0]
            signal = results['corrected_arr'][1]
            baseline_noise = results['baseline_noise']
            
            # Smooth the signal
            signal_smooth = savgol_filter(signal, 21, 3) if len(signal) > 21 else signal
            
            # Get all potential peaks for debugging
            all_peak_indices = find_all_significant_peaks(
                time, signal_smooth, baseline_noise, min_snr=10.0, min_prominence=0.01
            )
            
            # Plot debug information
            print(f"\n=== PEAK DETECTION DEBUG ===")
            plot_peak_detection_debug(results, all_peak_indices, peaks_info)
            
            # Compare manual vs automated detection
            print(f"\n=== COMPARISON: MANUAL vs AUTOMATED ===")
            manual_centers = []
            for start, end in manual_regions:
                if end - start >= 5:  # Only consider substantial regions
                    region_signal = signal[start:end]
                    peak_idx = start + np.argmax(region_signal)
                    manual_centers.append(time[peak_idx])
            
            auto_centers = [peak['fitted_center'] for peak in peaks_info]
            
            print(f"Manual detection found peaks at: {[f'{x:.3f}' for x in manual_centers]} min")
            print(f"Automated detection found peaks at: {[f'{x:.3f}' for x in auto_centers]} min")
            
            # Check for missed peaks
            missed_peaks = []
            for manual_rt in manual_centers:
                # Check if any automated peak is close to this manual peak
                closest_match = min(auto_centers, key=lambda x: abs(x - manual_rt)) if auto_centers else None
                if closest_match is None or abs(closest_match - manual_rt) > 0.05:  # 0.05 min tolerance
                    missed_peaks.append(manual_rt)
            
            if missed_peaks:
                print(f"WARNING: Potential missed peaks at: {[f'{x:.3f}' for x in missed_peaks]} min")
            else:
                print("All manual peaks were detected automatically")
            
            # Run diagnostics on the fitting
            print(f"\n=== PEAK FITTING DIAGNOSTICS ===")
            diagnose_peak_fitting(results, peaks_info)
            
            # Plot the Gaussian fits
            print(f"\n=== GAUSSIAN FITS ===")
            plot_gaussian_fits(results, peaks_info)
            
            # Plot comprehensive overview
            plot_comprehensive_chromatogram(results, peaks_info)
            
            # Plot fit quality vs SNR
            plot_fit_quality_vs_snr(peaks_info)
            
            # Print detailed information
            print(f"\n=== DETAILED PEAK INFORMATION ===")
            print_detailed_peak_info(peaks_info)
            
            # Summary statistics
            print(f"\n=== SUMMARY STATISTICS ===")
            print(f"Total peaks: {len(peaks_info)}")
            print(f"Retention time range: {min(auto_centers):.2f} - {max(auto_centers):.2f} min")
            print(f"Largest peak: {max(p['fitted_amplitude_pA'] for p in peaks_info):.2f} pA")
            print(f"Average SNR: {np.mean([p['signal_to_noise'] for p in peaks_info]):.1f}")
            print(f"Average reduced χ²: {np.mean([p['reduced_chi_squared'] for p in peaks_info]):.2f}")
            print(f"Total ion count: {sum(p['gaussian_area_ions'] for p in peaks_info):.0f} ions")
            
        else:
            print("No peaks found in this chromatogram")
            print("This suggests either:")
            print("1. The signal is very noisy")
            print("2. The baseline correction failed")
            print("3. The peak detection parameters are too strict")
            
            # Plot the signal to help diagnose
            plt.figure(figsize=(12, 6))
            time = results['corrected_arr'][0]
            signal = results['corrected_arr'][1]
            plt.plot(time, signal, 'b-', label='Corrected Signal')
            plt.axhline(y=results['noise_threshold'], color='r', linestyle='--', 
                       label=f'Noise Threshold ({results["noise_threshold"]:.3f} pA)')
            plt.xlabel('Retention Time (min)')
            plt.ylabel('Signal (pA)')
            plt.title(f'Chromatogram {i+1} - No Peaks Detected')
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.show()