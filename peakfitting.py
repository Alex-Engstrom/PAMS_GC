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



#%% Plotting functions


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