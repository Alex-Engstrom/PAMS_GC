# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 16:06:35 2025

@author: aengstrom
"""
from chromclass import Chromatogram
from pybaselines import Baseline
import matplotlib.pyplot as plt

def baseline(chromatogram: Chromatogram, lam=1e5, p=0.01):
    """
    Performs baseline correction using pybaselines.asls.
    Returns both the estimated baseline and the corrected signal.
    """
    baseline_fitter = Baseline()
    
    # Use ASLS for robust and smooth baseline estimation
    base, _ = baseline_fitter.asls(
        chromatogram.chromatogram[1],
        lam=lam,
        p=p
    )

    corrected = chromatogram.chromatogram[1] - base
    chromatogram.chromatogram[1] = corrected
    return base, corrected

if __name__ == "__main__":
    # Load chromatogram from separate class
    c1 = Chromatogram(
        "C:/Users/aengstrom/Downloads/rbch06cdat-Front Signal.cdf", "cdf"
    )

    rt = c1.chromatogram[0]
    raw_signal = c1.chromatogram[1].copy()

    # Perform baseline correction
    base, corrected = baseline(c1, lam=1e5, p=0.0001)

    # Plot results
    plt.figure(figsize=(14, 7))
    plt.plot(rt, raw_signal, label="Raw Chromatogram", color="blue", alpha=0.6)
    plt.plot(rt, base, label="Estimated Baseline", color="orange", linewidth=2)
    plt.plot(rt, corrected, label="Corrected Chromatogram", color="green", alpha=0.8)

    plt.xlabel("Retention Time (s)")
    plt.ylabel("Signal (FID)")
    plt.title("Baseline Correction Using ASLS")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

