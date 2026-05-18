#!/usr/bin/env python3 

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt

def read_pgm_p2(filename):
    """Reads a PGM file in text mode (P2) ignoring comments."""
    words = []
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.split('#')[0]
                words.extend(line.split())
    except FileNotFoundError:
        print(f"Error: The file '{filename}' does not exist.")
        sys.exit(1)

    if not words:
        raise ValueError("The file is empty.")
    if words[0] != 'P2':
        raise ValueError("The file is not in PGM P2 format (text mode).")

    width = int(words[1])
    height = int(words[2])
    max_val = int(words[3])

    pixels = np.array([int(w) for w in words[4:]], dtype=np.int32).reshape((height, width))
    return pixels, width, height, max_val

def extract_points(pixels, max_val):
    """Detects point coordinates assuming they are the minority of pixels."""
    unique, counts = np.unique(pixels, return_counts=True)
    background_val = unique[np.argmax(counts)]
    point_indices = np.argwhere(pixels != background_val)
    points = point_indices[:, [1, 0]].astype(float)
    return points

def main():
    parser = argparse.ArgumentParser(
        description="Analyzes hyperuniformity by calculating fluctuation exponents."
    )
    parser.add_argument("filename", help="Path to the PGM file in text mode")
    parser.add_argument("--samples", type=int, default=1000, 
                        help="Number of sampling windows per radius (default: 1000)")
    args = parser.parse_args()

    # 1. Read image
    print(f"Reading file: {args.filename}...")
    pixels, width, height, max_val = read_pgm_p2(args.filename)
    L = min(width, height)
    
    center_x, center_y = width / 2.0, height / 2.0
    print(f"Image detected: {width}x{height}. Center at: ({center_x}, {center_y})")

    # 2. Extract points
    points = extract_points(pixels, max_val)
    print(f"{len(points)} points detected in the distribution.")
    print(f"Density: {len(points)/(width * height):.4f}")

    # 3. Configure radii R (we avoid extremely small radii for a better fit)
    max_R = L / 2.0
    radii = np.linspace(5, max_R, 30) 
    
    var_N = []
    var_rho = []

    print("Calculating density fluctuations (Sampling with PBC)...")
    for R in radii:
        counts = []
        cx = np.random.uniform(0, L, args.samples)
        cy = np.random.uniform(0, L, args.samples)
        
        for x_c, y_c in zip(cx, cy):
            dx = np.abs(points[:, 0] - x_c)
            dy = np.abs(points[:, 1] - y_c)
            
            dx = np.minimum(dx, L - dx)
            dy = np.minimum(dy, L - dy)
            
            dist_sq = dx**2 + dy**2
            num_inside = np.sum(dist_sq <= R**2)
            counts.append(num_inside)
            
        counts = np.array(counts)
        v_N = np.var(counts)
        v_rho = v_N / ((np.pi * R**2) ** 2)
        
        var_N.append(v_N)
        var_rho.append(v_rho)

    radii = np.array(radii)
    var_N = np.array(var_N)
    var_rho = np.array(var_rho)

    # 4. Exponential Fit Calculation (Linear regression in Log-Log scale)
    log_R = np.log(radii)
    
    # Fit for var_N: log(var_N) = alpha * log(R) + C
    alpha, intercept_N = np.polyfit(log_R, np.log(var_N), 1)
    fit_N = np.exp(intercept_N) * (radii ** alpha)

    # Fit for var_rho: log(var_rho) = beta * log(R) + C'
    beta, intercept_rho = np.polyfit(log_R, np.log(var_rho), 1)
    fit_rho = np.exp(intercept_rho) * (radii ** beta)

    # --- CONSOLE OUTPUT ---
    print("\n" + "="*60)
    print("               EXPERIMENTAL FIT RESULTS")
    print("="*60)
    print(f"-> Exponent of var_N (alpha)  : {alpha:.3f}")
    print(f"   [Poisson = 2.00 | Hyperuniform < 2.00 (e.g. ~ 1.00)]")
    print(f"-> Exponent of var_rho (beta) : {beta:.3f}")
    print(f"   [Poisson = -2.00 | Hyperuniform < -2.00 (e.g. ~ -3.00)]")
    print("="*60 + "\n")

    # 5. Plot the results
    plt.figure(figsize=(13, 5.5))

    # Plot 1: Variance of Number of Points
    plt.subplot(1, 2, 1)
    plt.loglog(radii, var_N, 'o', color='blue', label='Measured data')
    plt.loglog(radii, fit_N, 'b-', label=fr'Experimental Fit ($\alpha$ = {alpha:.2f})')
    plt.loglog(radii, (radii**2) * (var_N[0] / (radii[0]**2)), '--', color='gray', label=r'Theoretical Poisson ($R^2$)')
    plt.loglog(radii, radii * (var_N[0] / radii[0]), ':', color='red', label=r'Theoretical Hyperuniform ($R^1$)')
    plt.xlabel('Radius (R)')
    plt.ylabel(r'Variance $\sigma_N^2$')
    plt.title('Point Number Variance Growth')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    # Plot 2: Density Variance
    plt.subplot(1, 2, 2)
    plt.loglog(radii, var_rho, 's', color='purple', label='Measured data')
    plt.loglog(radii, fit_rho, 'r-', label=fr'Experimental Fit ($\beta$ = {beta:.2f})')
    plt.loglog(radii, (radii**-2) * (var_rho[0] / (radii[0]**-2)), '--', color='gray', label=r'Theoretical Poisson ($R^{-2}$)')
    plt.loglog(radii, (radii**-3) * (var_rho[0] / (radii[0]**-3)), ':', color='red', label=r'Theoretical Hyperuniform ($R^{-3}$)')
    plt.xlabel('Radius (R)')
    plt.ylabel(r'Density Variance $\sigma_{\rho}^2$')
    plt.title('Density Variance Decay')
    plt.legend()
    plt.grid(True, which="both", ls="--")

    plt.tight_layout()
    print("Showing interactive plots...")
    plt.show()

if __name__ == "__main__":
    main()
