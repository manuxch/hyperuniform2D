# Hyperuniform2D

Generator for hyperuniform points in a 2D square grid

This repository contains a C++20 program that generates **binary hyperuniform distributions** on an $N \times N$ square lattice. The optimization is carried out via **Metropolis simulated annealing**, minimizing an energy based on the structure factor $S(\mathbf{k})$.

The output is a configuration of occupied/empty sites that suppresses large-scale density fluctuations — a fundamental property of disordered systems with hidden order, relevant in photonics, granular materials, and biology.

## Requirements

To build and run the simulation, you will need:
- A **C++20** compatible compiler (e.g., GCC 10+, Clang 10+, MSVC 19.29+)
- **CMake** 3.20 or newer

## Compilation

The project uses CMake as its build system. It is highly recommended to compile in `Release` mode to ensure maximum performance during the simulated annealing process.

```bash
# Generate build files
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# Compile the project
cmake --build build --parallel
```

## Usage

After compilation, the executable `hyper2D` will be located in the `build/` directory. The program takes a single optional argument: the path to the parameter configuration file. If no argument is provided, it defaults to `params.in`.

```bash
./build/hyper2D params.in
```

### Configuration File (`params.in`)

The parameter file uses a simple `key: value` format. Lines starting with `#` are considered comments. 

Example:
```yaml
# Grid size
N: 64
# Occupation fraction
density: 0.4
# Exclusion radius (in units of 2π/N)
K: 4.0
# Simulated annealing parameters
T_initial: 10.0
T_final: 0.01
cooling_rate: 0.999
steps_per_temp: 1000
# Random seed
seed: 1729
# Number of threads (0 = auto-detect, 1 = sequential strict)
num_threads: 0
# Output files
output_file: hyperuniform.pgm
annealing_file: annealing.dat
sk_file: sk.dat
```


### Multithreading

The code features a highly optimized, native C++20 thread pool using `std::barrier` to accelerate the evaluation of the structure factor modes.

- Set `num_threads: 0` to automatically detect and use all available hardware threads.
- Set `num_threads: 1` to strictly enforce sequential execution.
- **Smart Fallback**: Even if multiple threads are requested, the program will automatically fallback to sequential execution if the number of evaluated modes is small ($< 100$). This ensures maximum performance by avoiding synchronization overhead on small workloads.

### Output Files

The simulation generates three main output files (names can be configured in `params.in`):

1. **`hyperuniform.pgm`**: The final spatial configuration saved as a plain PGM (Portable GrayMap) image. The values represent grid occupation (0 or 1).
2. **`annealing.dat`**: Monitoring data from the simulated annealing process. It contains three columns: Temperature (`T`), Energy (`E`), and Acceptance Rate (`acc/rate`). This file is useful for tuning the cooling schedule.
3. **`sk.dat`**: The structure factor $S(\mathbf{k})$ evaluated at $k_y = 0$ for the final configuration. It contains two columns: $k_x$ and $S(k_x, 0)$.

#### Note 

The PGM image is stored in plain text format (magic number P2). Occupied sites (1) are shown as white pixels, and empty sites (0) as black pixels. If you prefer the opposite coloring (black for occupied, white for empty), you can invert the image using ImageMagick:

    magick input.pgm -negate -compress none output.pgm

The `-compress none` option keeps the output in human‑readable PGM text format.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## What is a hyperuniform distribution?

A point pattern (or binary occupation) is said to be **hyperuniform** if the variance of the number of particles within an observation window grows more slowly than the window volume as the window size increases. Equivalently, its **structure factor** $S(\mathbf{k})$ tends to zero as the wave vector $\mathbf{k}$ approaches zero:

$$
\lim_{|\mathbf{k}|\to 0} S(\mathbf{k}) = 0.
$$

A perfectly periodic system (like a crystal) has $S(\mathbf{k})=0$ away from the Bragg peaks, so a crystal is hyperuniform. However, hyperuniformity can also occur in **disordered** systems, where particles appear randomly distributed but possess long-range correlations that cancel large-scale fluctuations.

This program generates a special type of disordered hyperuniform system known as **“stealthy” hyperuniform**, in which the structure factor is **exactly zero** for all wave vectors with $|\mathbf{k}| \le K$, where $K$ is a user-defined exclusion radius.

## Mathematical foundation

### System representation

We consider an $N \times N$ square lattice with periodic boundary conditions. Each site $\mathbf{r} = (x, y)$ can be occupied (`1`) or empty (`0`). The occupation fraction is $\rho = M/N^2$, where $M$ is the total number of occupied sites.

The discrete Fourier transform of the occupation is:

$$
\hat{\rho}(\mathbf{k}) = \sum_{\mathbf{r}} s(\mathbf{r})  e^{-i \mathbf{k} \cdot \mathbf{r}},
\qquad \mathbf{k} = \frac{2\pi}{N}(k_x, k_y), \quad k_x,k_y \in \mathbb{Z}.
$$

The structure factor is defined as 

$$
S(\mathbf{k}) = \frac{1}{N^2} |\hat{\rho}(\mathbf{k})|^2
$$

### Stealthy hyperuniformity condition

To achieve $S(\mathbf{k}) \to 0$ as $\mathbf{k} \to 0$, we explicitly constrain all modes with $0 < |\mathbf{k}| \leq K$ (in units of $2\pi/N$). The condition $S(\mathbf{k}) = 0$ in that region is equivalent to requiring the Fourier coefficients 

$$
\hat{\rho}(\mathbf{k})
$$

to vanish.


### Optimization problem

We define an **energy** that penalizes the presence of Fourier modes within the selected mask:

$$
E = \sum_{\substack{\mathbf{k} \neq \mathbf{0} \\ |\mathbf{k}| \le K}} |\hat{\rho}(\mathbf{k})|^2.
$$

A configuration with $E = 0$ is perfectly hyperuniform within the suppressed range (up to finite-size effects). Because the search space is discrete and huge, we resort to **simulated annealing** to find near-zero-energy configurations.

## Algorithm

1. **Initialization**: $M$ particles are placed randomly on the lattice, respecting the target density.

2. **Mask construction**: All wave vectors $\mathbf{k} \neq \mathbf{0}$ with magnitude less than or equal to \(K\) are selected. For each one, phase factors are precomputed.

3. **Initial Fourier transform**: 

$$
\hat{\rho}(\mathbf{k})
$$ 

is evaluated for all modes in the mask.

4. **Simulated annealing**:
   - At each step, a random occupied site and a random empty site are chosen, and an exchange is proposed.
   - The change in energy $\Delta E$ is computed using an incremental update of the Fourier coefficients.
   - The swap is accepted according to the Metropolis criterion:

$$
P_{\text{accept}} = \min\left(1, e^{-\Delta E / T}\right),
$$
   
  where $T$ is the current temperature.
   - The temperature is gradually reduced by a user-specified cooling factor.

5. **Result**: After a sufficient number of steps, a binary configuration with near-zero energy is obtained — the long-wavelength modes are suppressed and the system exhibits hyperuniformity.

The program saves the final configuration as a PGM (grayscale) image and prints the residual structure factor along the $k_x$ axis to verify mode suppression.

## Parameter file

The simulation is controlled via an external text file. Its format is `key: value`, and lines starting with `#` are treated as comments.

### Available keys

| Key               | Type   | Description                                                                 | Default value       |
|-------------------|--------|-----------------------------------------------------------------------------|---------------------|
| `N`               | integer| Lattice size (N×N)                                                          | 64                  |
| `density`        | real   | Occupation fraction (between 0 and 1)                                       | 0.4                 |
| `K`               | real   | Exclusion radius in units of \(2\pi/N\) (defines the mask)                  | 4.0                 |
| `T_initial`       | real   | Initial temperature for simulated annealing                                 | 10.0                |
| `T_final`         | real   | Final temperature (stopping criterion)                                      | 0.01                |
| `cooling_rate`    | real   | Cooling multiplier (0 < factor < 1)                                         | 0.999               |
| `steps_per_temp`  | integer| Number of swap attempts at each temperature level                           | 1000                |
| `seed`         | integer| Random seed (for reproducibility)                                           | 42                  |
| `num_threads`  | integer| Number of threads for concurrent execution                           | 0                |
| `output_file`  | string | Output PGM file name                                                        | hiperuniform.pgm   |
| `annealing_file`  | string | Annealing algorithm progression                                                        | annealing.dat   |
| `sk_file`  | string | Resulting S(k)                                                         | sk.dat   |

### Example `params.in`

### References

- Torquato, S., & Stillinger, F. H. (2003). Local density fluctuations, hyperuniformity, and order metrics. Physical Review E, 68(4), 041113. [https://doi.org/10.1103/PhysRevE.68.041113](https://doi.org/10.1103/PhysRevE.68.041113).
- Robert D. Batten, Frank H. Stillinger, Salvatore Torquato; Classical disordered ground states: Super-ideal gases and stealth and equi-luminous materials. J. Appl. Phys. 1 August 2008; 104 (3): 033504. [https://doi.org/10.1063/1.2961314](https://doi.org/10.1063/1.2961314).
- Chase E Zachary and Salvatore Torquato J. Stat. Mech. (2009) P12015. [https://doi.org/10.1088/1742-5468/2009/12/P12015](https://doi.org/10.1088/1742-5468/2009/12/P12015).
