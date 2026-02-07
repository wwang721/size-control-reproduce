# Supplementary code & data - Controlling tissue size by active fracture

[![PRE Paper](https://img.shields.io/badge/DOI-10.1103%2Fdk15--hwzg-blue)](https://doi.org/10.1103/dk15-hwzg)
[![arXiv:2503.03126](https://img.shields.io/badge/arXiv-2503.03126-grey.svg?colorB=a42c25&logo=arxiv)](https://doi.org/10.48550/arXiv.2503.03126)
[![CI](https://github.com/wwang721/size-control-reproduce/actions/workflows/ci.yml/badge.svg)](https://github.com/wwang721/size-control-reproduce/actions/workflows/ci.yml)

These are the code and data required to reproduce the results in the paper:

- ***Controlling tissue size by active fracture***, Wei Wang (汪巍) and Brian A. Camley, [**Phys. Rev. E (2026)**](https://doi.org/10.1103/dk15-hwzg).

Preprint version available on **arXiv**: [arXiv:2503.03126](https://doi.org/10.48550/arXiv.2503.03126).

> [!CAUTION]
> This code is only intended to help in reproducing our results and explaining our methods fully - attempting to generalize without full understanding may lead to unpredictable issues or errors.


## Analytical results

Detailed analytical derivations of the generating functions are provided in the **Mathematica** notebook [generating_functions.nb](/generating_functions.nb).


## Requirements

The code was run with **Python 3.11.11** and the following packages:

| Package    | Version | Usage                                                                           |
| :--------- | :-----: | :------------------------------------------------------------------------------ |
| numpy      | 2.1.3   | Numerical computations                                                          |
| scipy      | 1.15.3  | Miscellaneous scientific functions                                              |
| matplotlib | 3.10.0  | Plotting and visualization                                                      |
| tqdm       | 4.67.1  | Progress bars during calculations                                               |
| mpi4py     | 4.0.3   | Parallel processing using MPI                                                   |
| lifelines  | 0.30.0  | Survival analysis                                                               |
| imageio    | 2.37.9  | Image and video reading/writing utilities (optional; only for video generation) |


## Usage

### Figure 2

For panel (b), the AOUP simulation code is located at [variances/NAOUPs.py](/variances/NAOUPs.py). Jobs are submitted using [script.slurm](/variances/script.slurm). The output data are written to [variances/data](/variances/data) and contain the final particle positions for each run. The figure is generated using [variance/plot.py](/variances/plot.py).

> [!NOTE]
> For small values of the persistence time $\tau$, the first few data points may exhibit increased noise because $\tau$ approaches the default integration time step $\Delta t$. In this regime, reducing the time step (using line 16 and commenting out line 13 in [NAOUPs.py](/variances/NAOUPs.py) for small $\tau$) and averaging over multiple runs is necessary to obtain reliable results.

For panel (c), the code used to compute the mean first-passage time is located in the directory [MFPT](/MFPT). Each subdirectory corresponds to varying a single parameter. Submit the `script.slurm` file within each subdirectory to run the simulations. The resulting escape times $\tau_\mathrm{esc}$ are saved as `res.npy` files. Running [MFPT/plot.py](/MFPT/plot.py) aggregates these results and generates the figure.


### Figure 4

Most curves in this figure are theoretical and do not rely on simulations.
We include only the rfKMC simulation code used to generate the numerical data points in panels (a) and (e).
The directory [distributions](/distributions) contains scripts for generating cluster-size distributions for the three models studied in the paper, with output saved to [distributions/data](/distributions/data). Running [distributions/plot.py](/distributions/plot.py) produces panel (a). For panel (e), the corresponding simulation scripts are located in [distributions/min_N](/distributions/min_N/).


### Figure 6

To perform the survival analysis for the 1D rfKMC simulations, run [survival/lineage_rupture.py](/survival/lineage_rupture.py) to generate rupture times, which are saved in [survival/data](/survival/data) (run different $k_b$ values manually). Subsequently, run [survival/plot.py](/survival/plot.py) to generate panel (a). The directory [survival/change_kd](/survival/change_kd) contains scripts in which the divisioin rate $k_d$ is varied (run different values of $k_d$ manually), corresponding to panel (b).


### Figure 7

The 2D AOUP simulations are located in [2D_simulations](/2D_simulations) directory.
[main.py](/2D_simulations/main.py) runs the main simulations and writes snapshot images to the [frames](/2D_simulations/frames/) directory. These snapshots can then be assembled into videos by running [video_generate.py](/2D_simulations/video_generate.py).

![2D simulation](/2D_simulations/2d_sim.gif)

Submit [2D_simulations/script.slurm](/2D_simulations/script.slurm) to run [scan_va.py](/2D_simulations/scan_va.py), which scans the parameter space over the active velocity $v_a$ (denoted $v_\textrm{rms}$ in the paper) and saves the results to [2D_simulations/data](/2D_simulations/data). The scripts in [2D_simulations/analysis](/2D_simulations/analysis) process these data to generate panels in Fig. 7, including fracture rates, histograms, cluster-size distributions, and survival probabilities.


### Other figures

All remaining figures are either purely theoretical or minor variants of the figures described above and can be generated following analogous procedures.


## More information

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18485014.svg)](https://doi.org/10.5281/zenodo.18485014)

- This repository is also archived and cross-listed on [Zenodo](https://doi.org/10.5281/zenodo.18485014).
- The authors acknowledge support from the National Institutes of Health (NIH) Grant No. R35GM142847.
- This work was carried out at the Advanced Research Computing at Hopkins (ARCH) core facility, which is supported by the National Science Foundation (NSF) Grant No. OAC1920103.


## License

This project is licensed under the [MIT License](/LICENSE), which essentially means you have the freedom to use, modify, and distribute this code for almost any purpose.
