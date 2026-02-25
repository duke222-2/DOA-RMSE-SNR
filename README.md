Basic DOA Estimation: Performance Analysis (RMSE vs. SNR)
This repository contains the Matlab implementation of fundamental Direction of Arrival (DOA) estimation algorithms and their statistical performance analysis.

🚀 Overview
The project focuses on evaluating the RMSE (Root Mean Square Error) performance against SNR (Signal-to-Noise Ratio) for classical array signal processing algorithms.

Algorithms Included:
MUSIC (Multiple Signal Classification)

ESPRIT (Estimation of Signal Parameters via Rotational Invariance Techniques)

CRB Calculation (Cramér-Rao Bound as a theoretical benchmark)

📂 Key Files
DOA1DSNR_RMSE.m: Main script to generate RMSE vs. SNR curves.

fun_MUSIC_1D.m: Implementation of the 1D MUSIC algorithm.

fun_esprit_1D.m: Implementation of the 1D ESPRIT algorithm.

DOA2D_CRB.m: Script for calculating the theoretical Cramer-Rao Bound.

🛠 Usage
Open Matlab and add the folder to your path.

Run DOA1DSNR_RMSE.m to start the simulation.

The results will be plotted automatically, comparing the algorithms with the CRB.
