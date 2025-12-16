# Multi-Beam-Methods-for-Increased-Throughput-and-Reliability

This is the final project code for Jamie Draper and Edgar Muniz for Modern Communication Systems Fall 2025

## dual_beam_scanning_plotting.m
This MATLAB file can be run with the phased array and 5G toolboxes and the two STL 3D model files provided in the repository. This script simulates the multipath profile of the 3D environment and searches for the two most constructive paths. The multipath profiles are visualized. 

A phased array is steered towards the two optimal constructive channels and the response is plotted in 2D and 3D. 

## n_beam_comparison.m
This script extends `dual_beam_scanning_plotting.m` to a higher number of beams, and sweeps across the number of beams. Then the SNR is plotted across number of beams.

For both the above scripts, the variable `scenefile` can be set to either `"env.stl"` for the indoor environment or `"env_outdoor.stl"`


## tracking.m: 

Orignal test simulation. This was where the orignal tracking was implemented. This simulation contains 2 paths (1 LOS, 1 Reflected path). Output plots show the first algorithm implementation's angle calculations/tracking. Also compares SNR of multi-beam (2 beam) tracking and single beam tracking

## py_mm_test.py 

Implements a multipath environment (max 4 constructive paths). In additon, it implements the power loss inference target tracking algorithm. It also incorporates the optimal beam forming algorithm. A general path loss of 78.29 db is applied to the paths. Output plots provide a comparison of single beam and multibeam tracking/SNR with simulated blockages. Also provides a polar plot showing the optimal beam forming.

Running python 3.14.2
required packages: numpy, matplotlib, pandas
