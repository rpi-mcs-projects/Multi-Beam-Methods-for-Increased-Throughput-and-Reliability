# Multi-Beam-Methods-for-Increased-Throughput-and-Reliability

## dual_beam_scanning_plotting.m
This MATLAB file can be run with the phased array and 5G toolboxes and the two STL 3D model files provided in the repository. This script simulates the multipath profile of the 3D environment and searches for the two most constructive paths. The multipath profiles are visualized. 

A phased array is steered towards the two optimal constructive channels and the response is plotted in 2D and 3D. 

## n_beam_comparison.m
This script extends `dual_beam_scanning_plotting.m` to a higher number of beams, and sweeps across the number of beams. Then the SNR is plotted across number of beams.

For both the above scripts, the variable `scenefile` can be set to either `env.stl` for the indoor environment or `env_outdoor.stl`

