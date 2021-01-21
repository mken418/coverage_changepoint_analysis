# coverage_changepoint_analysis
This code is designed to run on the output of Sentieon's coverageMetrics algorithm for a structural variant, which generates per-base coverage data for an input bed range. In order for these scripts to work correctly, organize the output of Sentieon's coverageMetrics into 1 file per sample and place them into a single directory per structural variant. 

1. Use the python wrapper, *calculate_signals_from_input_loci.py*, to call the matlab script, *SignalDetection_CL_args.m*. Matlab will need to installed, along with the matlab module for python. 
* will output a single file containing the changepoint detection analysis, one sample per line, for each input sample's coverageMetrics output in the input directory.

2. Plot the coverage data along with output from the changepoint detection analysis with *plot_coverage.py*
* will output a plot with a maximum of 25 individuals with bolded names of individuals where the *SignalDetection_CL_args.m* determined there was a variant based on change point detection. 
