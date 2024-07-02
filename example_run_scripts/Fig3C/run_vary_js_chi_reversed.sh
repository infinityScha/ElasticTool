#!/bin/bash

# Define the range of arguments
js1=(0.25 0.2 0.15 0.1 0.05 0 -0.05)
js1_prev=(0.3 0.25 0.2 0.15 0.1 0.05 0)

# Loop through the arguments of js1 and js2 and bash run_vary_chi.sh
for i in {0..8}
do   
    # first run for chi=0
    echo "Running ForHeatmap_reversed.jl with js=${js1[$i]} and chi=3.0"
    nohup julia --threads=2 ForHeatmap_reversed.jl ${js1_prev[$i]} ${js1[$i]} 3.0 3.0 &>vary_curr_chi_rev.log
    echo "Running run_vary_chi_reversed.sh for js=${js1[$i]}"
    nohup bash run_vary_chi_reversed.sh ${js1[$i]} ${js1[$i]} &>vary_curr_chi_rev.log
done

