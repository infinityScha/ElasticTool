#!/bin/bash

# Define the range of arguments
js1=(-0.05 0)
js1_prev=(-0.05 -0.05)

# Loop through the arguments of js1 and js2 and bash run_vary_chi.sh
for i in {0..1}
do
    # first run for chi=0
    echo "Running ForHeatmap.jl with js=${js1[$i]} and chi=0"
    nohup julia --threads=2 ForHeatmap_domain2.jl ${js1_prev[$i]} ${js1[$i]} 2.6 2.6 &>vary_curr_chi.log
    echo "Running run_vary_chi.sh for js=${js1[$i]}"
    nohup bash run_vary_chi_domain2.sh ${js1[$i]} ${js1[$i]} &>vary_curr_chi.log
done

