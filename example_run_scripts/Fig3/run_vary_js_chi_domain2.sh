#!/bin/bash

# Define the range of arguments
js1=(-0.1 -0.2 0)
js1_prev=(basic_structure -0.1 -0.1)

# Loop through the arguments of js1 and js2 and bash run_vary_chi.sh
for i in {0..2}
do
    nohup julia --threads 4 ForHeatmap_domain2.jl ${js1_prev[$i]} ${js1[$i]} 3.0 3.0 &>curr.log
    echo "Running run_vary_chi_domain.sh for js=${js1[$i]}"
    nohup bash run_vary_chi_domain2.sh ${js1[$i]} ${js1[$i]} &>vary_curr_chi.log
done

