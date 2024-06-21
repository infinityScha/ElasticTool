#!/bin/bash

# Define the range of arguments
js1=(0 0.2 0.3 0.4 0.5 0.6 -0.1 -0.2)
js1_prev=(0.1 0.1 0.2 0.3 0.4 0.5 0 -0.1)

# Loop through the arguments of js1 and js1_prev and run ForHeatmap.jl
for i in {0..7}
do
    echo "Running ForHeatmap_domain.jl with js1=${js1[$i]} and js1_prev=${js1_prev[$i]}"
    nohup julia --threads 8 ForHeatmap_domain.jl ${js1_prev[$i]} ${js1[$i]} 3.0 3.0 &>curr.log
done
