#!/bin/bash

# Define the range of arguments
js1=(0 0.1 0.2 0.3 0.4 0.5 0.6 -0.1 -0.2)
js1_prev=(basic_structure 9 0.1 0.2 0.3 0.4 0.5 0 -0.1)

# Loop through the arguments of js1 and js1_prev and run ForHeatmap.jl
for i in {7..8}
do
    echo "Running ForHeatmap.jl with js1=${js1[$i]} and js1_prev=${js1_prev[$i]}"
    nohup julia --threads 8 ForHeatmap.jl ${js1_prev[$i]} ${js1[$i]} 0 0 &>curr.log
done
