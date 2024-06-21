#!/bin/bash

# get js and prev js from the command line
js=($1)
prev_js=($2)

# Define the range of arguments
chi1=(2.8 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2 0)
chi1_prev=(3.0 2.8 2.6 2.4 2.2 2.0 1.8 1.6 1.4 1.2 1.0 0.8 0.6 0.4 0.2)

# Loop through the arguments of chi1 and chi1_prev and run ForHeatmap.jl
for i in {0..14}
do
    echo "Running ForHeatmap_reversed.jl with chi1=${chi1[$i]} and chi1_prev=${chi1_prev[$i]}"
    nohup julia --threads 2 ForHeatmap_reversed.jl ${prev_js} ${js} ${chi1_prev[$i]} ${chi1[$i]}
done
