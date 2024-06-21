#!/bin/bash

# Define the range of arguments
v1=(0.11 0.12 0.09 0.08)
v1_prev=(0.1 0.11 0.1 0.09)

# Loop through the arguments of js1 and js1_prev and run ForHeatmap.jl
for i in {0..3}
do
    echo "Running MakePath.jl with v1=${v1[$i]} and v1_prev=${v1_prev[$i]}"
    nohup julia --threads 2 MakePath.jl ${v1_prev[$i]} ${v1[$i]} 0 0 &>curr.log
done
