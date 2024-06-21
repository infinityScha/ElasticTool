#!/bin/bash

# Loop through the arguments of js1 and js2 and bash run_vary_chi.sh
echo "Running ForHeatmap.jl with chi=1.2"
nohup julia --threads 12 ForHeatmap.jl 0.6 0.6 1.2 1.2 &>curr.log
echo "Running ForHeatmap.jl with chi=1.4"
nohup julia --threads 12 ForHeatmap.jl 0.6 0.6 1.4 1.4 &>curr.log
echo "Running ForHeatmap.jl with chi=1.6"
nohup julia --threads 12 ForHeatmap.jl 0.6 0.6 1.6 1.6 &>curr.log
#echo "Running ForHeatmap_reversed.jl with chi=1.2"
#nohup julia --threads 12 ForHeatmap_reversed.jl 0.6 0.6 1.2 1.2 &>curr.log
#echo "Running ForHeatmap_reversed.jl with chi=1.4"
#nohup julia --threads 12 ForHeatmap_reversed.jl 0.6 0.6 1.4 1.4 &>curr.log
#echo "Running ForHeatmap_reversed.jl with chi=1.6"
#nohup julia --threads 12 ForHeatmap_reversed.jl 0.6 0.6 1.6 1.6 &>curr.log
#echo "Running ForHeatmap_domain.jl with chi=1.2"
#nohup julia --threads 12 ForHeatmap_domain.jl 0.6 0.6 1.2 1.2 &>curr.log
#echo "Running ForHeatmap_domain.jl with chi=1.4"
#nohup julia --threads 12 ForHeatmap_domain.jl 0.6 0.6 1.4 1.4 &>curr.log
#echo "Running ForHeatmap_domain.jl with chi=1.6"
#nohup julia --threads 12 ForHeatmap_domain.jl 0.6 0.6 1.6 1.6 &>curr.log