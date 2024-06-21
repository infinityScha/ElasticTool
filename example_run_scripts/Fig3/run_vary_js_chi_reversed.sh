#!/bin/bash

# Define the range of arguments
js1=(0.6 0.5 0.4 0.3 0.2 0.1 0 -0.1 -0.2)

# Loop through the arguments of js1 and js2 and bash run_vary_chi.sh
for i in {0..8}
do
    echo "Running run_vary_chi_reversed.sh for js=${js1[$i]}"
    nohup bash run_vary_chi_reversed.sh ${js1[$i]} ${js1[$i]} &>vary_curr_chi.log
done

