#!/bin/bash

# Define the range of arguments
js1=(0 0.1 0.2 0.3 0.4 0.5 0.6 -0.1 -0.2)

# for js = 0 just copy the results_fixed/0/0 to results_fixed/0/x for x=0.2, 0.4, 0.6, ... , 4.0
# chis=(0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0)
# for chi in ${chis[@]}
# do
#   cp -r results_fixed/0/0 results_fixed/0/${chi}
# done

# Loop through the arguments of js1 and js2 and bash run_vary_chi.sh
for i in {0..8}
do
    echo "Running run_vary_chi.sh for js=${js1[$i]}"
    nohup bash run_vary_chi.sh ${js1[$i]} ${js1[$i]} &>vary_curr_chi.log
done

