echo "Running PathWithDiffRes.jl for 0.5"
nohup julia --threads=4 PathWithDiffRes.jl > curr.log
echo "Running PathWithDiffM.jl for 5"
nohup julia --threads=4 PathWithDiffM.jl 5 > curr.log
echo "Running PathWithDiffM.jl for 15"
nohup julia --threads=4 PathWithDiffM.jl 15 > curr.log