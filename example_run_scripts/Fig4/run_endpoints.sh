echo "Running MinimizeEndPoints.jl for mixed and initial"
nohup julia --threads=4 MinimizeEndPoints.jl mixed initial > curr.log
echo "Running MinimizeEndPoints.jl for mixed and final"
nohup julia --threads=4 MinimizeEndPoints.jl mixed final > curr.log