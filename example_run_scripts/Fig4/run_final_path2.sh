#echo "Running MakePath.jl for 20"
#nohup julia --threads=4 MakePath3.jl 20 > curr.log
#echo "Running MakePath.jl for 30"
#nohup julia --threads=4 MakePath3.jl 30 > curr.log
#echo "Running MakePath.jl for 40"
#nohup julia --threads=4 MakePath3.jl 40 > curr.log
echo "Running MakePath.jl for 15"
nohup julia --threads=4 MakePath3.jl 15 > curr.log
echo "Running MakePath.jl for 5"
nohup julia --threads=4 MakePath3.jl 5 > curr.log