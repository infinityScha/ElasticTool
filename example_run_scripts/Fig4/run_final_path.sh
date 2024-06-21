echo "Running MakePath.jl for mixed"
nohup julia --threads=8 MakePath2.jl mixed > curr.log
echo "Running MakePath.jl for upper"
nohup julia --threads=8 MakePath2.jl upper > curr.log
echo "Running MakePath.jl for both"
nohup julia --threads=8 MakePath2.jl both > curr.log
echo "Running MakePath.jl for 50"
nohup julia --threads=8 MakePath2.jl 50 > curr.log
echo "Running MakePath.jl for none"
nohup julia --threads=8 MakePath2.jl none > curr.log