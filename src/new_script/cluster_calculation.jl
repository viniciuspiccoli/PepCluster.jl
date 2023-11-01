import Pkg; Pkg.activate("ComputeDistances", shared=true)

# loading packages
using LinearAlgebra: norm
using StaticArrays
using PDBTools
using Plots
import CellListMap: wrap_relative_to
import Chemfiles

# function that computes the distance between atoms for each frame
function compute_distance(;pdbname, trajectory_name, selection1, selection2)
    pdb = readPDB(pdbname)
    dcd = Chemfiles.Trajectory(trajectory_name)
    atom1 = select(pdb, selection1)
    atom2 = select(pdb, selection2)
    distances = Float64[]
    for frame in dcd
        # read unit cell matrix from trajectory
        matrix_read = Chemfiles.matrix(Chemfiles.UnitCell(frame))
        unit_cell = SMatrix{3,3}(transpose(matrix_read))
        # read coordinates, convert them to small static vectors
        coordinates = Chemfiles.positions(frame)
        x = SVector{3}(@view(coordinates[:,atom1[1].index]))
        y = SVector{3}(@view(coordinates[:,atom2[1].index]))
        # wrap coordinates according to PBC. Note that here we use an
        # internal function of `CellListMap`. 
        y_wrapped = wrap_relative_to(y,x,unit_cell)
        # compute distance
        d = norm(y_wrapped-x)
        # add data to distance array
        push!(distances, d)
    end
    return distances
end

# Gromacs trajectory file
function run_gromacs()
    distances = compute_distance(
        pdbname = "./processed.pdb", 
        trajectory_name ="./processed.xtc",
        selection1 = "index 3",
        selection2 = "index 9200",)
    Plots.default(fontfamily="Computer Modern")
    plt = plot(distances, xlabel = "frame", ylabel = "distance / Å",
        framestyle=:box, linewidth=2, label=nothing)
    
         savefig("./plot.svg") # to save the figure
    display(plt)
    return plt, distances
end

## Criteria for cluster analysis
#=
distance between index 3   and index 92  between: 3 and 7 Angs
                 index 118 and index 207 between: 3 and 7 Angs
=#


