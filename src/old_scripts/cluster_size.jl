
  using PDBTools 
  using FortranFiles
  using Printf
  using StructTypes
  using StaticArrays
  using CellListMap
  using LinearAlgebra
  using LsqFit
  using Polynomials
  using EasyFit 

#  export distance 

  include("Selection.jl")
  include("FileOperations.jl")  
  include("viewmol.jl") 
  include("distance.jl")

  include("Trajectory.jl")
  include("ChemFiles.jl")
  include("NamdDCD.jl")
  include("PDBTraj.jl")
#end
