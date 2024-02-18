module PepCluster

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
  using Distances  # For distance calculation
  using Clustering  # For k-means clustering
  using StatsBase
  using Statistics

  include("Selection.jl")
  include("FileOperations.jl")  
  include("viewmol.jl") 
  include("distance.jl")

  include("Trajectory.jl")
  include("ChemFiles.jl")
  include("NamdDCD.jl")
  include("PDBTraj.jl")
   
  include("Cluster_func.jl")


  function trajectory_data(solute_sel::String, pdbfile::String, xtcfile::String, nmols::Int64)
    atoms = PDBTools.readPDB(pdbfile)         # selection of all atoms from the pdb
    sel0 = PDBTools.select(atoms, solute_sel) #
    solute = Selection(sel0, nmols=nmols)       # nmols = number of entities of the selection (=1 -> one molecule)

    trajectory = Trajectory(xtcfile, solute, format="XTC" ) # traj loading

    nprot = length(trajectory.x_solute) # coordinates of all peptide atoms
    nframes = trajectory.nframes          # number of frames
    sides =  length(getsides(trajectory, 1))  # "number of sides"

    data_traj = Array{Vector{Float64}}(undef, nframes, nprot)

    #data_traj = SVector{Float64}(undef, nframes, nprot) # saving all data
    #data_traj = Array{SVector{3, Float64}}(undef, nframes, nprot) 

    box_sides = zeros(nframes, sides)     # vector to save all sides

    for iframe in 1:nframes
        nextframe!(trajectory) # access coordinates of the next frame
        data_traj[iframe, :] = trajectory.x_solute  # matrix with protein coordinates
        box_sides[iframe, :] = getsides(trajectory, 1)  # box sides
    end

    return data_traj, box_sides, nframes
end
export trajectory_data




end
#
#  # functions
#  
#  # functions to add proper keys to a Dictionary - this is for the sake
## of the code organization
#function add_keys!(n, data::Dict{String, Vector{Int64}})
#  for key in 1:n
#    data["$key"] = []
#  end  
#end
#
#export add_keys!
#
## concatenating two vectors or numbers
#function concat(x,y)
#  return [x;y]
#end         
#export concat
#
#
## findind clusters from the matrix nodes and saving into the dict data
#function find_clusters!(nodes, data::Dict{String, Vector{Int64}})
#  for col in 1:length(nodes[:,1]) - 1 # olhando para as colunas primeiro
#    n = 0 
#    for row in col + 1:length(nodes[1,:])
#      if nodes[row,col] == 1
#        n = n + 1
#        if n == 1 
#          push!(data["$col"], col)
#        end 
#        push!(data["$col"],row)
#      end
#    end
#  end
#end 
#
#export find_clusters!
#
## function that check if a value is inside of a 1D array
#check(x, vector) = x in vector
#
#export check
#
## check function for comparison between two vectors
#function check(vector1::Vector, vector2::Vector)
#  for i in 1:length(vector1)
#    if vector1[i] in vector2
#        return true
#    end
#  end
#  return false
#end
#
## function to look through the Dictionary and find the clusters
#export get_clusters
#function get_clusters(data::Dict{String, Vector{Int64}})
#  empty    = []
#  clusters = []
#  for i in 1:Data.count             
#    if isempty(Data["$i"])
#      push!(empty, i)
#    else  
#      push!(clusters, Data["$i"])
#    end 
#  end
#  
## REMOVING PEPTIDES FROM THE EMPTY VECTOR THAT ARE IN SOME CLUSTERS AND WAS NOT ASSINGED DUE TO THE FACT THAT IS LABELED BY A HUGE INDEX
# # println("all peptides = $empty")
#  cempty = copy(empty)
#
#  for i in 1:length(cempty)
#    check_val = cempty[i]
#   # println("checking for $(cempty[i])") 
#    for j in 1:length(clusters) 
#      if check(check_val, clusters[j])
#        deleteat!(empty, findall(x->x==check_val, empty))
#       # println("must remove $check_val from empty vector")        
#      end
#    end
#   # println("   ")
#  end
#
#  # adjusting clusters to join cluster that contains the same peptide
#  return empty, clusters
#end
#
#export append_vectors
#function append_vectors(vector1, vector2)
#  new_vec=Vector{Int64}()
#  for i in 1:length(vector2)  
#    if  check(vector2[i], vector1)              
#      new_vec = union(concat(vector2, vector1))
#    end
#  end
#  return new_vec
#end
#
#teste = [[1, 11], [2, 10], [9, 12, 3, 6], [4, 13, 14], [9, 12, 6], [9, 12, 7]]
## function requeried to organize the matrix of vectors that contains the clusters
#
#export adjust
#
#function adjust(clusters)
#  # loop for unite the clusters
#  data_str = []
#  for i in 1:length(clusters)
#      new = []
#      n = 0
#      for j in i+1:length(clusters)
#        # vector to check if there is a common element between two vectors
#        if check(clusters[i], clusters[j])
#            new = append_vectors(clusters[i], clusters[j])
#            n += 1          
#        elseif check(clusters[i], clusters[j]) == false
#            continue
#        end
#
#        if n == 1 || n > 1
#          if check(new, clusters[j])
#            new = append_vectors(new, clusters[j])
#            n = n + 1
#          end
#        end
#      end
#
#      if n == 0
#        value = false
#        for k in 1:length(data_str)
#          if check(clusters[i], data_str[k])
#            value = true
#          end  
#        end    
#        if value ==false  
#          push!(data_str, clusters[i])
#        end
#      else
#        push!(data_str, new)
#      end
#  end
#  return data_str
#end
#
### calculation of distances
#
#export find_min
## calculating the minimum-distance between two peptides
#function find_min(x1::Vector{SVector{3, Float64}}, x2::Vector{SVector{3, Float64}}, sides::SVector{3, Float64}, natoms)
#  mindist = Vector{Float64}(undef, natoms) # vector to save all minimum-distances beetween the peptide1 atoms and the others peptides
#  for k in 1:natoms     # loop through arbitrarily chosen peptide 1
#    dist = +Inf
#    for i in 1:natoms   # loop through the second peptide
#      d = distance(x1[k], x2[i], sides)
#      if dist > d
#        dist = d
#      end
#    end
#    mindist[k] = dist 
#  end
#  return minimum(mindist)
#end
#
#export determine
## function to check if the distance between peptides i and j is smaller than the cutoff
#function determine(matrix1::Matrix{Float64}, cutoff)
#  new = Matrix{Any}(undef,length(matrix1[:,1]),length(matrix1[1,:]))
#  for i in 1:length(matrix1[1,:])
#    for j in 1:length(matrix1[:,1])
#      if i == j || j > i
#          new[i,j] = Nothing
#      else
#          if matrix1[i,j] < cutoff
#              new[i,j] = true
#          else
#              new[i,j] = false
#          end
#      end
#    end
#  end
#  return new
#end
#
#
#export find_neigh!
## storage = vector to save the number of neighbors for each ptptide
#function find_neigh!(storage, data)
#  for i in 1:length(data[1,:])
#    for j in 1:length(data[:,1])
#      if data[i,j] == true    
#        storage[i] += 1
#      #  println("$i is in contact with $j")  
#      end
#    end
#  end
#end 
#
#export get_indexes
## function to get the indexes of a peptide    
#get_indexes(indexes::Vector{Int64}, ipep::Int64, natoms) = indexes[ (ipep-1) * natoms + 1 : ipep*natoms ]
#get_indexes(indexes::Vector{Int64}, natoms;ipep=1) = indexes[1:natoms] 
#
#export create_dict
#function create_dict()
#  return Data = Dict{String, Vector{Int64}}() # creating dict to save data 
#end
#
#
#
#export cluster_calculation
#
## 0.6 mol / L - BGL = 227 ---- AGL == 127 
#function cluster_calculation(peptideSel::String, pdbfile::String, nmols::Int64, cutoff::Float64) 
#
#  SELsolute     = "protein"
#  atoms      = PDBTools.readPDB("./10/processed.pdb") # selection - all atoms of the PDB file.
#  
#  sel0    = PDBTools.select(atoms, SELsolute)        # Selection of which atoms from the PDB are the protein atoms.
#  protein = Selection(sel0, nmols=15)             # Selection nmols = 1 (one protein)
#  
#  trajectory = Trajectory("./10/processed.xtc", protein, format="XTC" )
#  
#  # number of protein atoms
#  nprot = length(trajectory.x_solute)
#  
#  # number of frames
#  nframes = trajectory.nframes
#                             
#  dict_new() = Dict{String, Vector{Int64}}() 
# 
#  # vector to storage the data
#  nclusters = zeros(Float64, nframes) 
#  maxclusters = zeros(Float64, nframes)
#
#  for iframe in 1:nframes #:nframes
#
#    global Data = dict_new()
#    nextframe!(trajectory) # acess coordinates of the next frame 
#    solute    = trajectory.solute    # selection  
#    x_solute  = trajectory.x_solute  # matrix with protein coordinates
#    sides = getsides(trajectory, 1)  # box sides
#  
#    # tests to define the estrategy to get coordinates for each peptide
#    natoms= solute.natomspermol
#    peptides = Vector{Vector{Int64}}(undef,solute.nmols) # vector to save indexes for all peptides  
#    
#    # saving each peptide index
#    peptides[1] =  get_indexes(solute.index, natoms) 
#    for i in 2:solute.nmols
#      peptides[i] = get_indexes(solute.index, i, natoms)
#    end
#  
#    # getting coordinates of the peptides
#    coord_matrix = []   #  Vector{SVector{3, Float64}}(undef,solute.nmols)    
#    for i in 1:solute.nmols
#      push!(coord_matrix, x_solute[peptides[i]])  
#    end
#  
#    min_dist = zeros(solute.nmols, solute.nmols)
#    
#    for i in 1:solute.nmols
#      for j in 1:solute.nmols 
#        if i == j
#          min_dist[i,j] = 0
#        elseif j > i
#          min_dist[i,j] = 0
#        else
#          min_dist[i,j] = find_min(coord_matrix[i], coord_matrix[j], sides, natoms)  
#        end   
#      end
#    end
#  
#    # info[i,j] = true, means that peptide i and j are part of the same cluster
#    cutoff = 4
#    info = determine(min_dist, cutoff)
#    number_neigh = zeros(Int64, natoms)
#    find_neigh!(number_neigh, info) # finding peptides that are close  to the ith peptide
#    add_keys!(15,Data) # adding 1 to 15 as the keys of the dict
#    find_clusters!(info, Data) # finding clusters and assigining the contacts for each peptide
#    empty, clusters = get_clusters(Data) # getting clusters of 1 elemtent (empty) and matrix disordered
#     
#    # to be fixed : sometimes it is required to perform the calculation more than one time 
#    new = adjust(clusters) # organizing the matrix
#    for i in 1:10
#      new = adjust(new)
#    end
#
#    nclusters[iframe] = length(new) + length(empty)
#  
#  end   # end of the main loop though the frames
#  
#  return nclusters 
#end
#
#
#nclusters = cluster_calculation()
#
#
#
#
#
#
#
#
#
#
#
#end
#
