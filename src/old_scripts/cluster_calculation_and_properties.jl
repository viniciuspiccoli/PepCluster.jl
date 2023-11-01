#include("./cluster_size.jl") # loading all scripts for cluster calculation

# first step:  Retrieve the trajectory data from the molecular dynamics simulations.
# The vector will include the coordinates of all peptides at each time step

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






function trajectory_data(solute_sel::String, pdbfile::String, xtcfile::String, nmols::Int64)
    atoms  = PDBTools.readPDB(pdbfile)         # selection of all atoms from the pdb
    sel0   = PDBTools.select(atoms, solute_sel) #
    solute = Selection(sel0, nmols=nmols)       # nmols = number of entities of the selection (=1 -> one molecule)

    trajectory = Trajectory(xtcfile, solute, format="XTC" ) # traj loading

    nprot   = length(trajectory.x_solute) # coordinates of all peptide atoms
    nframes = trajectory.nframes          # number of frames
    sides   =  length(getsides(trajectory, 1))  # "number of sides"
   
    data_traj = Array{SVector{3, Float64}}(undef,nframes, nprot) # saving all data


    box_sides = zeros(nframes, sides)     # vector to save all sides

    for iframe in 1:nframes
        nextframe!(trajectory) # acess coordinates of the next frame 
        data_traj[iframe,:]  = trajectory.x_solute  # matrix with protein coordinates
        box_sides[iframe,:]  = getsides(trajectory, 1)  # box sides
    end

    return data_traj, box_sides, nframes
end


#= second step: calculate pairwise distances: for each time step in the trajectory
 calculate pairwise distances between all pairs of peptides. There are various distances
 metrics: Here we will use the minimum image convention distance in periodic systems.
 
 The cutoff distance will defines when two peptides are considered to be in the
 same cluster. Each cutoff is related to the expected range of intermolecular interactions
 or the size of the peptides themselves. Here, however, we will use two distances to define the cluster_size
 in a more accurate way.
 (the first script will use ONE distance for the sake of my mental health)

 For the clustering itself, at each time step, clusters of the peptides bases on the pairwise distances.
 A regular approach is to use clustering algorithms such as hierarchical clustering, DBSCAN (density-Based Spatial Clustering of
 applications with Noise), or k-means clustering.

=#

using Distances  # For distance calculation
using Clustering  # For k-means clustering
using StaticArrays  # For SVector
using StatsBase

function calculate_cluster_size(traj_data, num_peptides, box_sides, distance_cutoff)
    num_frames = size(traj_data, 1)
    max_clusters = nmols
    
    #println("Clustering calculation...")
    cluster_sizes = []
    number_of_clusters = []
    for frame in 1:num_frames
        #println("Frame # $frame") 
        coordinates = [Array(coord) for coord in traj_data[frame, :]]

        # Apply minimum image convention
        coordinates = minimum_image_convention(coordinates, box_sides[frame, :])

        # Calculate pairwise distances using Euclidean distance
        pairwise_distances = pairwise(Euclidean(), coordinates)

        #println("performing k-means clustering")        
        # Perform k-means clustering
        k = 1  # Initial number of clusters
        largest_cluster_size = 0
        n = 0
        while k <= max_clusters 
            assignments = kmeans(pairwise_distances, k).assignments

            # Find the largest cluster size
            cluster_counts = StatsBase.countmap(assignments)
            return cluster_counts
            largest_cluster_size = maximum(values(cluster_counts))
            println(cluster_counts)
            # Increase the number of clusters until the largest cluster size becomes smaller than the cutoff
            if largest_cluster_size <= distance_cutoff
                break
            end
            n = length(cluster_counts)
            k += 1
        end
        #println("  ")
        push!(cluster_sizes, largest_cluster_size)
        push!(number_of_clusters, n)
    end

    return cluster_sizes, number_of_clusters
end

# Helper function to apply minimum image convention
function minimum_image_convention(coordinates, box_sides)
    modified_coordinates = copy(coordinates)  # Create a new array to store modified coordinates
    for i in 1:size(coordinates, 1)
        for j in 1:size(coordinates, 2)
            modified_coordinates[i][j] -= box_sides[j] * round(coordinates[i][j] / box_sides[j])
        end
    end

    return modified_coordinates
end






## testing the code

# general parameters
solute_sel   = "protein" 
nmols        =  15
pdbfile      = "processed.pdb"
xtcfile      = "processed.xtc"

# getting peptide coordinates and box sides
traj_data, box_sides, num_frames = trajectory_data(solute_sel, pdbfile, xtcfile, nmols)

# parameters for cluster calculation
num_dimensions        = 3
distance_cutoff       = 5.0


# asdasd
num_atoms = size(traj_data, 2)
atom_per_peptide = num_atoms ./ 15


#cluster calculation
cluster_sizes = calculate_cluster_size(traj_data, nmols, box_sides, distance_cutoff)


#=
Certainly! Here's an explanation of what each vector stores:

1. `cluster_max_sizes`: This vector stores the maximum cluster size for each frame in the trajectory. It represents the largest number of peptides within a single cluster in each frame. For example, if `cluster_max_sizes[3]` is 10, it means that the largest cluster in the third frame of the trajectory contains 10 peptides.

2. `cluster_total_counts`: This vector stores the total number of clusters for each frame in the trajectory. It represents the count of distinct clusters present in each frame. For example, if `cluster_total_counts[5]` is 3, it means that there are three separate clusters in the fifth frame of the trajectory.

3. `cluster_lifetimes`: This vector stores the cluster lifetime for each frame in the trajectory. It represents the duration or time span for which each cluster exists in each frame. It is calculated as the count of peptides belonging to the most abundant cluster in each frame. For example, if `cluster_lifetimes[2]` is 20, it means that the most abundant cluster in the second frame of the trajectory persisted for 20 time steps.

t(). `cluster_dynamics`: This vector stores the cluster dynamics for each frame in the trajectory. It represents the ratio of the number of clusters to the maximum cluster size in each frame. It provides an indication of the diversity or spread of clusters compared to the largest cluster. Higher values indicate a more dynamic or diverse cluster distribution, while lower values indicate a dominance of a few large clusters.

5. `cluster_rg`: This vector stores the cluster radius of gyration for each frame in the trajectory. It represents the average spatial extent or size of the clusters in each frame. It is calculated as the mean radius of gyration across all clusters present in each frame. A larger value indicates more dispersed or extended clusters, while a smaller value indicates more compact clusters.

6. `cluster_interactions`: This vector stores the cluster inter-interactions for each frame in the trajectory. It represents the total number of interactions between different clusters in each frame. It is calculated by considering pairwise interactions between all pairs of clusters and counting the number of interactions. This metric provides insights into the degree of interplay or interactions between clusters in the system.

These vectors allow you to analyze and compare various properties of the clusters across different frames in your trajectory. You can plot them, calculate statistics, or perform further analysis based on your research goals and requirements.



















=#

