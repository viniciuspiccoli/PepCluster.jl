



  using Distances
  using Clustering
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




function calculate_clusters(data_traj, box_sides, nframes, distance_threshold)
    biggest_clusters = []
    total_clusters   = []

    for iframe in 1:nframes
        frame_coords = data_traj[iframe, :]
        box_sides_frame = box_sides[iframe, :]

        # Apply periodic boundary conditions
        for i in 1:size(frame_coords, 1)
            for j in 1
                frame_coords[i, j] -= box_sides_frame[j] * round.(frame_coords[i, j] ./ box_sides_frame[j], RoundNearest)
            end
        end

        println("periodic boundary conditions applied")    
        # Calculate pairwise distances between all peptides
        dist_matrix = pairwise(Euclidean(), frame_coords)

        # Initialize clustering variables
        clusters = Vector{Vector{Int}}()
        assigned = zeros(size(dist_matrix, 1))
        cluster_id = 1

        println("initialized vectors of the clusters")
        # Cluster peptides based on distance threshold
        for i in 1:size(dist_matrix, 1)
            if assigned[i] == 0
                cluster = [i]
                assigned[i] = cluster_id

                for j in (i+1):size(dist_matrix, 2)
                    if dist_matrix[i, j] < distance_threshold
                        push!(cluster, j)
                        assigned[j] = cluster_id
                    end
                end

                push!(clusters, cluster)
                cluster_id += 1
            end
        end

        # Find the biggest cluster in the frame
        biggest_cluster = argmax(length, clusters)
        push!(biggest_clusters, clusters[biggest_cluster])

        # Count the total number of clusters
        total_cluster_count = length(clusters)
        push!(total_clusters, total_cluster_count)
    end

    return biggest_clusters, total_clusters
end

# Example usage
data_traj, box_sides, nframes = trajectory_data("protein", "processed.pdb", "processed.xtc", 15)
distance_threshold = 5.0  # Distance threshold for clustering

biggest_clusters, total_clusters = calculate_clusters(data_traj, box_sides, nframes, distance_threshold)

# Print the results
for iframe in 1:nframes
    println("Frame $iframe:")
    println("Biggest Cluster: ", biggest_clusters[iframe])
    println("Total Clusters: ", total_clusters[iframe])
    println()
end
