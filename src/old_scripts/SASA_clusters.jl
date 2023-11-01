
#=
using PDBTools
using Statistics

function write_cluster_pdb(pdbfile, coords)
    atom_lines = []
    for (i, coord) in enumerate(coords)
        atom_line = @sprintf("ATOM  %5d  CA  UNK A%4d    %8.3f%8.3f%8.3f  1.00  0.00\n", i, i, coord[1], coord[2], coord[3])
        push!(atom_lines, atom_line)
    end

    open(pdbfile, "w") do f
        for line in atom_lines
            write(f, line)
        end
    end
end

function calculate_mean_cluster_surface_area(traj_data, cluster_assignments, pdbfile)
    num_frames = size(traj_data, 1)
    num_clusters = maximum(cluster_assignments)

    cluster_surface_areas = []

    for frame in 1:num_frames
        coordinates = traj_data[frame, :]
        assignments = cluster_assignments[frame, :]

        # Calculate surface area for each cluster
        cluster_surface_area = zeros(Float64, num_clusters)

        for cluster in 1:num_clusters
            cluster_indices = findall(assignments .== cluster)
            cluster_coords = coordinates[cluster_indices, :]

            # Create a temporary PDB file with cluster coordinates
            write_cluster_pdb(pdbfile, cluster_coords)

            # Calculate surface area using PDBTools
            atoms = readPDB(pdbfile)
            selected_atoms = selectAtoms(atoms, atom -> atom.name == "CA")
            if isempty(selected_atoms)
                error("Could not find any CA atom in the PDB file.")
            end
            surface_area = getSurfaceArea(selected_atoms)

            cluster_surface_area[cluster] = surface_area
        end

        # Calculate mean surface area for the frame
        mean_surface_area = mean(cluster_surface_area)
        push!(cluster_surface_areas, mean_surface_area)
    end

    return cluster_surface_areas
end

# Usage example
#traj_data = ...  # Replace with your trajectory data
#cluster_assignments = ...  # Replace with your cluster assignments
#pdbfile = "temp.pdb"  # Provide a temporary PDB file name

#mean_cluster_sasa = calculate_mean_cluster_surface_area(traj_data, cluster_assignments, pdbfile)

=#




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
  using StaticArrays  # For SVector
  using StatsBase
  using Statistics
  using Printf
  using LinearAlgebra
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
    atoms = PDBTools.readPDB(pdbfile)         # selection of all atoms from the pdb
    sel0 = PDBTools.select(atoms, solute_sel) #
    solute = Selection(sel0, nmols=nmols)       # nmols = number of entities of the selection (=1 -> one molecule)

    trajectory = Trajectory(xtcfile, solute, format="XTC" ) # traj loading

    nprot = length(trajectory.x_solute) # coordinates of all peptide atoms
    nframes = trajectory.nframes          # number of frames
    sides =  length(getsides(trajectory, 1))  # "number of sides"

    data_traj = Array{SVector{3, Float64}}(undef, nframes, nprot) # saving all data

    box_sides = zeros(nframes, sides)     # vector to save all sides

    for iframe in 1:nframes
        nextframe!(trajectory) # access coordinates of the next frame
        data_traj[iframe, :] = trajectory.x_solute  # matrix with protein coordinates
        box_sides[iframe, :] = getsides(trajectory, 1)  # box sides
    end

    return data_traj, box_sides, nframes
end


# Example usage

#cluster_sizes = calculate_cluster_size(traj_data, nmols, box_sides, distance_cutoff)

##cluster_max_sizes, cluster_total_counts, cluster_lifetimes, cluster_assignments = calculate_cluster_properties2(traj_data, nmols, box_sides, distance_cutoff)




using Distances
using Clustering
using LinearAlgebra

function calculate_clusters(data_traj, box_sides, nframes, distance_threshold)
    biggest_clusters = []
    total_clusters = []

    for iframe in 1:nframes
        frame_coords = data_traj[iframe, :]
        box_sides_frame = box_sides[iframe, :]


       # Apply periodic boundary conditions
       for i in 1:size(frame_coords, 1)
           for j in 1:size(frame_coords, 2)
               frame_coords[i, j] -= box_sides_frame[j] * round.(frame_coords[i, j] ./ box_sides_frame[j], RoundNearest)
           end
       end




        # Calculate pairwise distances between all peptides
        dist_matrix = pairwise(Euclidean(), frame_coords)

        # Initialize clustering variables
        clusters = []
        assigned = zeros(size(dist_matrix, 1))
        cluster_id = 1

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

function calculate_clusters(data_traj, box_sides, nframes, distance_threshold)
    biggest_clusters = []
    total_clusters = []

    for iframe in 1:nframes
        frame_coords = data_traj[iframe, :]
        box_sides_frame = box_sides[iframe, :]

        # Apply periodic boundary conditions
        for i in 1:size(frame_coords, 1)
            for j in 1:size(frame_coords, 2)
                frame_coords[i, j] -= box_sides_frame[j] * round.(frame_coords[i, j] ./ box_sides_frame[j], RoundNearest)
            end
        end


        # Calculate pairwise distances between all peptides
        dist_matrix = pairwise(Euclidean(), frame_coords)

        # Initialize clustering variables
        clusters = Vector{Vector{Int}}()
        assigned = zeros(size(dist_matrix, 1))
        cluster_id = 1

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
        #biggest_cluster = argmax(length, clusters)
        #push!(biggest_clusters, clusters[biggest_cluster])

        # Count the total number of clusters
        total_cluster_count = length(clusters)
        push!(total_clusters, total_cluster_count)
    end

    return biggest_clusters, total_clusters
end











# data
solute_sel = "protein"
nmols = 15
pdbfile = "processed.pdb"
xtcfile = "processed.xtc"

# trajectory data
traj_data, box_sides, num_frames = trajectory_data(solute_sel, pdbfile, xtcfile, nmols)

#
distance_threshold1 = 2.0
distance_threshold2 = 2.0


biggest_clusters, total_clusters = calculate_clusters(traj_data, box_sides, num_frames, distance_threshold1)

## Print the results
#for iframe in 1:nframes
#    println("Frame $iframe:")
#    println("Biggest Cluster: ", biggest_clusters[iframe])
#    println("Total Clusters: ", total_clusters[iframe])
#    println()
#end
#









