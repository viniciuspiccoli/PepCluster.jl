
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

function calculate_cluster_size(traj_data, num_peptides, box_sides, distance_cutoff)
    num_frames = size(traj_data, 1)
    max_clusters = num_peptides
    
    cluster_sizes = Vector{Int64}(undef, num_frames)
    for frame in 1:num_frames
        coordinates = [Array(coord) for coord in traj_data[frame, :]]

        # Apply minimum image convention
        coordinates = minimum_image_convention(coordinates, box_sides[frame, :])

        # Calculate pairwise distances using Euclidean distance
        pairwise_distances = pairwise(Euclidean(), coordinates)

        # Perform k-means clustering
        k = 1  # Initial number of clusters
        largest_cluster_size = 0
        while k <= max_clusters
            assignments = kmeans(pairwise_distances, k).assignments

            # Find the largest cluster size
            cluster_counts = countmap(assignments)
            largest_cluster_size = maximum(values(cluster_counts))

            # Increase the number of clusters until the largest cluster size becomes smaller than the cutoff
            if largest_cluster_size <= distance_cutoff
                break
            end

            k += 1
        end

        cluster_sizes[frame] = largest_cluster_size
    end

    return cluster_sizes
end

function minimum_image_convention(coordinates, box_sides)
    modified_coordinates = copy(coordinates)  # Create a new array to store modified coordinates
    for i in 1:size(coordinates, 1)
        for j in 1:size(coordinates, 2)
            modified_coordinates[i][j] -= box_sides[j] * round(coordinates[i][j] / box_sides[j])
        end
    end

    return modified_coordinates
end




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




function calculate_cluster_properties(traj_data, num_peptides, box_sides, distance_cutoff)
    num_frames = size(traj_data, 1)
    cluster_max_sizes = Vector{Int64}(undef, num_frames)
    cluster_total_counts = Vector{Int64}(undef, num_frames)
    cluster_lifetimes = Vector{Int64}(undef, num_frames)
     
    # vectors for SASA calculation
    cluster_assignments  = []  # Array to store SASA values for each cluster

    for frame in 1:num_frames
        println("For frame number $frame")
        coordinates = [Array(coord) for coord in traj_data[frame, :]]
        println("passed coord")
        # Apply minimum image convention
        coordinates = minimum_image_convention(coordinates, box_sides[frame, :])
        println(" minimum image convention")
        # Calculate pairwise distances using Euclidean distance
        pairwise_distances = pairwise(Euclidean(), coordinates)
        println("pairwise distances")
        # Perform k-means clustering
        k = 1  # Initial number of clusters
        largest_cluster_size = 0
        cluster_counts = countmap(fill(1, size(coordinates, 1)))  # Initialize cluster counts
        assignments = []
        while k <= num_peptides
            assignments = kmeans(pairwise_distances, k).assignments
            println("begin while loop")
            # Find the largest cluster size
            cluster_counts = countmap(assignments)
            largest_cluster_size = maximum(values(cluster_counts))

            # Increase the number of clusters until the largest cluster size becomes smaller than the cutoff
            if largest_cluster_size <= distance_cutoff
                break
            end

            k += 1

        end

        cluster_max_sizes[frame] = largest_cluster_size
        println("cluster sizes for frame $frame")
        #println(unique(collect(assignments)))
        cluster_total_counts[frame] = length(unique(collect(assignments)))
        println("total conts for frame $frame")
        cluster_lifetimes[frame] = maximum(values(cluster_counts))
        println("lifetime for frame $frame")
        push!(cluster_assignments, assignments)
        


    end

    return cluster_max_sizes, cluster_total_counts, cluster_lifetimes, cluster_assignments
end

#new function

function calculate_cluster_properties2(traj_data, num_peptides, box_sides, distance_cutoff)
    num_frames = size(traj_data, 1)
    cluster_max_sizes = Vector{Int64}(undef, num_frames)
    cluster_total_counts = Vector{Int64}(undef, num_frames)
    cluster_assignments = []  # Array to store SASA values for each cluster

    for frame in 1:num_frames
        coordinates = [Array(coord) for coord in traj_data[frame, :]]
        # Apply minimum image convention
        coordinates = minimum_image_convention(coordinates, box_sides[frame, :])
        # Calculate pairwise distances using Euclidean distance
        pairwise_distances = pairwise(Euclidean(), coordinates)

        println(pairwise_distances)
        # Perform k-means clustering
        k = 1  # Initial number of clusters
        largest_cluster_size = 0
        cluster_counts = countmap(fill(1, size(coordinates, 1)))  # Initialize cluster counts
        assignments = []
        while k <= num_peptides
            assignments = kmeans(pairwise_distances, k).assignments
            # Find the largest cluster size
            cluster_counts = countmap(assignments)
            largest_cluster_size = maximum(values(cluster_counts))

            # Increase the number of clusters until the largest cluster size becomes smaller than the cutoff
            if largest_cluster_size <= distance_cutoff
                break
            end

            k += 1
        end


        cluster_max_sizes[frame] = largest_cluster_size
        # Calculate cluster lifetime and total counts
        cluster_counts_sorted = sort(collect(cluster_counts), rev=true)

        # Filter cluster counts to consider only clusters with more than one peptide
        #filtered_counts = filter(x -> x > 1, values(cluster_counts))
        #cluster_total_counts[frame] = length(filtered_counts)
 
        println(cluster_counts) 
        filtered_counts = [count for count in values(cluster_counts) if count > 1]
        cluster_total_counts[frame] = length(filtered_counts)

        #cluster_total_counts[frame] = length(cluster_counts_sorted)
        push!(cluster_assignments, assignments)
    end

    return cluster_max_sizes, cluster_total_counts, cluster_assignments
end
















# Usage example
solute_sel = "protein"
nmols = 15
pdbfile = "processed.pdb"
xtcfile = "processed.xtc"

traj_data, box_sides, num_frames = trajectory_data(solute_sel, pdbfile, xtcfile, nmols)

distance_cutoff = 2.0

#cluster_sizes = calculate_cluster_size(traj_data, nmols, box_sides, distance_cutoff)

cluster_max_sizes, cluster_total_counts, cluster_lifetimes, cluster_assignments = calculate_cluster_properties2(traj_data, nmols, box_sides, distance_cutoff)









