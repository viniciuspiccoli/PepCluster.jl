
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
using BioStructures # for SASA and print of PDBS


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

# Usage example
solute_sel = "protein"
nmols = 15
pdbfile = "processed.pdb"
xtcfile = "processed.xtc"

traj_data, box_sides, num_frames = trajectory_data(solute_sel, pdbfile, xtcfile, nmols)



using Distances
using Clustering

function calculate_clusters(data_traj, box_sides)
    nframes = size(data_traj, 1)
    npeptides = size(data_traj, 2) ÷ 115  # Assuming each peptide has 115 atoms
    
    cluster_sizes = Vector{Int}(undef, nframes)
    num_clusters = Vector{Int}(undef, nframes)
    
    for frame in 1:nframes
        # Create a distance matrix for the pairwise distances between peptides
        distances = pairwise(Euclidean(), reshape(data_traj[frame, :], 3, npeptides)')
        
        # Apply periodic boundary conditions to distances
        distances = apply_pbc(distances, box_sides[frame, :])
        
        # Initialize an array to store the cluster assignments
        clusters = zeros(Int, npeptides)
        cluster_count = 0
        
        for i in 1:npeptides
            if clusters[i] == 0  # Not assigned to any cluster yet
                cluster_count += 1
                clusters[i] = cluster_count
                
                # Check distances with other peptides to assign them to the same cluster
                for j in (i + 1):npeptides
                    if clusters[j] == 0 && distances[i, j] < 5.0
                        clusters[j] = cluster_count
                    end
                end
            end
        end
        
        # Calculate cluster sizes
        cluster_sizes[frame] = maximum([count(x -> x == c, clusters) for c in 1:cluster_count])
        num_clusters[frame] = cluster_count
    end
    
    return cluster_sizes, num_clusters
end

function apply_pbc(distances, box_sides)
    # Compute periodic distances using minimum image convention
    distances_pbc = copy(distances)
    
    for i in 1:size(distances, 1)
        for j in 1:size(distances, 2)
            dx = distances[i, j, 1]
            dy = distances[i, j, 2]
            dz = distances[i, j, 3]
            
            dx -= box_sides[1] * round(dx / box_sides[1])
            dy -= box_sides[2] * round(dy / box_sides[2])
            dz -= box_sides[3] * round(dz / box_sides[3])
            
            distances_pbc[i, j, 1] = dx
            distances_pbc[i, j, 2] = dy
            distances_pbc[i, j, 3] = dz
        end
    end
    
    return distances_pbc
end


function calculate_clusters(data_traj, box_sides)
    nframes = size(data_traj, 1)
    npeptides = size(data_traj, 2) ÷ 115  # Assuming each peptide has 115 atoms
    
    cluster_sizes = Vector{Int}(undef, nframes)
    num_clusters = Vector{Int}(undef, nframes)
    
    for frame in 1:nframes
        # Extract peptide coordinates for the current frame
        peptide_coords = data_traj[frame, 1:npeptides*115, :]
        
        # Calculate pairwise distances between peptides
        distances = zeros(npeptides, npeptides)
        
        for i in 1:npeptides
            for j in (i + 1):npeptides
                # Calculate distances between the first two alpha carbons and the last two alpha carbons
                dist_1 = norm(periodic_diff(peptide_coords[(i-1)*115+3, :], peptide_coords[(j-1)*115+3, :], box_sides[frame, :]))
                dist_2 = norm(periodic_diff(peptide_coords[(i-1)*115+92, :], peptide_coords[(j-1)*115+92, :], box_sides[frame, :]))
                
                if dist_1 < 5.0 && dist_2 < 5.0
                    # Assign peptides to the same cluster
                    distances[i, j] = 1
                    distances[j, i] = 1
                end
            end
        end
        
        # Apply transitive closure to determine cluster assignments
        cluster_assignments = transitive_closure(distances)
        
        # Count clusters and calculate cluster sizes
        cluster_count = maximum(cluster_assignments)
        cluster_sizes[frame] = maximum([count(x -> x == c, cluster_assignments) for c in 1:cluster_count])
        num_clusters[frame] = cluster_count
    end
    
    return cluster_sizes, num_clusters
end



function round_custom(x::Float64)
    if x < 0
        return Int(floor(x))
    else
        return Int(ceil(x))
    end
end

function periodic_diff(coord1, coord2, box_sides)
    # Calculate the periodic difference between two coordinates

    diff = coord1 .- coord2
    box_sides_broadcasted = reshape(box_sides, (1, 1, length(box_sides)))
    diff = diff .- box_sides_broadcasted .* map(round_custom, diff ./ box_sides_broadcasted)

    return diff
end














# Calculate clusters
cluster_sizes, num_clusters = calculate_clusters(traj_data, box_sides)

# Print cluster information for each frame
for frame in 1:nframes
    println("Frame $frame: Number of clusters = $(num_clusters[frame]), Maximum cluster size = $(cluster_sizes[frame])")
end















#distance_cutoff = 2.0

#cluster_sizes = calculate_cluster_size(traj_data, nmols, box_sides, distance_cutoff)

#cluster_max_sizes, cluster_total_counts, cluster_lifetimes = calculate_cluster_properties(traj_data, nmols, box_sides, distance_cutoff)











#=
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



# function of properties
function calculate_cluster_properties(traj_data, num_peptides, box_sides, distance_cutoff)
  num_frames = size(traj_data, 1)
  cluster_max_sizes = Vector{Int64}(undef, num_frames)
  cluster_total_counts = Vector{Int64}(undef, num_frames)
  cluster_lifetimes = Vector{Int64}(undef, num_frames)
  cluster_rg = Vector{Float64}(undef, num_frames)
  cluster_assignments = Vector{Vector{Int64}}(undef, num_frames)
  cluster_counts = Dict{Int64, Int64}()  # To store cluster counts for each frame

  max_clusters = num_peptides  # Maximum number of clusters to try

  for frame in 1:num_frames
      coordinates = [Array(coord) for coord in traj_data[frame, :]]

      # Apply minimum image convention
      coordinates = minimum_image_convention(coordinates, box_sides[frame, :])

      # Calculate pairwise distances using Euclidean distance
      pairwise_distances = pairwise(Euclidean(), coordinates)

      # Initialize assignments outside the while loop
      assignments = Vector{Int64}()

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

      # Calculate cluster properties
      num_clusters = length(unique(assignments))
      cluster_max_sizes[frame] = largest_cluster_size
      cluster_total_counts[frame] = num_clusters
      cluster_lifetimes[frame] = countmap(assignments)[findmax(assignments)[2]]
      cluster_rg[frame] = calculate_cluster_rg(coordinates, assignments)

      cluster_assignments[frame] = assignments
      cluster_counts[frame] = num_clusters
  end

  return cluster_max_sizes, cluster_total_counts, cluster_lifetimes, cluster_rg, cluster_assignments, cluster_counts
end


# testing


# general parameters
solute_sel   = "protein" 
nmols        =  15
pdbfile      = "processed.pdb"
xtcfile      = "processed.xtc"

# getting peptide coordinates and box sides
traj_data, box_sides, num_frames = trajectory_data(solute_sel, pdbfile, xtcfile, nmols)

# parameters for cluster calculation
num_dimensions        = 3
distance_cutoff       = 2.0


# asdasd
num_atoms = size(traj_data, 2)
atom_per_peptide = num_atoms ./ 15





cluster_max_sizes, cluster_total_counts, cluster_lifetimes, cluster_rg, cluster_assignments, cluster_counts = calculate_cluster_properties(traj_data, nmols, box_sides, distance_cutoff)




















using Distances
using Clustering
using LinearAlgebra
using PDBTools

function minimum_image_convention(coordinates, box_sides)
  num_atoms = size(coordinates, 1)
  corrected_coordinates = similar(coordinates)

  for i in 1:num_atoms
      for j in 1:3
          corrected_coordinates[i, j] = coordinates[i, j] - box_sides[j] * round(coordinates[i, j] / box_sides[j])
      end
  end

  return corrected_coordinates
end

function calculate_cluster_properties(traj_data, num_peptides, box_sides, distance_cutoff)
  num_frames = size(traj_data, 1)
  cluster_max_sizes = Vector{Int64}(undef, num_frames)
  cluster_total_counts = Vector{Int64}(undef, num_frames)
  cluster_lifetimes = Vector{Int64}(undef, num_frames)
  cluster_rg = Vector{Float64}(undef, num_frames)
  cluster_assignments = Vector{Vector{Int64}}(undef, num_frames)

  max_clusters = num_peptides  # Maximum number of clusters to try

  for frame in 1:num_frames
      coordinates = [Array(coord) for coord in traj_data[frame, :]]

      # Apply minimum image convention
      coordinates = minimum_image_convention(coordinates, box_sides[frame, :])

      # Calculate pairwise distances using Euclidean distance
      pairwise_distances = pairwise(Euclidean(), coordinates)

      # Initialize assignments outside the while loop
      assignments = Vector{Int64}()

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

      # Calculate cluster properties
      num_clusters = length(unique(assignments))
      cluster_max_sizes[frame] = largest_cluster_size
      cluster_total_counts[frame] = num_clusters
      cluster_lifetimes[frame] = countmap(assignments)[findmax(assignments)[2]]
      cluster_rg[frame] = calculate_cluster_rg(coordinates, assignments)

      cluster_assignments[frame] = assignments
  end

  return cluster_max_sizes, cluster_total_counts, cluster_lifetimes, cluster_rg, cluster_assignments
end

function calculate_cluster_rg(coordinates, assignments)
  num_clusters = length(unique(assignments))
  rg_sum = 0.0

  for cluster in 1:num_clusters
      cluster_coords = coordinates[assignments .== cluster]
      rg_sum += radius_of_gyration(cluster_coords)
  end

  return rg_sum / num_clusters
end

function radius_of_gyration(coordinates)
  center_of_mass = mean(coordinates, dims=1)
  sq_distances = sum((coordinates .- center_of_mass).^2, dims=1)
  return sqrt(mean(sq_distances))
end
=#