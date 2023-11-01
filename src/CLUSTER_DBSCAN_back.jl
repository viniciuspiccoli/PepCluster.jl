# SPDX-License-Identifier: MIT

using PepCluster, Statistics, Clustering, Distances, Plots

using PDBTools

#function save_cluster_to_pdb(atoms, indices::Vector{Int}, output_pdb::String)
#    cluster_atoms = []
#    for idx in indices
#        push!(cluster_atoms, atoms[idx])
#    end
#    writePDB(cluster_atoms, output_pdb)
#end

function save_cluster_to_pdb(atoms::Vector{Atom}, indices::Vector{Int64}, output_pdb::String)
    cluster_atoms = [atoms[i] for i in indices]  # Ensure this is a Vector{Atom}
    writePDB(cluster_atoms, output_pdb)
end


### some function to extract data from the dbscan clsuter struct final_result
function cluster_parameters(data::DbscanResult)
    nclusters = length(data.counts)
    max_size  = maximum(data.counts)
    size = data.counts ./ 115    
    return nclusters, max_size, size
end
###


# Assuming a maximum of 15 clusters and maximum size of 15
max_clusters = 15
max_size = 15
cluster_size_freq = zeros(Int, max_clusters, max_size)

# data
solute_sel = "protein"
nmols = 15
pdbfile = "processed.pdb"
xtcfile = "processed.xtc"

# trajectory data
peptide_coords, box_sides, nframes = trajectory_data(solute_sel, pdbfile, xtcfile, nmols)
distance_threshold = 3.5  # In Ångströms
natoms_total = size(peptide_coords, 2)    # 3 for x, y, z coordinates

if natoms_total % nmols != 0
    error("Total number of atoms is not divisible by the number of peptides. Please check your input data.")
end
natoms_per_peptide = natoms_total ÷ nmols

atoms = PDBTools.readPDB(pdbfile)
# Inside your loop over frames
#for frame in 1:nframes
 frame = 1
    println("\nFrame $frame")
    coords = hcat(peptide_coords[frame, :]...)

    # Apply DBSCAN
    eps = 3.5  # Example value, you might need to adjust this
    minpts = 2  # Example value, you might need to adjust this
    #db_result = dbscan(distance_matrix, eps, minpts)
    db_result = dbscan(coords, eps, min_neighbors=1, min_cluster_size=2)
    
    # Identify the biggest cluster
    biggest_cluster_idx = argmax(db_result.counts)
    biggest_cluster_points = findall(x -> db_result.assignments[x] == biggest_cluster_idx, 1:natoms_total)
  
    # Extract the atom coordinates corresponding to this cluster
    cluster_coords = coords[:, biggest_cluster_points]
  
    # Save to PDB
    if frame == 173
        output_pdb = "biggest_cluster_frame_$(frame).pdb"
        save_cluster_to_pdb(atoms, biggest_cluster_points, output_pdb)
    end
    nc, max, sizes = cluster_parameters(db_result)
    println("Total number of cluster = $nc")
    println("maximum size = $max") 
#end

