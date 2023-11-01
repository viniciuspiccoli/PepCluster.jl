# SPDX-License-Identifier: MIT

using PepCluster, Statistics, Clustering, Distances, Plots


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


# Iterate over frames in the trajectory
for frame in 1:nframes
    #frame = 1
    println("    ") 
    println("Frame $frame") 
    coords = hcat(peptide_coords[frame, :]...)  # Ensure that coords is a 3 x N matrix
   
    # Calculate the centers of mass for each peptide
    coms = [mean(coords[:, j:j+natoms_per_peptide-1]; dims=2) for j in 1:natoms_per_peptide:natoms_total]
   
    # Convert list of vectors to a matrix
    coms_matrix = hcat(coms...)
   
    # If coms_matrix is not 3 x N, transpose it
    if size(coms_matrix, 1) >= size(coms_matrix, 2)
       coms_matrix = coms_matrix'
    end
    # Calculate distance matrix
    distance_matrix = pairwise(Euclidean(),coms_matrix, dims=2)
    # for k in k_range
    optimal_k = 5
    final_result = kmeans(distance_matrix, optimal_k)
    # Sizes of clusters
    cluster_sizes = counts(final_result)
    # Smallest and Largest Clusters
    smallest_cluster_size = minimum(cluster_sizes)
    largest_cluster_size = maximum(cluster_sizes)
    #Update cluster size frequencies matrix
    for (i, size) in enumerate(cluster_sizes)
       if size <= max_size
           cluster_size_freq[i, size] += 1
       else
           @warn "Cluster size is larger than maximum size, it's ignored."
       end
    end
   
    println("Smallest cluster size: ", smallest_cluster_size)
    println("Largest cluster size: ", largest_cluster_size)
   
    # To get the indices of peptides in each cluster
    for i in 1:optimal_k
      println("Peptides in Cluster $i: ", findall(x -> x == i, assignments(final_result)))
    end
    
    # Plot to find the elbow point
    #plot(k_range, wcss, xlabel="Number of Clusters", ylabel="WCSS", title="Elbow Method", legend=false)
        # Analyze clusters for this frame
    #    println("Frame: ", frame)
    #    println("Number of Clusters: ", length(clusters))
        # Additional analysis can be added here
end
    
    
#end
# Normalize the matrix to convert counts to probabilities
cluster_size_freq_normalized = cluster_size_freq ./ sum(cluster_size_freq)

# Find the minimum and maximum frequencies for the color bar limits
min_freq, max_freq = extrema(cluster_size_freq_normalized)
#cluster_size_freq_normalized = cluster_size_freq / sum(cluster_size_freq)

# Plot the heatmap
heatmap(cluster_size_freq_normalized,
        c=:plasma,
        xlabel="Cluster Size",
        ylabel="Cluster Index",
        title="Cluster Size Distribution",
        color=:auto,
        xlims=(0,8),
        ylims=(0,6),  
        clims=(min_freq, max_freq),
        colorbar_title="Frequency")

savefig("cluster.png")
    

