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
## frame = 1
#    println("\nFrame $frame")
#    coords = hcat(peptide_coords[frame, :]...)
#
#    # Calculate the centers of mass for each peptide
#    #coms = [mean(coords[:, j:j+natoms_per_peptide-1]; dims=2) for j in 1:natoms_per_peptide:natoms_total]
#    #coms_matrix = hcat(coms...)
#
#    #if size(coms_matrix, 1) >= size(coms_matrix, 2)
#    #    coms_matrix = coms_matrix'
#    #end
#
#
#    # Ensure that coms_matrix is a 3 x N matrix
#    #if size(coms_matrix, 1) > size(coms_matrix, 2)
#    #    coms_matrix = coms_matrix'
#    #end
#
#
#    # Apply DBSCAN
#    eps = 3.5  # Example value, you might need to adjust this
#    minpts = 2  # Example value, you might need to adjust this
#    #db_result = dbscan(distance_matrix, eps, minpts)
#    db_result = dbscan(coords, eps, min_neighbors=1, min_cluster_size=2)
#    
#    # Identify the biggest cluster
#    biggest_cluster_idx = argmax(db_result.counts)
#    biggest_cluster_points = findall(x -> db_result.assignments[x] == biggest_cluster_idx, 1:natoms_total)
#  
#    # Extract the atom coordinates corresponding to this cluster
#    cluster_coords = coords[:, biggest_cluster_points]
#  
#    # Save to PDB
#    if frame == 173
#        output_pdb = "biggest_cluster_frame_$(frame).pdb"
#        save_cluster_to_pdb(atoms, biggest_cluster_points, output_pdb)
#    end
#    nc, max, sizes = cluster_parameters(db_result)
#    println("Total number of cluster = $nc")
#    println("maximum size = $max") 
#end


using DataFrames, Plots, Plots.Measures, LaTeXStrings

function analyze_clusters(peptide_coords, box_sides, nframes, nmols, distance_threshold)
    natoms_total = size(peptide_coords, 2)
    if natoms_total % nmols != 0
        error("Total number of atoms is not divisible by the number of peptides. Please check your input data.")
    end
    natoms_per_peptide = natoms_total ÷ nmols

    # Initialize variables to store results
    all_cluster_sizes = Int[]
    #cluster_size_freq = Dict{Int,Int}()
    cluster_size_freq = zeros(Int, 15, 15)
    total_clusters_per_frame = Int[]

    for frame in 1:nframes
      #  println("\nFrame $frame")
        coords = hcat(peptide_coords[frame, :]...)
        db_result = dbscan(coords, distance_threshold, min_neighbors=1, min_cluster_size=2)
        push!(total_clusters_per_frame, length(db_result.counts))

        for cluster_size in db_result.counts
            push!(all_cluster_sizes, cluster_size)
            cluster_size = cluster_size ./ 115
            for (i, size) in enumerate(cluster_size)
                if size <= 15
                    cluster_size_freq[i, Int64(size)] += 1
                else
                    @warn "Cluster size is larger than maximum size, it's ignored."
                end
            end
            #cluster_size_freq[cluster_size] = get(cluster_size_freq, cluster_size, 0) + 1
        

        end
    end

    return all_cluster_sizes, cluster_size_freq, total_clusters_per_frame
end



function plot_clusters(all_cluster_sizes, cluster_size_freq, total_clusters_per_frame)
    # Plot total number of clusters in each frame


    plot_font="Computer Modern"
    default(legend=:topright,
       fontfamily=plot_font,
       linewidth=2, framestyle=:box,
       grid=false)
    scalefontsizes(1.2)

    plot(layout=(2,2), size=(1000,1000))
    plot!(subplot=1, 1:length(total_clusters_per_frame), total_clusters_per_frame, title="Total Clusters per Frame", xlabel="Frame", ylabel="Total Clusters", legend=false)

    # Plot distribution of cluster sizes
    histogram!(subplot=2, all_cluster_sizes, bins=maximum(all_cluster_sizes)-1, title="Distribution of Cluster Sizes", xlabel="Cluster Size", ylabel="Frequency", legend=false)

    # Prepare data for heat map
      # Prepare data for bar plot
      cluster_sizes = collect(keys(cluster_size_freq))
      frequencies = [cluster_size_freq[size] for size in cluster_sizes]
  
      # Create a bar plot for cluster size vs frequency
#bar!(subplot=3, cluster_sizes, frequencies, title="Cluster Size vs Frequency", xlabel="Cluster Size", ylabel="Frequency", legend=false)
  
    
# Normalize the matrix to convert counts to probabilities
cluster_size_freq_normalized = cluster_size_freq ./ sum(cluster_size_freq)

# Find the minimum and maximum frequencies for the color bar limits
min_freq, max_freq = extrema(cluster_size_freq_normalized)
#cluster_size_freq_normalized = cluster_size_freq / sum(cluster_size_freq)

# Plot the heatmap
heatmap!(subplot=3, cluster_size_freq_normalized,
        c=:plasma,
        xlabel="Cluster Size",
        ylabel="Cluster Index",
        title="Cluster Size Distribution",
        color=:auto,
        xlims=(0,8),
        ylims=(0,6),  
        clims=(min_freq, max_freq),
        colorbar_title="Frequency")



## aqui eu preciso adicionar esta conta no dbconuts

# Assuming a maximum of 15 clusters and maximum size of 15
#max_clusters = 15
#max_size = 15
#cluster_size_freq = zeros(Int, max_clusters, max_size)

#Update cluster size frequencies matrix
#for (i, size) in enumerate(cluster_sizes)
#    if size <= max_size
#        cluster_size_freq[i, size] += 1
#    else
#        @warn "Cluster size is larger than maximum size, it's ignored."
#    end
#end


    
   # heatmap_data = zeros(Int, maximum(keys(cluster_size_freq)), 2)
   # for (size, freq) in cluster_size_freq
   #     heatmap_data[size, 1] = size
   #     heatmap_data[size, 2] = freq
  #  end
    
    #max_cluster_size = maximum(keys(cluster_size_freq))
    #heatmap_data = zeros(Int, max_cluster_size)
    #for (size, freq) in cluster_size_freq
    #    heatmap_data[size] = freq
   # end
   
    #heatmap_data_2d = reshape(heatmap_data[:, 2], 1, :)

    #heatmap!(subplot=3, heatmap_data_2d, title="Cluster Size vs Frequency", xlabel="Cluster Size", ylabel="Frequency", color=:viridis, legend=false)
   
    #heatmap!(subplot=3, heatmap_data[:, 1]', 1:1, heatmap_data[:, 2]', title="Cluster Size vs Frequency", xlabel="Cluster Size", ylabel="Frequency", color=:viridis, legend=false)

   
   
#    heatmap_data = zeros(Int, maximum(keys(cluster_size_freq)), 2)
    for (size, freq) in cluster_size_freq
        heatmap_data[size, 1] = size
        heatmap_data[size, 2] = freq
    end

#    println(heatmap_data)  
    # Plot heat map of cluster size vs frequency
#    heatmap!(subplot=3,[heatmap_data[:, 1]], [heatmap_data[:, 2]], title="Cluster Size vs Frequency", xlabel="Cluster Size", ylabel="Frequency", color=:viridis, legend=false)
    #plot(p1, p2, p3, layout=(3,1))   #, size=(600, 800))
    savefig("results.png")
end

# Call the analyze function
all_cluster_sizes, cluster_size_freq, total_clusters_per_frame = analyze_clusters(peptide_coords, box_sides, nframes, nmols, distance_threshold)

# Call the plot function
#plot_clusters(all_cluster_sizes, cluster_size_freq, total_clusters_per_frame)




















    #using Plots, StatsPlots

    #function plot_dbscan_results(coords, db_result)
        # Create a scatter plot for the points
    #    scatter(coords[1, :], coords[2, :], zcolor=db_result.assignments, label="", markerstrokewidth=0, markersize=5)
    #    xlabel!("X")
    #    ylabel!("Y")
    #    zlabel!("Z")
    #    title!("DBSCAN Clustering Results")
    #end


    #plot_dbscan_results(coords, db_result)
    #savefig("dbscan.png")
 
    # Analyze clusters for this frame
    #num_clusters = length(unique(db_result.assignments)) - (0 in db_result.assignments ? 1 : 0) 
    #println("Number of Clusters (including noise): ", num_clusters)

    # Loop through clusters
#    for (i, cluster) in enumerate(db_result)
#        if i == 1
#            println("Noise Points: ", cluster)
#        else
#            println("Cluster $i: ", cluster)
#        end
#    end

#for cluster_id in 1:num_clusters
#    points_in_cluster = findall(x -> x == cluster_id, db_result.assignments)
#    println("Cluster $cluster_id: $points_in_cluster")
#end

    
#end



