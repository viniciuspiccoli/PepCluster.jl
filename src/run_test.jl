using PepCluster, Plots, Plots.Measures, LaTeXStrings

Data = data_info()

max_sizes, total_numbers, cluster_size_freq , natoms_per_peptide = cluster_dbscan(Data) 

# tst for plot
cluster_size_freq_normalized = cluster_size_freq ./ sum(cluster_size_freq)
min_freq, max_freq = extrema(cluster_size_freq_normalized);
tN_max = maximum(max_sizes) / natoms_per_peptide 
sz_max = maximum(total_numbers)

# color
custom_gradient = cgrad([:white, :lightpink, :red])
#c:=plasma


# heat map plot
heatmap(cluster_size_freq_normalized',
               c=custom_gradient,
               ylabel="Cluster Size",
               xlabel="Number of clusters",
               title="Cluster Size Distribution",
               color=:auto,
               xlims=(0,tN_max),
               ylims=(0,sz_max),  
               clims=(min_freq, max_freq),
               colorbar_title="Frequency")

savefig("cluster_dbscan.png")