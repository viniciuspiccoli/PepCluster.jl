# script test - requires further optm.
using PepCluster, Plots, Plots.Measures, LaTeXStrings, ColorSchemes

#plot headings
plot_font="Computer Modern"
default(legend=:topright,
       fontfamily=plot_font,
       linewidth=2, framestyle=:box,
       grid=false)
scalefontsizes(1.2)
gr(dpi=450)
plot(layout=(2,2))

# color
custom_gradient = cgrad([:white, :lightpink, :red])



# data struct for water
Data_wt = data_info(
                     solute_sel="protein",
                     pdbfile="/home/viniciusp/Documents/GLUCOSE/water/1/protein_only.pdb",
                     xtcfile="/home/viniciusp/Documents/GLUCOSE/water/1/protein_only.xtc",
                     nmols=15,
                     distance_threshold=3.5)


max_sizes, total_numbers, cluster_size_freq , natoms_per_peptide = cluster_dbscan(Data_wt) 

# tst for plot
cluster_size_freq_normalized = cluster_size_freq ./ sum(cluster_size_freq)
min_freq, max_freq = extrema(cluster_size_freq_normalized);
tN_max = maximum(max_sizes) / natoms_per_peptide 
sz_max = maximum(total_numbers)


# heat map plot for water
heatmap!(cluster_size_freq_normalized',
               c=custom_gradient,
               xlabel="Number of Clusters",
               ylabel="Cluster Size",
              # title="Cluster Size Distribution",
               color=:auto,
               xlims=(0,15),
               ylims=(0,15),  
               clims=(0, 0.1),
               colorbar_title="Frequency",
               subplot=1)

annotate!( 10, 10, text(L"\textrm{Water}", :center, 10), subplot=1)

plot!(subplot=3, total_numbers, xlabel="Frame", ylabel="NC")



# data for glucose solutions
Data_gl = data_info(
                     solute_sel="protein",
                     pdbfile="/home/viniciusp/Documents/GLUCOSE/gl_1M/01/protein_only.pdb",
                     xtcfile="/home/viniciusp/Documents/GLUCOSE/gl_1M/01/protein_only.xtc",
                     nmols=15,
                     distance_threshold=3.5)


max_sizes, total_numbers, cluster_size_freq , natoms_per_peptide = cluster_dbscan(Data_gl) 

# tst for plot
cluster_size_freq_normalized = cluster_size_freq ./ sum(cluster_size_freq)
min_freq, max_freq = extrema(cluster_size_freq_normalized);
tN_max = maximum(max_sizes) / natoms_per_peptide 
sz_max = maximum(total_numbers)



# heat map plot for glucose
heatmap!(cluster_size_freq_normalized',
               c=custom_gradient,
               ylabel="Cluster Size",
               xlabel="Number of clusters",
             #  title="Cluster Size Distribution",
               color=:auto,
               xlims=(0,15),
               ylims=(0,15),  
               clims=(0, 0.1),
               colorbar_title="Frequency",
               subplot=2)

annotate!( 10, 10, text(L"\textrm{Glucose - 1.0}\ \mathrm{mol\ L^{-1}}", :center, 10), subplot=2)
plot!(subplot=4, total_numbers, xlabel="Frame", ylabel="NC")


plot!(margin=(5mm), size=(800,500))
savefig("Cluster_water_gl1M.png")