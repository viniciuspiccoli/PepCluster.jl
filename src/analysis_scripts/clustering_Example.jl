using Clustering, Plots, Random, Distributions, LinearAlgebra

# naive
Random.seed!(123)
points = randn(3, 1000)
#clusters = dbscan(points, 0.05, min_neighbors = 1, min_cluster_size = 3)

clusters = dbscan(points, 0.05, min_neighbors=10)

labels = assignments(clusters)
n_clusters = maximum(labels)
plot(legend = false, xlabel="Dimension 1", ylabel="Dimension 2", zlabel="Dimension 3")

# Plot each cluster with a different color
for i in 1:n_clusters
    cluster_points = points[:, labels .== i]
    scatter!(cluster_points[1, :], cluster_points[2, :], cluster_points[3, :], 
             marker_z=cluster_points[3, :], label="Cluster $i")
end

# Plot noise points
noise_points = points[:, labels .== 0]
scatter!(noise_points[1, :], noise_points[2, :], noise_points[3, :], color=:grey, 
         alpha=0.5, label="Noise")
savefig("example.png")