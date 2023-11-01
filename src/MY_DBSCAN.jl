using Distances

mutable struct Point
    coords::Vector{Float64}
    cluster::Int  # -1 for noise, 0 for unassigned
end

function dbscan(points::Vector{Point}, eps::Float64, min_points::Int)
    cluster_id = 0
    for point in points
        if point.cluster == 0  # unassigned
            if expand_cluster(points, point, cluster_id, eps, min_points)
                cluster_id += 1
            end
        end
    end
end

function expand_cluster(points::Vector{Point}, point::Point, cluster_id::Int, eps::Float64, min_points::Int)
    seeds = region_query(points, point, eps)
    if length(seeds) < min_points
        point.cluster = -1  # mark as noise
        return false
    else
        point.cluster = cluster_id  # assign to cluster
        setdiff!(seeds, [point])  # remove point from seeds
        
        while !isempty(seeds)
            current_point = popfirst!(seeds)
            if current_point.cluster <= 0  # unassigned or noise
                if current_point.cluster == 0  # unassigned
                    push!(seeds, current_point)  # add to seeds if it's not already there
                end
                current_point.cluster = cluster_id  # assign to cluster
            end
        end
        
        return true
    end
end

function region_query(points::Vector{Point}, query_point::Point, eps::Float64)
    neighbors = []
    for point in points
        if euclidean(point.coords, query_point.coords) < eps
            push!(neighbors, point)
        end
    end
    return neighbors
end

function run_example()
    # Example data points
    data = [1 1; 2 2; 3 3; 8 8; 9 9; 25 25]
    points = [Point(data[i, :], 0) for i in 1:size(data, 1)]
    
    # Run DBSCAN
    dbscan(points, 3.0, 2)
    
    # Show results
    for (i, point) in enumerate(points)
        println("Point $i: coords = $(point.coords), cluster = $(point.cluster)")
    end
end

run_example()
