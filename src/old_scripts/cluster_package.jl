
# script using a clustering package

# 0.6 mol / L - BGL = 227 ---- AGL == 127 
include("./cluster_size.jl")

solute     = "protein"
cosolvent  = "resname BGL and resname AGL"
atoms      = PDBTools.readPDB("./10/processed.pdb") # selection - all atoms of the PDB file.

sel0    = PDBTools.select(atoms,solute)        # Selection of which atoms from the PDB are the protein atoms.
protein = Selection(sel0, nmols=15)             # Selection nmols = 1 (one protein)

#sel1    = PDBTools.select(atoms, cosolvent)        # Selection of which atoms from the PDB are the protein atoms.
#cosv    = Selection(sel0, nmols=354)             # Selection nmols = 1 (one protein)

trajectory = Trajectory("./10/processed.xtc", protein, format="XTC" )

# number of protein atoms
nprot = length(trajectory.x_solute)

# number of frames
nframes = trajectory.nframes

#for iframe in 1:nframes
  nextframe!(trajectory) # acess coordinates of the next frame 
  solute    = trajectory.solute    # selection  
  x_solute  = trajectory.x_solute  # matrix with protein coordinates
  sides = getsides(trajectory, 1)  # box sides

  # tests to define the estrategy to get coordinates for each peptide
  natoms= solute.natomspermol
  peptides = Vector{Vector{Int64}}(undef,solute.nmols) # vector to save indexes for all peptides  
  
  # function to get the indexes of a peptide    
  get_indexes(indexes::Vector{Int64}, ipep::Int64, natoms) = indexes[ (ipep-1) * natoms + 1 : ipep*natoms ]
  get_indexes(indexes::Vector{Int64}, natoms;ipep=1) = indexes[1:natoms] 

  # saving each peptide index
  peptides[1] =  get_indexes(solute.index, natoms) 
  for i in 2:solute.nmols
    peptides[i] = get_indexes(solute.index, i, natoms)
  end

  #pep1 = solute.index[1:natoms] # indexes for the first peptide    
  #pep2 = solute.index[natoms+1:natoms+natoms] # indexes for the second peptide
  #pep3 = solute.index[2*natoms+1:3*natoms]
  #pep8 = get_indexes(solute.index, 8, natoms)
  #pep15 = get_indexes(solute.index, 15, natoms) 
  
  # getting coordinates of the peptides
  coord_matrix = []   #  Vector{SVector{3, Float64}}(undef,solute.nmols)    
  for i in 1:solute.nmols
    push!(coord_matrix, x_solute[peptides[i]])  
  end
  
  #x1 = x_solute[pep1]
  #x2 = x_solute[pep2]
  #x3 = x_solute[pep3]
  #x8  = x_solute[pep8] 
  #x15 = x_solute[pep15]   
    
  
  # calculating the minimum-distance between two peptides
  
  function find_min(x1::Vector{SVector{3, Float64}}, x2::Vector{SVector{3, Float64}}, sides::SVector{3, Float64}, natoms)
    mindist = Vector{Float64}(undef, natoms) # vector to save all minimum-distances beetween the peptide1 atoms and the others peptides
    for k in 1:natoms     # loop through arbitrarily chosen peptide 1
      dist = +Inf
      for i in 1:natoms   # loop through the second peptide
        d = distance(x1[k], x2[i], sides)
        if dist > d
          dist = d
        end
      end
      mindist[k] = dist 
    end
    return minimum(mindist)
  end

  min_dist = zeros(solute.nmols, solute.nmols)
  
  for i in 1:solute.nmols
    for j in 1:solute.nmols 
      if i == j
        min_dist[i,j] = 0
      elseif j > i
        min_dist[i,j] = 0
      else
        min_dist[i,j] = find_min(coord_matrix[i], coord_matrix[j], sides, natoms)  
      end   
    end
  end

  #matrix with the minimum distances

  min_dist


##  # function to check if the distance between peptides i and j is smaller than the cutoff
##  
##  function determine(matrix1::Matrix{Float64}, cutoff)
##    new = Matrix{Any}(undef,length(matrix1[:,1]),length(matrix1[1,:]))
##    for i in 1:length(matrix1[1,:])
##      for j in 1:length(matrix1[:,1])
##        if i == j || j > i
##            new[i,j] = Nothing
##        else
##            if matrix1[i,j] < cutoff
##                new[i,j] = true
##            else
##                new[i,j] = false
##            end
##        end
##      end
##    end
##    return new
##  end
##  
##  
##  # info[i,j] = true, means that peptide i and j are part of the same cluster
##  cutoff = 15.0
##  info = determine(min_dist, cutoff)
##  
##  
##  number_neigh = zeros(Int64, natoms)
##  
##  
##  # storage = vector to save the number of neighbors for each ptptide
##  function find_neigh!(storage, data)
##    for i in 1:length(data[1,:])
##      for j in 1:length(data[:,1])
##        if data[i,j] == true    
##          storage[i] += 1
##          println("$i is in contact with $j")  
##        end
##      end
##    end
##  end 
##  
##  
##  find_neigh!(number_neigh, info) 
##  
##  # matrix to say if the clusters are connected
##  
###  nodes = Matrix{Bool}(undef,15,15)
##  
###  function get_nodes!(nodes, data)
###    for i in 1:length(data[1,:])
###      for j in 1:length(data[:,1])
###        if data[i,j] == true
###          nodes[i,j] = true
###        else 
###          nodes[i,j] = false
###        end
###      end
###    end
###  end
##  
### get_nodes!(nodes, info) 
##  
##  # aqui eu preciso pensar em um algoritmo que percorre a matriz e checa quantos clusteres tem
##  # tem que levar em conta que peptideos distantes podem estar no mesmo cluster devido ao fato de
##  # que podemos ter 1 --- 3 ---6 --- 8
##  
##  # assim o algortimo precisa ser esperto o suficiente para classificar as coisas em clusteres corretamente.
##  
###function find_clusters(nodes)
###  n_clusters    = 0
###  name_clusters = []
###  size_clusters = []
###
###  for j in 1:length(nodes[:,1]) # olhando para as colunas primeiro
###     
###    cluster = []  
###    for i in j + 1:length(nodes[1,:])
###      n = 0
###      if nodes[i,j] == 1
###        n = n + 1
###        if n == 1 
###          push!(cluster, j)
###        end 
###        push!(cluster,i)  
###      end
###    end
###   
###    if length(cluster) == 0
###      push!(name_clusters, j)
###      push!(size_clusters, 1)   
###    else
###      push!(name_clusters, cluster) 
###      push!(size_clusters, length(cluster))
###    end 
###
###  n_clusters = length(name_clusters)
###
###  end
###
###  return n_clusters, name_clusters, size_clusters
###end  
##
### n, names, size = find_clusters(nodes)
##
##Data = Dict{String, Vector{Int64}}()
##
##function add_keys!(n, data::Dict{String, Vector{Int64}})
##  for key in 1:n
##    data["$key"] = []
##  end  
##end
##        
###findall(x->x in A,B) # to find common values between two values          
##
##function get_together(x,y)
##  return [x;y]
##end         
##
##
##function find_clusters!(nodes, data::Dict{String, Vector{Int64}})
##  for col in 1:length(nodes[:,1]) - 1 # olhando para as colunas primeiro
##    n = 0 
##    for row in col + 1:length(nodes[1,:])
##      if nodes[row,col] == 1
##        n = n + 1
##        if n == 1 
##          push!(data["$col"], col)
##        end 
##        push!(data["$col"],row)
##      end
##    end
##  end
##end 
##
##                             
###  https://people.revoledu.com/kardi/tutorial/VB/tips/Symmetric-Matrix.html
##add_keys!(15,Data)
##find_clusters!(info, Data)
##
### function that check if a value is inside of a 1D array
##check(x, vector) = x in vector
##
##
### function to look through the Dictionary and find the clusters
##
##function get_clusters(data::Dict{String, Vector{Int64}})
##  empty    = []
##  clusters = []
##  for i in 1:Data.count             
##    if isempty(Data["$i"])
##      push!(empty, i)
##    else  
##      push!(clusters, Data["$i"])
##    end 
##  end
##  
### REMOVING PEPTIDES FROM THE EMPTY VECTOR THAT ARE IN SOME CLUSTERS AND WAS NOT ASSINGED DUE TO THE FACT THAT IS LABELED BY A HUGE INDEX
##
##  for i in empty
##    check_val = i 
##    for j in 1:length(clusters) 
##      if check(check_val, clusters[j])
##        println("entrou")
##        deleteat!(empty, findall(x->x==check_val, empty))
##        println("passou do delete")        
##        continue
##      end
##    end
##  end
##
##  return empty, clusters
##end
##
##
##
##empty, clusters = get_clusters(Data)
##
##  
##  # function that find clusters based on three different conditions
###  
###  function max_size_cluster(nodes)
### 
###    for i in 1:length(nodes[1,:])
###      for j in 1:length(nodes[:,1])
###        if i == j || j > i
###            new[i,j] = Nothing
###        else
###            if matrix1[i,j] < cutoff
###                new[i,j] = true
###            else
###                new[i,j] = false
###            end
###        end
###      end
###    end
###
###  
###  end
###  
###  
###  function n_clusters(nodes)
###  
###  
###  
###  
###  
###  
###  end
##
##
##
##
##
##
### min1_3  =  find_min(x1,x3, simes, natoms)
## # min1_2  =  find_min(x1,x2, sides, natoms) 
## # min2_3  =  find_min(x2,x3, sides, natoms) 
## # min1_8  =  find_min(x1,x8, sides, natoms) 
## # min1_15 =  find_min(x1,x15, sides, natoms) 
## # min8_15 =  find_min(x8,x15, sides, natoms) 
##
##  # plotting the result to evaluate the calculation
##  #using Plots
##
##  #scatter(min1_3, label="min1_3", color="red")
##  #scatter!(min1_2, label="min1_2", color="blue")
##  #scatter!(min2_3, label="min2_3", color="orange")
##  #scatter!(min1_8, label="min1_8", color="black")
##  #scatter!(min1_15, label="min1_15", color="gray") 
##  #scatter!(min8_15, label="min8_15", color="green")
##  #
##  #savefig("peptides_dist.png")
##
##
###https://julia.school/julia/arrays/
##
##
### by the position
##
###deleteat!(astronauts, 2)
##
##
##
##
##
### by name or value
###rockets = ["Apollo", "Saturn", "Falcon Heavy"]
##
###deleteat!(rockets, findall(x->x=="Saturn",rockets))
##
##
##
##
##
##  # D = distance(x1[1], x3[1], sides) 
##  #dist = +inf
##  #for i in 1:natomspermol
##  #  d = distance(x1[1], x3[i], sides)
##  #  if dist > d
##  #    dist = d
##  #  end
##  #end
##
##
##
##
##
##
##
###end   # end of the main loop though the frames





