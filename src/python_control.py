# script to calculate the same thing and see if it is ok!

import MDAnalysis as mda
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt

# Replace with your actual file paths
pdbfile = 'processed.pdb'
xtcfile = 'processed.xtc'

# Load MD trajectory
u = mda.Universe(pdbfile, xtcfile)

# Assuming nmols is the number of peptides
nmols = 15
distance_threshold = 3.5  # In Angstroms

# Analyze each frame
for ts in u.trajectory:
    print(f"Frame: {ts.frame}")

    # Calculate centers of mass of each peptide
    coms = np.array([u.select_atoms(f"resid {i}:{i+nmols-1}").center_of_mass() for i in range(1, len(u.residues)+1, nmols)])

    # Compute pairwise distances
    distance_matrix = squareform(pdist(coms))

    # Apply DBSCAN clustering
    db = DBSCAN(eps=distance_threshold, min_samples=2, metric='precomputed').fit(distance_matrix)
    labels = db.labels_
    
    # Number of clusters in labels, ignoring noise if present.
    n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
    print(f"Estimated number of clusters: {n_clusters}")

    # Plotting (optional)
    unique_labels = set(labels)
    colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = [0, 0, 0, 1]

        class_member_mask = (labels == k)

        xy = coms[class_member_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col), markeredgecolor='k', markersize=14)

    plt.title('Estimated number of clusters: %d' % n_clusters)
    plt.show()

