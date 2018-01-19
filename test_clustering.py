
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree, to_tree
import numpy as np
from matplotlib import pyplot as plt

flat_list = [0.9, 0.3, 0.1, 0.1,
                   0.2, 0.5, 0.6,
                       0.95, 0.8,
                             0.99]

flat_list = [0.9, 0.2, 0.9, 0.9, 0.1, 0.9]
Z = linkage(flat_list, 'ward')
clusters = cut_tree(Z)
print clusters

tree = to_tree(Z) #, rd=True)
print tree.get_id()
# print(dir(Z))

plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
plt.show()
