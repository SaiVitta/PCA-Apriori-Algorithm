"""
Description: Generate scatter plots using PCA, SVD and TSNE algorithms for a given dataset 
Group Members: Mitali Bhiwande | Sumedh Ambokar | Tejasvi Sankhe
"""

import sys
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE


# Function to perform PCA algorithm on the given data set.
def plot_pca(original_attrs, disease_array, filename):
    attrs_mean = original_attrs.mean(axis=0)
    adjusted_attrs = original_attrs - attrs_mean
    covariance = np.dot(np.transpose(adjusted_attrs),adjusted_attrs)/len(adjusted_attrs)
    w, v = LA.eig(covariance)
    top_eigen_vectors = v[:,0:2]
    new_coordinates = np.dot(original_attrs,top_eigen_vectors)
    draw_scatter_plot(new_coordinates[:,0:1],new_coordinates[:,1:2],disease_array,'pca', filename)

	
# Function to perform SVD on the given data set using numpy linalg SVD library
def plot_svd(original_attrs, disease_array, filename):
    U, s, V = LA.svd(original_attrs, full_matrices=True)
    newdata = U[:,:2]
    draw_scatter_plot(newdata[:,0:1],newdata[:,1:2], disease_array,'svd', filename)


# Function to perform TSNE on the given data set using sklearns manifold TSNE library
def plot_tsne(original_attrs, disease_array, filename):
    tsne = TSNE(n_components=2, init='pca', n_iter=1000, learning_rate=100)
    attr_tsne = tsne.fit_transform(original_attrs)
    draw_scatter_plot(attr_tsne[:,0:1],attr_tsne[:,1:2], disease_array,'tsne', filename)


# Function to differentiate attributes and diseases from the given data file
def fetch_attributes(original_data):
    original_attrs = original_data[:,0:original_data.shape[1]-1]
    original_attrs = np.array(original_attrs, dtype=float)
    disease_array = original_data[:,original_data.shape[1]-1:original_data.shape[1]]
    return original_attrs, disease_array


# Function to generate scatter plots for passed data parameters
def draw_scatter_plot(x,y,disease_array,plot_type,filename):
    categories = np.unique(disease_array)
    colors = np.linspace(0, 1, len(categories))
    use_colours = dict(zip(categories, colors))
    unique = list(use_colours.values())
    plot_size = plt.rcParams["figure.figsize"]
    plot_size[0] = 12
    plot_size[1] = 9
    plt.rcParams["figure.figsize"] = plot_size
    colors = [plt.cm.jet(float(i)/max(unique)) for i in unique]
    for u in use_colours.keys():
        xi = [x[j] for j  in range(len(x)) if disease_array[j] == u]
        yi = [y[j] for j  in range(len(x)) if disease_array[j] == u]
        plt.scatter(xi, yi, c=colors[unique.index(use_colours.get(u))], label=str(u))
    if plot_type=='pca':
        plt.xlabel("PC 1")
        plt.ylabel("PC 2")
        title = "PCA: "+filename
        plt.title(title)
    elif plot_type=='svd':
        plt.xlabel("Component 1")
        plt.ylabel("Component 2")
        title = "SVD: "+filename
        plt.title(title)
    elif plot_type=='tsne':
        title = "TSNE: "+filename
        plt.title(title)
    plt.legend()
    plt.show()


# Funtion to call individual functions for generating plots of PCA, SVD and TSNE
def generate_plots(filename):
    with open(filename) as textFile:
        lines = [line.replace("\n","").split("\t") for line in textFile]
    original_attrs, disease_array = fetch_attributes(np.array(lines));
    print("PCA FOR FILE: ",filename)
    plot_pca(original_attrs, disease_array, filename)
    print("SVD FOR FILE: ",filename)
    plot_svd(original_attrs, disease_array, filename)
    print("TSNE FOR FILE: ",filename)
    plot_tsne(original_attrs, disease_array, filename)

# Reads the data file path as the command line input
if len(sys.argv) > 1:
	generate_plots(sys.argv[1])
else:
	print("Not sufficient Inputs")



