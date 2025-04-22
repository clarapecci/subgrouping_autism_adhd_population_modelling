# Subgrouping autism and ADHD based on structural MRI population modelling centiles

This repository contains an implementation to cluster MRI population modelling centiles using k-medoids. Prior to clustering, the dimensionality of the data is reduced using UMAP (or tSNE). 
To use the code, you will need a csv file containing the participants ID, the features to be used, and a 'group' column where you specify whether the participant belongs to the control group (labelled -1), or not (labelled 1). 

## Repository structure
* `alluvial.R`: Code to plot alluvials comparing the distribution of participants across different solutions of the same technique, or across techniques.
* `medoids_clustering.R`: Code to cluster the data using k-medoids with a dimensionality reduction technique (UMAP or tSNE). Outputs 1) the participant assigments for k=2 to k-10 clusters, 2) the silhouette score for each solution, 3) a plot of the embeddings coloured in based on the optimal clustering solution.
* `preprocessing.R`: Code to preprocess the data including what diagnosis to use, select participants within given age window.
* `significance_analysis.R` / `significance_kmedoids.R`: This code can be used post-clustering to compare the different subgroups against the controls for the features required. The first script is specifically for HYDRA results, and the second for k-medoids.
* `util.R`: Contains utility functions used in the other scripts. 
