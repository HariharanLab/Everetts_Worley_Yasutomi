# -*- coding: utf-8 -*-
"""
Example Code for Running scVI (Ver. 0.4.1)

Author: Nicholas Everetts

NOTE: Since this publication, scVI has been updated significantly, and is now 
called scvi-tools). I recommend learning how to run the newest version of
scvi-tools, although this code can help understand what inputs you will need 
and parameters to consider. scVI version 0.4.1 contained functions for loading
data from local files, but these functions have since been removed from 
scvi-tools, which now utilizes Scanpy and its AnnData object.
"""

#%%
#Define paths to files

#variable_genes_file = path_to_seurat_variable_genes
#.csv file with two columns, corresponding to the gene names and their
#corresponding indices in the dataset.
#NOTE: Be sure to account for zero-based indexing in Python and one-based indexing in R.

#dataset_file = path_to_dataset
#.csv file of the expression count matrix, with genes as rows and cells as columns
#include gene names and cell barcodes.

#batch_id_file
#.csv file with a single column corresponding to batches.

#%%
#Loading important packages
import torch
import csv
import numpy as np
from scvi.dataset import CsvDataset
from scvi.models import *
from scvi.inference import UnsupervisedTrainer

#%%
#Parse the variable_genes_file to assemble two lists corresponding to variable
#gene names and their indices in the dataset.

seurat_gene_list_num = []
seurat_gene_list_names = []

with open(variable_genes_file) as csvfile:
    readCSV = csv.reader(csvfile, delimiter = ",")
    for row in readCSV:
        seurat_gene_list_num.append(int(row[1]))
        seurat_gene_list_names.append(row[0])
        
print("Number of variable genes: " + str(len(seurat_gene_list_names)))

#%%
#Loading the expression data from dataset_file, and subset the data by the
#variable genes provided in variable_genes_file.

dataset = CsvDataset(filename = dataset_file, gene_by_cell = True, save_path="")
dataset.subsample_genes(subset_genes = seurat_gene_list_num)

#%%
#Load batch identities from the batch_id_file
batch_identities = np.loadtxt(batch_id_file, delimiter=",", dtype = np.float64)
batch_identities = np.reshape(batch_identities, (len(batch_identities), 1))
print("Batch Identities: " + str(np.unique(batch_identities)))

#Provide batch information to the scVI model
dataset.n_batches = num_of_batchs
dataset.batch_indices = batch_identities

#%%
#Define the output files
#key outputs from training the scVI model include the model itself, the 
#latent space (used for clustering and UMAP), and the imputed expression matrix 
#(which we used for DistMap). The "scVI_gene_names.csv" should match the genes
#in the variable_genes_file

directory_path = "./" #folder location where you want to save scVI outputs
scVI_model_save_file = directory_path + "scVI_model.pkl"
latent_save_file = directory_path + "scVI_latent.csv"
imputation_save_file = directory_path + "scVI_imputed_exprs_mat.csv"
gene_names_save_file = directory_path + "scVI_gene_names.csv"
print(scVI_model_save_file,
      latent_save_file,
      imputation_save_file,
      gene_names_save_file,
      sep = "\n")

#%%
#Set parameters for model training

#below are the parameters used in our Everetts, Worley, et al. 2021 paper.
n_epochs = 400
lr = 1e-3

#For use_cuda: If you have an NVIDIA graphics card and have installed the
#NVIDIA CUDA toolkit, set this as True. Otherwise, set this as False
use_cuda = True
num_of_layers = 2

#For reconstrution loss: Use "nb" for negative binomial gene modeling, 
#or "zinb" to include zero-inflation component (if you think your data is
#zero-inflated)
reconstruction_loss = "nb"
train_size = 0.75

#%%
#Initialize the scVI model
vae = VAE(n_input = dataset.nb_genes,
          n_batch = dataset.n_batches * batch_boolean,
          n_latent = latent_dims,
          n_layers = num_of_layers,
          reconstruction_loss = reconstruction_loss)

#Train the scVI model
#depending on the size of your data and if you have an NVIDIA GPU, this could
#take 10 minutes to 1+ hours. If you'd like to make some tea or coffee, now
#would be an appropriate time to do so.
trainer = UnsupervisedTrainer(vae,
                              dataset,
                              train_size = train_size,
                              use_cuda = use_cuda,
                              frequency = 5)
trainer.train(n_epochs = n_epochs,
              lr = lr)
print("Model training finished!")

#Create the posterior representation of the data, and extract the latent space
#and imputed data
downsampled_gene_names = dataset.gene_names
full_posterior = trainer.create_posterior(vae, dataset, indices=np.arange(len(dataset)))
scVI_latent = full_posterior.sequential().get_latent()[0]
scVI_imputed = full_posterior.sequential().imputation()

#Save the relevant output files
np.savetxt(latent_save_file, scVI_latent, fmt='%s', delimiter = ",")
np.savetxt(imputation_save_file, scVI_imputed, fmt='%s', delimiter = ",")
np.savetxt(gene_names_save_file, downsampled_gene_names, fmt='%s', delimiter = ",")
torch.save(trainer.model.state_dict(), scVI_model_save_file)
