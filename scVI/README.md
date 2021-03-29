# scVI
Processing of single-cell data: 
## Myoblasts / AMPs

![alt text](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/scVI/AMP_sc.jpg?raw=true)

> Data subsets were harmonized into collective AMP or epithelium datasets using scVI. The scVI VAE model consisted of 2 layers (n_layers=2) and 20 latent dimensions (n_latent=20), with a negative-binomial reconstruction loss (reconstruction_loss=‘nb’). The model was trained on variable genes selected by Seurat’s variance-stabilizing transformation method; 1,000 (for epithelial subsets) or 2,000 (for AMP subsets) variable genes were calculated for each inputted batch, and then the union of these genes was supplied to scVI. The following parameters were used for model training: train_size=0.75, n_epochs=400, and lr=1e-3 (other parameters were left as default). Cell clustering and UMAP was performed using [Seurat](https://satijalab.org/seurat/index.html) on the latent space derived from the scVI model.

* [Running scVI](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/scVI/running_scVI.py)
* [Correcting for cell sex and cell-cycle](https://github.com/HariharanLab/Everetts_Worley_Yasutomi/blob/master/scVI/cell_sex_and_cell_cycle_correction.R)
