# BIOl_project
## Please Ensure the following dependencies are installed
```
install.packages(c("shiny", "shinythemes", "Seurat", "ggplot2", "patchwork", "dplyr", "hdf5r", "visNetwork", "viridis", "ConsensusClusterPlus", "devtools", "ggpubr" , "pheatmap"))
library(devtools)
install_github("navinlabcode/CellTrek")
```
## Please also make sure the data is stored in the appropriate folder and create a folder for the rds data
# Especially the Anterior and Rds_data folder , so there will be no need for changing the script
![Screenshot 2022-11-29 at 21 49 06](https://user-images.githubusercontent.com/84302343/204695777-bba5cc46-ccf5-454f-9e32-897a9b63cdb0.png)

## when you first time create an conda environment
Log on your pace-ice.pace.gatech.edu and enter your login and password

If you want to work on Jupyter NoteBook
```
git clone https://github.com/haojunsong-08/BIOl_project.git && cd BIOL_project
module load anaconda3/2021.05
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create --name project_seurat r-base jupyter python=3.6
conda activate project_seurat
conda install r-seurat
conda install -c conda-forge r-irkernel
conda install -c conda-forge r-hdf5r

```
<img width="1206" alt="Screen Shot 2022-11-01 at 13 39 14" src="https://user-images.githubusercontent.com/84302343/199300837-d660be34-b85a-4c72-b964-875597f77bda.png">

## From now on, You can access the package can submit PBS script with the following at the beginning
By actiavting conda environment
```
module load anaconda3/2021.05
conda activate project_seurat
```
