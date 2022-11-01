# BIOl_project
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

```
<img width="1206" alt="Screen Shot 2022-11-01 at 13 39 14" src="https://user-images.githubusercontent.com/84302343/199300837-d660be34-b85a-4c72-b964-875597f77bda.png">

## From now on, You can access the package can submit PBS script with the following at the beginning
By actiavting conda environment
```
module load anaconda3/2021.05
conda activate project_seurat
```
