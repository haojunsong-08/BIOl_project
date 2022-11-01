# BIOl_project
## when you first time create an conda environment
Log on your pace-ice.pace.gatech.edu and enter your login and password

If you want to work on Jupyter NoteBook
```
git clone https://github.com/haojunsong-08/BIOl_project.git && cd BIOL_project
module load anaconda3/2021.05
conda create --name project r-base jupyter
conda activate project
conda install -c r r-irkernel

```
## From now on, You can access the package can submit PBS script with the following at the beginning
By actiavting conda environment
```
module load anaconda3/2021.05
conda activate project
```
