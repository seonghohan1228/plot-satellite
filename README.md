# Steps
1. Install Python 3
2. Install Visual Studio Code (Add Python extension after installation)
3. Install Anaconda3

# After Visual Studio Code installation
In 'Extension' in the left panel, install the Python extension.

# After Anaconda3 installation (Windows)
In
```
Control Panel > System and Security > System > Advanced system settings > Environmetn variables > System variables > find 'Path' > New > Browse
```
add path to
```
path/to/anaconda3
path/to/anaconda3/Scripts
path/to/anaconda3/Library/bin
```

# Packages
### basemap
In VS Code terminal, type
```
conda init
conda create -p ./venv python=<python-version>
conda activate ./venv
conda install -c conda-forge proj4
conda install -c conda-forge basemap
```
For IGRF, type
```
pip install IGRF
```


