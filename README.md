# Anaconda
  Install: https://www.anaconda.com/download

# .condarc
  Change the channel priority to avoid problems per https://conda-forge.org/docs/conda-forge_gotchas.html using

  `conda config --prepend channels conda-forge`

# NB Conda:
  In order to use an environment's kernel in jupyter, you must install this package in your global environment.

  Install: `conda install nb_conda`

# Virtual Environment:
  Create: `conda create --name myenv`

  or

  Create with GeoPandas: `conda create -n myenv python=3.6 geopandas -c conda-forge`

  or

  Create Environment Complete with Packages: `conda create --name myenv --file spec-file.txt` (created with `conda list --explicit > spec-file.txt`)

  Activate: `source activate myenv`

  Deactivate: `source deactivate`

# Package Info:
## NB Conda (Allows running jupyter notebook in an environment):
  Install: `conda install nb_conda`

## Jupyter Notebook (needs to be installed in environment):
  Install: `conda install jupyter`

## GeoPandas:
  Install: `conda install -c conda-forge geopandas`

  Reference: http://geopandas.org/index.html

## ArcGIS API for Python:
  Install: `conda install -c esri arcgis`

  Reference: http://esri.github.io/arcgis-python-api/apidoc/html/index.html
