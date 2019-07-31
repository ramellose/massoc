# massoc

A platform for inference and analysis of microbial association networks. V0.1.0.

Currently, the following features are available:
* Preprocessing of count files for network inference
* Clustering + splitting files by sample properties
* Batch network inference
* Network + BIOM file storage in a graph database
* Multi-network logic operations
* Taxonomy-dependent edge agglomeration
More features are currently under construction.

Welcome to <i>massoc</i>! Contact the author at lisa.rottjers (at) kuleuven.be. Your feedback is much appreciated!
This is version 0.3.0, and therefore still in early alpha. Encountering bugs is highly likely!

## Getting Started

You can access <i>massoc</i> through a CLI and a GUI.
Detailed documentation of these interfaces, combined with a  manual describing all inputs, will be available soon.

## Installation

To run <i>massoc</i> on Windows or Ubuntu, you only need to run the appropriate executable for your platform.
These executables are stand-alone and do not need the other files in the repository. However, some of the features
require you to have installed specific software, for example Neo4j Server. Read the sections below to figure out which
dependencies you need to install and how to do so.
<a href="https://github.com/ramellose/massoc/releases">Find the latest releases here</a>.

If you do not want to use these executables or have macOS system, you can install <i>massoc</i> directly from the repository.
This can result in some issues with dependencies.
```
pip install git+https://github.com/ramellose/massoc.git
```

There is no executable available for macOS systems, but you can set up your macOS system to run <i>massoc</i> as follows.
First install macports, then run:
```
sudo port install python35
sudo port select --set python python35
sudo port install py35-pip
sudo port select --set pip pip35
pip install numpy==1.13.1 scipy h5py nose Cython psutil neo4j matplotlib networkx multiprocess sklearn scikit-bio biom-format --user
pip install -U wxPython --user
```

To reset to the default Python version, run:
```
sudo port select --set python python37
sudo port select --set pip pip37
```

## Network inference in <i>massoc</i>

To run network inference, you will need to ensure that your system can run your chosen network inference tools.
If you want to run CoNet, SparCC or SPIEC-EASI from <i>massoc</i>, the instructions below will help you get started. 

### Running SPIEC-EASI

To run SPIEC-EASI from massoc, you will need the following:
* Rscript
* biomformat
* docopt
* SpiecEasi

The required libraries can be installed in your R environment with the script below.
```
source("https://bioconductor.org/biocLite.R")
biocLite("biomformat")
library(devtools)
install_github("zdk123/SpiecEasi")
install_packages("docopt")
```

### Running CoNet

Download the CoNet file from [here](http://psbweb05.psb.ugent.be/conet/download.php).
<i>massoc</i> will ask you for the location of the CoNet3 folder. 

### Running SparCC

Download the SparCC file from [here](https://bitbucket.org/yonatanf/sparcc).
<i>massoc</i> will ask you for the location of the SparCC folder. 
To run SparCC, you will need to have Python 2.6 or 2.7 installed.

On Windows, it will be necessary to set your PATH so both installations can be found.
Find your Python 2 installation, rename the python.exe file to python2.exe and add it to your PATH.

### Neo4j

The graph database utilities are built on Neo4j.
Download the Neo4j Server community edition from [here](https://neo4j.com/download-center/#releases).

Extract the files from the distribution.
On Windows,  run <b>bin\neo4j.bat console</b> to initialize the server.

On Linux, first set your default Java version to 1.8.0:
```
sudo apt install openjdk-8-jk
sudo update-java-alternatives --set java-1.8.0-openjdk-amd64
```
For macOS, first install the appropriate Java version:
```
brew cask install java8
```
After downloading and extracting the Neo4j Server folder, initialize the server with <b>bin\neo4j console</b>.


### Tutorials

Check out the demo at https://github.com/ramellose/massoc/raw/master/massoc/docs/massoc_demo.pdf

## Development version and command line

For development purposes, the GUI can be executed from massocGUI.py.

massoc runs on Python 3.5. All its dependencies are listed in the requirements.txt file.

### Contributions

This software is still in early alpha. Any feedback or bug reports will be much appreciated!

## Authors

* **Lisa Röttjers** - [ramellose](https://github.com/ramellos)

See also the list of [contributors](https://github.com/ramellose/massoc/contributors) who participated in this project.

## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

massoc is largely built on and inspired by the [BIOM format](http://biom-format.org/).


