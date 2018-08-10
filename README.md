# massoc

A platform for inference and analysis of microbial association networks. V0.1.0.

## Getting Started

To run <i>massoc</i>, you only need to run the appropriate executable for your platform.
These executables are stand-alone and do not need the other files in this repository.
<a href="https://github.com/ramellose/massoc/releases">Find the latest releases here</a>.

## Development version and command line

For development purposes, the GUI can be executed from run_massoc.py.
The command line interface is accessible through main.py, but I am working on a stand-alone version that does not require you to manage dependencies.

massoc runs on Python >= 3.5. All its dependencies are listed in the requirements.txt file.

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
Extract the files from the distribution and run <b>bin\neo4j.bat console</b> to initialize the server.
On Linux, set your default Java version to 1.8.0:
sudo apt install openjdk-8-jk
sudo update-java-alternatives --set java-1.8.0-openjdk-amd64

### Tutorials

Check out the demo at https://github.com/ramellose/massoc/raw/master/massoc/docs/massoc_demo.pdf

### Contributions

This software is still in early alpha. Any feedback or bug reports will be much appreciated!

## Authors

* **Lisa Röttjers** - [ramellose](https://github.com/ramellos)

See also the list of [contributors](https://github.com/ramellose/massoc/contributors) who participated in this project.

## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

massoc is largely built on and inspired by the [BIOM format](http://biom-format.org/).


