# massoc

A platform for inference and analysis of microbial association networks.

## Getting Started

Simply download the repository and run the GUI through run_massoc.py, or the command interface with main.py.

## Prerequisites

massoc runs on Python >= 3.5. All its dependencies are listed in the requirements.txt file.
Install them with:
```
pip install -r requirements.txt
```

To run network inference, you will need to ensure that your system can run your chosen network inference tools.
Copy the CoNet files and SparCC files in the execs folder and massoc will run them for you.

We cannot guarantee that the network inference tools will work for you, but the tips below should help you get started.

### Installing SPIEC-EASI

To run SPIEC-EASI from massoc, you will need the following:
* Rscript
* biomformat
* docopt
* SpiecEasi

Install Rscript on Ubuntu:
```
sudo apt-get littler
```
For Windows, add the location of Rscript.exe to PATH if you cannot call Rscript from command line.

The required libraries also need to be installed in your R environment:
```
source("https://bioconductor.org/biocLite.R")
biocLite("biomformat")
library(devtools)
install_github("zdk123/SpiecEasi")
install_packages("docopt")
```

### Installing CoNet

Download the CoNet file from [here](http://psbweb05.psb.ugent.be/conet/download.php).
Unzip it in the massoc/execs folder. To run CoNet, the Java Runtime environment will need to be installed.

### Installing SparCC

Download the SparCC file from [here](https://bitbucket.org/yonatanf/sparcc) and unzip it in the massoc/execs folder.
To run SparCC, you will need to have Python 2.6 or 2.7 installed.
On Windows, it will be necessary to set your PATH so both installations can be found.
Find your Python 2 installation, rename the python.exe file to python2.exe and add it to your PATH.

### Tutorials

The massoc documentation contains a demo. Give the example BIOM files in the data folder a try!

### Contributions

This software is still in early alpha. Any feedback or bug reports will be much appreciated!

## Authors

* **Lisa Röttjers** - [ramellose](https://github.com/ramellos)

See also the list of [contributors](https://github.com/ramellose/massoc/contributors) who participated in this project.

## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

massoc is largely built on and inspired by the [BIOM format](http://biom-format.org/).


