# massoc

A platform for inference and analysis of microbial association networks.

## Getting Started

To run <i>massoc</i>, you only need to run the appropriate executable for your platform.
These executables are stand-alone and do not need the other files in this repository.

If you run into any bugs or the program appears to be stuck, the massoc.log file should tell you what went wrong.
You can find this file in your specified user directory.

## Development version and command line

For development purposes, the GUI can be executed from run_massoc.py.
The command line interface is accessible through main.py, but we are working on a stand-alone version that does not require you to manage dependencies.

massoc runs on Python >= 3.5. All its dependencies are listed in the requirements.txt file.
Install them with:
```
pip install -r requirements.txt
```

## Network inference in <i>massoc</i>

To run network inference, you will need to ensure that your system can run your chosen network inference tools.
If you want to run CoNet, SparCC or SPIEC-EASI from <i>massoc</i>, the instructions below will help you get started. 
However, the next update will also let you import networks from separate files that were generated outside <i>massoc</i>.

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

### Tutorials

The <i>massoc</i> documentation contains a demo. Give the example BIOM files in the data folder a try!

### Contributions

This software is still in early alpha. Any feedback or bug reports will be much appreciated!

## Authors

* **Lisa Röttjers** - [ramellose](https://github.com/ramellos)

See also the list of [contributors](https://github.com/ramellose/massoc/contributors) who participated in this project.

## License

This project is licensed under the Apache License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

massoc is largely built on and inspired by the [BIOM format](http://biom-format.org/).


