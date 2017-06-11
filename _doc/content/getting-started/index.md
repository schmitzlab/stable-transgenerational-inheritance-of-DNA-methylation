+++
date = "2017-06-02T00:11:02+01:00"
title = "Getting started"
weight = 20
+++

## Installation
 To install, either clone this repositiory or download the source code.
 
 ```
    # Using git clone
    git clone https://github.com/schmitzlab/stable-transgenerational-inheritance-of-DNA-methylation.git
 ```

## Set-Up
 All python scripts were created for Phyton 3.4+ and specifically tested with Python 3.5.2. All scripts import `sys` and `os`. Many import `math`, `glob`, `bisect`, and `random`. Several import `pandas`, `numpy`, `scipy`, `functools`, and `sklearn`. Any script ending in `_pe.py` has the ability to be run using multiple processors and imports `multiprocessing` and `subprocess`. See [appendix](/appendix) to determine which scripts need which packages
 
 All required imports, as well as Python itself, are available via [anaconda](https://www.continuum.io/downloads). This the preferred way to install package dependencies. See anaconda for instructions on how to install and how to download required packages.
 
## Software Versions
 
Listed below are the software versions used for analysis in the paper.

- Python: 3.5.2
- Anaconda: 4.1.6
    - Numpy: 1.11.0
    - Scipy: 0.17.1
	- Scikit: 0.17.1
	- Pandas: 0.17.0

Most scripts should work with other versions, but it is not guaranteed. Scripts will **not*
* work for Python2.

## Required Scripts

These scripts are required by several of the other scripts and contain utility functions. These must be in the same folder as the other scripts or on the python path. 

To add to python path, `export PYTHONPATH=$PYTHONPATH:/path/to/scripts`.

* *bioFiles.py*: for reading common file types
* *bth_util.py*: extra utility functions such as converting an integer to a human-readable size
