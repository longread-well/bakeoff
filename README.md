# Bakeoff
This repository contains code for the long-read bakeoff

# Contributing

This repository is publicly viewable!  Please include here any code or docs relevant to the project.  But please do not include any non-public project data, or other project-sensitive information.

# Using this code

I suggest using an environment variable to specify the root of the bakeoff project folder in scripts:
```
BAKEOFF_ROOT=<path to root folder of bakeoff dir>
```
This can be set (for example) in `~/.bashrc` to make it always available.

This has the advantages that a) we don't need to store absolute paths in scripts, and b) it's useful to run the same scripts on other folders, e.g. test folders, and c) it's useful to port code to other machines where the mount point differs. 

To use an environment variable in R:
```R
BAKEOFF_ROOT = Sys.getenv( "BAKEOFF_ROOT" )
```
or in python:
```python
import os
BAKEOFF_ROOT = os.environ[ 'BAKEOFF_ROOT' ]
```

