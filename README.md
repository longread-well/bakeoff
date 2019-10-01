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

# A possible way of working

Suppose you are conducting analysis in `${BAKEOFF_ROOT}/your/folder`.  I suggest the following workflow (given as an example for snakemake):

* Clone this repo as the scripts/ folder.
* symlink the relevant snakefile: `ln -s scripts/[...]/conduct_analysis.snakemake ./Snakefile.symlink`
* Now run your analysis as `snakemake -s Snakefile.symlink`

(If you're not using snakemake, an analogous process no doubt applies using e.g. a shell or python script.)

The advantages of this are:

* we get to keep all scripts / snakefiles in one place (this github repo), so we can see them all.
* but you also have the code right there in the folder you're analysing (including any code shared between analyses).
* any changes you make can be tracked using git.

