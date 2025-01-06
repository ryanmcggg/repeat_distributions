# McGinty, Balick, Lyskova, Mirkin and Sunyaev – Supplementary Code - README

## For Computational model:

System Requirements:

Miniconda installation of Python 3.10.14	64bit   
			              	- https://docs.anaconda.com/miniconda/
Additional Python libraries:  
	Biopython SeqIO			- https://biopython.org/wiki/SeqIO
  FastaPy             - https://github.com/aziele/fastapy
	Plotly			      	- https://plotly.com/python/getting-started/  
	Liftover	      		- https://pypi.org/project/liftover/
  Pandas              - https://pandas.pydata.org/
  Numpy               - https://numpy.org/
  Scipy.stats         - https://docs.scipy.org/doc/scipy/reference/stats.html
  
Tested on the following system:
	Windows 11  
	128 Gb RAM  
	AMD 5800X processor  
	Python 3.10.14; pandas 2.2.2; numpy 1.26.4; biopython 1.83;  
	plotly 5.22.0; liftover 1.2.2; scipy 1.13.1  

Installation guide:

Follow above instructions to install Python and additional dependencies.
Typical install time: ~15 minutes

Instructions for use:

Open Jupyter Notebooks in listed order.
Further instructions are included within the notebook, including instructions for downloading required public datasets and creating directories.
Run each notebook cell in order.

Scripts provided for running the full grid of computational models on a Slurm-based cluster.
Run time will vary based on cluster performance.

Expected output:
All plots (.pdf images) included in the manuscript.
(see: https://www.biorxiv.org/content/[link here])
Various temporary or intermediate files and databases in compressed .csv or .pickle format.
Note: temporary files can be used to resume running from various points.
