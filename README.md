# Supplementary Code - README

# Inherent instability of simple DNA repeats shapes an evolutionarily stable distribution of repeat lengths

## Ryan J. McGinty, Daniel J. Balick, Sergei M. Mirkin & Shamil Sunyaev

#### For questions, contact ryan_mcginty@hms.harvard.edu

## For Bioinformatic analysis and computational model:

System Requirements: <br>

Miniconda installation of Python 3.10.14	64bit    - https://docs.anaconda.com/miniconda/ <br>
Additional Python libraries: <br>
Biopython SeqIO		- https://biopython.org/wiki/SeqIO <br>
FastaPy			- https://github.com/aziele/fastapy <br>
Plotly			- https://plotly.com/python/getting-started/ <br>
Liftover		- https://pypi.org/project/liftover/ <br>
Pandas			- https://pandas.pydata.org/ <br>
Numpy			- https://numpy.org/ <br>
Scipy.stats		- https://docs.scipy.org/doc/scipy/reference/stats.html <br>
   <br>
Tested on the following system: <br>
	Windows 11   <br>
	128 Gb RAM   <br>
	AMD 5800X processor   <br>
	Python 3.10.14; pandas 2.2.2; numpy 1.26.4; biopython 1.83;   <br>
	plotly 5.22.0; liftover 1.2.2; scipy 1.13.1   <br>
 <br>
Installation guide: <br>
Follow above instructions to install Python and additional dependencies. <br>
Typical install time: ~15 minutes <br>
 <br>
Instructions for use: <br>
Open Jupyter Notebooks in listed order. <br>
Further instructions are included within the notebook, including instructions for downloading required public datasets and creating directories. <br>
Run each notebook cell in order. <br>
 <br>
Scripts provided for running the full grid of computational models on a Slurm-based cluster. <br>
Run time will vary based on cluster performance. <br>
 <br>
Expected output: <br>
All plots (.pdf images) included in the manuscript. <br>
(see: https://www.biorxiv.org/content/10.1101/2025.01.09.631797v1) <br>
Various temporary or intermediate files and databases in compressed .csv or .pickle format. <br>
Note: temporary files can be used to resume running from various points. <br>
