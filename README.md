# quantpy
quantpy is a link-wise gaussian and quantics parser, with a focus on Ehrenfest-dynamics simulations

### Features
Built-in flexible link-wise gaussian parser (glogpy)
Quantics input + output parsers

Extraction script drops numpy array binaries - GWP scaled and raw - for plotting as needed
Some simple plotting scripts
- CSF population
- CI State population
- Spin denities
- Muilliken charges
- Bond lengths
- Normal modes

Some other things are in various stages of development
- Dipole moments [NYI]
- QTAIM paramters

 
### Running extractor
Recommended to set up a virtual python enviroment (conda or venv)
Install requirements : 

```
pip install -r requirements.txt
```

Add the folder dir to PYTHONPATH : 

```
export PYTHONPATH=$PYTHONPATH:/path/to/dir
```

Read the analysis script help with 

```
python3 extract.py --h
```

### Visulatations
These scripts are found in analysis folder

- pairplot_sd_mq - Plots up the spin densities or mulliken charges in pairs
- geom - plotting of bond angles/lengths and dihedral angles
- plotter - *new* all in one plotter