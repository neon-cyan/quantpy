# quantpy
quantpy is a link-wise gaussian and quantics parser, with a focus on Ehrenfest-dynamics simulations

### Features
Built-in flexible link-wise gaussian parser (glogpy) - may wish to factor this out as a submodule in future
Rudimentary quantics input + output parsers

Extraction script drops numpy array binaries - GWP scaled and raw - for manipulating & plotting as needed

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

- pairplot_sd_mq - Plots up all the spin densities (or mulliken charges) in pairs
- plotter - *new* all in one plotter. Plot up a range a properties

## TODO
Missing a few things from old scripts - these are all on the todo list:
- [ANALYSIS] Single GWP analysis
- [PLOTTER] Plot up dipoles. Hopefully coord system is OK
- [ANALYSIS] New debugging plot suite? CAS DE / MaxForce / CAS Energy (with hints?)