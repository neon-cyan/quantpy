# quantpy
quantpy is a link-wise gaussian and quantics parser, with a focus on Ehrenfest-dynamics simulations using either QUANTICS or L118

### Features
Built-in flexible link-wise gaussian route parser (glogpy)
A few link parsers have been written (Focus on L405, L510)
Rudimentary quantics input + output parsers

### Setup
Recommended to set up a virtual python enviroment (conda or venv)

Install requirements : 

```
pip install -r requirements.txt
```

Add the folder dir to PYTHONPATH : 

```
export PYTHONPATH=$PYTHONPATH:/path/to/dir
```

### Running extractor
Extraction yields a maniest and set of numpy binaries with extracted data

For a quantics job use

```
python3 extract.py /path/to/quantics.inp
```

For a L118 job use

```
python3 extractl118new.py /path/to/l118.log
```

### Analysis tools
These scripts are found in analysis folder

- pairplot_sd_mq - Plots up all the spin densities (or mulliken charges) in pairs
- plotter - All in one plotter. Plot up a range a properties
- gwp_dbg - per-GWP analysis tool to identify problem GWPs & inspect
- to_xyz - Generate a trajectory xyz file for VMD

## TODO
Still missing a few things from old scripts - these are all on the todo list:

- [GWPDBG] Like to refactor bulk per-GWP plotter into single one
- [QEXTRACT] Like to move quantics extractor to use new IO jobs
- [PLOTTER] Plot up dipoles. Hopefully coord system is OK
- [PLOTTER] Plot up Hole Survival Probability
