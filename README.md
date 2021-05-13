# quantpy
quantpy is a link-wise gaussian and quantics parser

### Features
Extract and plot up data scaled by GWP populations
- Normal mode evolution
- Bond lengths
- CI / CSF State population evolution
- Spin densities (with H summed into heavy atoms)
- Mulliken charges (with H summed into heavy atoms)
 
### Running
Recommended to set up a virtual python enviroment (conda or venv)
Install requirements : 

```
pip install -r requirements.txt
```
Add the folder dir to PYTHONPATH : 

```
export PYTHONPATH=$PYTHONPATH:/path/to/dir
```
Start the analysis script with 

```
python3 quantpy.py [args ...]
```