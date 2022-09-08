from glogpy.freqency_job import frequency_job
import sys
import numpy as np
import copy
import os
import mathutils

try:
    _, ffile, nm, min_nm, max_nm, steps, *template = sys.argv

except:
    print(f'Use: python {sys.argv[0]} freqjob.log NM START STOP STEP (optional template)')
    sys.exit(-1)

try:
    with open(ffile, 'r') as f:
        data = f.read()
    data = data.split('Initial command:\n')[-1] # Assume the freq job is the last
    job = frequency_job(data)
    # job.prettyprint()
    data = job.parse()
except:
    print('Invalid frequency file')
    sys.exit(-1)

if len(template) == 0:
    template = None
else:
    try:
        outpath = os.path.dirname(template[0])
        with open(template[0], 'r') as f:
            template = f.read()
    except:
        raise Exception('Inalid template!')

nm = int(nm)-1 # Zero indexing!
min_nm, max_nm = float(min_nm), float(max_nm)
steps = int(steps)

ATOMICLABELS = ['H',  'He',  'Li',  'Be',  'B',  'C',  'N',  'O',  'F',  'Ne',  'Na',  'Mg',  'Al',  'Si',  'P',  'S', 
    'Cl',  'Ar',  'K',  'Ca',  'Sc',  'Ti',  'V',  'Cr',  'Mn',  'Fe',  'Co',  'Ni',  'Cu',  'Zn',  'Ga',  'Ge',  'As',  'Se',  
    'Br',  'Kr',  'Rb',  'Sr',  'Y',  'Zr',  'Nb',  'Mo',  'Tc',  'Ru',  'Rh',  'Pd',  'Ag',  'Cd',  'In',  'Sn',  'Sb',  'Te', 
    'I',  'Xe',  'Cs',  'Ba',  'La',  'Ce',  'Pr',  'Nd',  'Pm',  'Sm',  'Eu',  'Gd',  'Tb',  'Dy',  'Ho',  'Er',  'Tm',  'Yb', 
    'Lu',  'Hf',  'Ta',  'W',  'Re',  'Os',  'Ir',  'Pt',  'Au',  'Hg',  'Tl',  'Pb',  'Bi',  'Po',  'At',  'Rn',  'Fr',  'Ra', 
    'Ac',  'Th',  'Pa',  'U',  'Np',  'Pu',  'Am',  'Cm',  'Bk',  'Cf',  'Es',  'Fm',  'Md',  'No',  'Lr',  'Rf',  'Db',  'Sg', 
    'Bh',  'Hs',  'Mt',  'Ds',  'Rg',  'Cn',  'Nh',  'Fl',  'Mc',  'Lv',  'Ts',  'Og'] # Yes *all* of the elements

atomnumbers = sorted(data['geom'].keys())
geom = np.array([data['geom'][i][1] for i in atomnumbers])
elems = np.array([ATOMICLABELS[data['proton_nums'][i]-1] for i in atomnumbers])
# nm_oi = data['vibdisps'][nm]

# Cosntrict NM Matrix & maniopulate to get to needed shape
nm2xyz, xyz2nm = mathutils.NormModeUtils.nm_matrix(data['atommasses'],  data['vibfreqs'], data['vibdisps'])
newnm = nm2xyz.T[nm]
newnm = np.reshape(newnm, (-1, 3))

# print(geom)
# print(elems)
# print(nm_oi)

print(f'{min_nm} --> {steps} steps --> {max_nm}')
# print(ffile, nm, min_nm, max_nm, steps)

n = 0
for x in np.linspace(min_nm,max_nm,steps):
    xyz = geom + x * newnm
    g = '\n'.join([f'{i[0]}   {i[1][0]:10.7f}   {i[1][1]:10.7f}   {i[1][2]:10.7f}' for i in zip(elems, xyz)])
    if template == None:
        print(f'{x:.5f}')
        print(g)
    else:
        text = copy.deepcopy(template)
        text = text.replace('#XYZ#', g)
        text = text.replace('#DISP#', f'{x:.5f}')
        out = os.path.join(outpath, f'{n}.com')
        print(out)
        with open(out, 'w') as f:
            f.write(text)
        n+=1