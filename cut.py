import sys, os, json
import numpy as np
 
try:
    _, manifest_dir, cutoff = sys.argv
except:
    print(f'Use f{sys.argv[0]} /path/to/analysis/manifest.json Cutoff[100 | 100fs]')

basepath = os.path.dirname(manifest_dir)
if 'fs' in cutoff:
    ts = float(cutoff.replace('fs', ''))
    times = np.load(os.path.join(basepath, 'times'))
    cutoff=sum(1 if x <= ts else 0 for x in times)
else:
    cutoff=int(cutoff)

with open(manifest_dir, 'r') as f:
    manifest = json.load(f)

print(f'Cutting at {cutoff}')
outputdir=os.path.join(basepath, f'{os.path.dirname(manifest_dir)}_CUT_{cutoff}')
print(f'Writing to {outputdir}')
os.makedirs(outputdir, exist_ok=True)

for f in os.listdir(basepath):
    if '.' in f: continue # Ignore files with extensions
    data = np.load(os.path.join(basepath, f))
    try:
        dimx = data.shape.index(manifest['steps'])
        transposevector = [dimx] + [i if i!=dimx else 0 for i in range(1,len(data.shape))]
        tr = data.transpose(*transposevector)
        assert(tr.shape[0] == manifest['steps'])
        tr = tr[:cutoff]
        data = tr.transpose(*transposevector)
    except:
        print(f'{f} is uncutable! Saving as is')
    with open(os.path.join(outputdir, f), 'wb') as fx:
        np.save(fx, data)
# Fix up and write manifest
manifest['steps'] = cutoff
with open(os.path.join(outputdir, 'manifest.json'),'w') as fx:
    json.dump(manifest, fx)
