import numpy as np
import sys, os

try:
    folder = sys.argv[1]
except:
    print(f'use: {sys.argv[0]} /path/to/outfolder')
    sys.exit(-1)
    
assert(os.path.exists(os.path.join(folder, 'times')))
times = np.load(os.path.join(folder, 'times'))
nts = times.shape[0]
print(f'Found a times file with {nts} steps')

if os.path.exists(os.path.join(folder, 'weights')):
    wshape = list(np.load(os.path.join(folder, 'weights')).shape)
    wshape.remove(nts)
    assert(len(wshape) == 1)
    nws = wshape[0]
    print(f'Identified {nws} GWPs')
else: nws = 1

for x in os.listdir(folder):
    # print(x)
    if x == 'times': continue
    try: 
        data = np.load(os.path.join(folder, x))
        assert(nts in data.shape)
    except: continue
            
    print(f'Found a valid {len(data.shape)}D data file {x} - writing ascii')
    tidx = data.shape.index(nts)
    trpv = list(range(len(data.shape)))
    trpv.remove(tidx)

    if nws in data.shape:
        widx = data.shape.index(nws)
        trpv.remove(widx)
        new_tv = [widx, tidx] + trpv
        data = data.transpose(new_tv)
        for i in range(nws):
            with open(os.path.join(folder,f'{x}.{i}.ascii'), 'w') as f:
                for ii, j in zip(times, data[i]):
                    text = '\t'.join([f'{xx: 8.10F}' for xx in j.flatten()])
                    f.write(f'{ii: 8.10F}\t{text}\n')
    else:
        new_tv = [tidx] + trpv
        data = data.transpose(new_tv)
        with open(os.path.join(folder,f'{x}.ascii'), 'w') as f:
            for ii, j in zip(times, data):
                text = '\t'.join([f'{xx: 8.10F}' for xx in j.flatten()])
                f.write(f'{ii: 8.10F}\t{text}\n')
