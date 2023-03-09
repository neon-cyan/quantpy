import sys, os
import numpy as np
import matplotlib.pyplot as plt

try: 
    _,  trunc, nm1, nm2, *manifest = sys.argv
except:
    print(f'Use: {sys.argv[0]} TRUNCATION=None|10fs|1233 NM1 NM2 LABEL___/path/to/manifest1.json ...')
    sys.exit(-1)

print (
'''
      ::::    :::   :::   :::                               :::::::: 
     :+:+:   :+:  :+:+: :+:+:   :+:     :+:   :+:     :+: :+:    :+: 
    :+:+:+  +:+ +:+ +:+:+ +:+    +:+ +:+       +:+ +:+         +:+   
   +#+ +:+ +#+ +#+  +:+  +#+ +#++:++#++:++ +#++:++#++:++    +#+      
  +#+  +#+#+# +#+       +#+    +#+ +#+       +#+ +#+     +#+         
 #+#   #+#+# #+#       #+#  #+#     #+#   #+#     #+#  #+#           
###    #### ###       ###                            ##########      
''') 

nm1, nm2 = int(nm1)-1, int(nm2)-1
assert(nm1 >= 0)
assert(nm2 >= 0)

for m in manifest:

    label, path = m.split('___')
    path = os.path.dirname(path)
    nms = np.load(os.path.join(path,'nm_ave'))
    x , y = nms[nm1], nms[nm2]
    if trunc == 'None':
        pass
    elif 'fs' in trunc:
        ts = float(trunc.replace('fs', ''))
        times = np.load(os.path.join(path, 'times'))
        cutoff=sum(1 if x <= ts else 0 for x in times)
        x, y = x[:cutoff], y[:cutoff]
    else: raise Exception('Invalid truncation specified!')


    plt.plot(x,y, label=label)

if len(manifest) > 1: plt.legend(loc='upper right')
plt.title(f'2D Normal mode - cutoff {trunc}' if trunc != 'None' else '2D Normal mode') # Title of the plot
plt.xlabel(f'NM{nm1+1}')
plt.ylabel(f'NM{nm2+1}')
plt.show()
