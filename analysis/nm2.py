import sys, os
import numpy as np
import matplotlib.pyplot as plt

try: 
    _, manifest, nm1, nm2 = sys.argv
except:
    print(f'Use: {sys.argv[0]} /path/to/manifest.json NM1 NM2')
    sys.exit(-1)

print ('''
      ::::    :::   :::   :::                               :::::::: 
     :+:+:   :+:  :+:+: :+:+:   :+:     :+:   :+:     :+: :+:    :+: 
    :+:+:+  +:+ +:+ +:+:+ +:+    +:+ +:+       +:+ +:+         +:+   
   +#+ +:+ +#+ +#+  +:+  +#+ +#++:++#++:++ +#++:++#++:++    +#+      
  +#+  +#+#+# +#+       +#+    +#+ +#+       +#+ +#+     +#+         
 #+#   #+#+# #+#       #+#  #+#     #+#   #+#     #+#  #+#           
###    #### ###       ###                            ##########      
    ''') 

nm1, nm2 = int(nm1)-1, int(nm2)-1
path = os.path.dirname(manifest)
nms = np.load(os.path.join(path,'nm_ave'))

plt.plot(nms[nm1], nms[nm2])
plt.title('2D Normal mode') # Title of the plot
plt.xlabel(f'NM{nm1+1}')
plt.ylabel(f'NM{nm2+1}')
plt.show()
