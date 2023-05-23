import sys
import numpy as np

nbf = int(sys.argv[1])

ans = np.zeros(nbf)

while True:
    i = input()
    if i == '': break
    try:
        p, v = i.split()
        ans[int(p)-1] = float(v)
    except:
        print('Invalid input')
    
print(','.join((str(x) for x in ans)))
