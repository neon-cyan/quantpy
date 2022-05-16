import numpy as np
import copy

class Stitcher:
    def compute(m, thresh=0.88, quiet=False, doublecheck=True):
        if not quiet:
            print('''
            .dP"Y8 888888 88 888888  dP""b8 88  88 888888 88""Yb     
            `Ybo."   88   88   88   dP   `" 88  88 88__   88__dP     
            o.`Y8b   88   88   88   Yb      888888 88""   88"Yb      
            8bodP'   88   88   88    YboodP 88  88 888888 88  Yb    
            ''')
            print(f'{m.shape[0]} TRAJECTORIES // {m.shape[1]} CI STATES // {m.shape[2]} STEPS // {m.shape[3]} CSF/SD STATES')
            if not doublecheck: print('Will use the best candidate - use with care as you may get active space drift!')

        ans = []

        for ntraj, traj in enumerate(m):
            if not quiet: print(f'Start traj {ntraj + 1}/{m.shape[0]}')
            for nci_state, _ in enumerate(traj):
                if not quiet: print(f'-> Following CI state {nci_state + 1}/{m.shape[1]}')
                i = 0
                current_ci = nci_state
                while i < m.shape[2] -1:
                    # Walk through the whole simulatin step-by-step
                    d1 = np.linalg.norm(m[ntraj][current_ci][i+1] - m[ntraj][current_ci][i])
                    d2 = np.linalg.norm(m[ntraj][current_ci][i+1] + m[ntraj][current_ci][i])
                    if min(d1,d2) > thresh: # States have re-ordered - need a stitch
                        cans = []
                        for j in range(m.shape[1]):
                            d1 = np.linalg.norm(m[ntraj][current_ci][i] - m[ntraj][j][i+1])
                            d2 = np.linalg.norm(m[ntraj][current_ci][i] + m[ntraj][j][i+1])
                            cans.append((j, min(d1,d2)))
                        best, best_value = min(cans, key=lambda x : x[1])
                        if not quiet: print(f'--> Need a stitch @ ({ntraj+1},{i+1}) for states {nci_state+1} [{current_ci+1}] <-> {best+1} [VALUE={best_value:.5f}]')
                        if doublecheck and best_value > thresh : raise Exception("BEST STITCH TOO LARGE {best_value} > {thresh}!")
                        if best_value > thresh and not quiet: print(f'LARGE STITCH {best_value} > {thresh}! DOUBLE CHECK - ACTIVE SPACE DRIFT!')

                        ans.append((ntraj, i, nci_state, best))
                        current_ci = best
                    i += 1

        if not quiet: print(f'STITCHER FINISHED - Suggest {len(ans)} stitches')
        return ans

    def run(m, stitches, quiet=True):
        # M is a tensor with dims traj, ci_state, steps, quantity
        # stitches are a tuple of form traj,step_num,state1,state2
        if not quiet: print(f'''EXECTUING {len(stitches)} STITCHES ON {m.shape}''')
        ans = copy.deepcopy(m)
        for s in stitches:
            if not quiet: print(f"Start on stitch {s}")
            ans[s[0],s[2],s[1]+1:] = m[s[0],s[3],s[1]+1:]
        return ans

class MathUtils:

    def bond_angle(a1, a2, a3, mode="rad"):
        v1 = a1-a2
        v2 = a3-a2
        v1 /= np.linalg.norm(v1)
        v2 /= np.linalg.norm(v2)
        anglerad = np.arccos(np.dot(v1, v2))
        if mode == "rad":
            return anglerad
        elif mode == "deg":
            return np.rad2deg(anglerad)
        else: raise Exception("Illegal mode for bond angle")

    def bond_length(pos1, pos2):
        # Expect a tensor shape 2,3
        d = np.linalg.norm(pos1 - pos2)
        return d

    def dihedral(a_pos, mode="rad"):
        # Expect tensor shape 4,3 A-B-C-D
        # Construct 2 planes containing {AB, BC} and {BC, CD}
        bond_vec1 = a_pos[0] - a_pos[1]
        bond_vec2 = a_pos[3] - a_pos[2]
        common_vec = a_pos[2] - a_pos[1]

        n1 = np.cross(bond_vec1, common_vec)
        n2 = np.cross(bond_vec2, common_vec)

        # Normalize to avoid messy factors
        n1 = n1/np.linalg.norm(n1)
        n2 = n2/np.linalg.norm(n2)

        # Use dot product def. n1.n2=cos(theta)
        anglerad = np.arccos(np.dot(n1, n2))
        if mode == "rad":
            return anglerad
        elif mode == "deg":
            return np.rad2deg(anglerad)
        else: raise Exception("Illegal mode for dihedral angle")


    def gs_ortho(m):
        matrix = copy.deepcopy(m)
        #Diagonalisation
        for i1 in range(matrix.shape[0]):
            for i2 in range(i1):
                x=np.dot(matrix[i1,:],matrix[i2,:])
                matrix[i1,:]=matrix[i1,:]-x*matrix[i2,:]
            x=np.dot(matrix[i1,:],matrix[i1,:])
            x=1/np.sqrt(x)
            matrix[i1,:]=x*matrix[i1,:]
        return matrix

    def dict_to_list(d):
        i = 1
        res = []
        while True:
            if i not in d : break
            res.append(d[i])
            i += 1
        return res

    def moving_avg(x, n):
        cumsum = np.cumsum(np.insert(x, 0, 0)) 
        return (cumsum[n:] - cumsum[:-n]) / float(n)


class NormModeUtils:
    def nm_matrix(AtomMass, freq, vibmatrix):
        #  Some unit coversions
        AtoB = 1/0.529177249  # convert Angstrom to Bohr
        CMtoEV = 0.00012398  # convert from cm-1 to ev
        EVtoAU = 1/27.2116 # convert from ev to Hartree
        AMUtoAU = 1836

        nbmode = np.shape(vibmatrix)[0]
        Natom = len(AtomMass)
        smassvibmatrix = np.zeros_like(vibmatrix)

        #Generate real orthogonal vector by using mass-weighted displacement
        for i in range(nbmode):
            for j in range(Natom):
                smassvibmatrix[i,j,:]=vibmatrix[i,j,:]*np.sqrt(AtomMass[j+1])
            smassvibmatrix[i] /= np.linalg.norm(smassvibmatrix[i,:,:])

        #Reorthogonalize the matrix
        smatrix = np.reshape(smassvibmatrix,(nbmode,Natom*3))
        smatrix = MathUtils.gs_ortho(smatrix)

        nm2xyz = np.zeros((Natom*3, nbmode))
        for j in range(Natom):
            for dim in range(3):
                for i in range(nbmode):
                    idof = j*3+dim
                    div = (AtoB*np.sqrt(AtomMass[j+1]*freq[i]*CMtoEV*AMUtoAU*EVtoAU))
                    nm2xyz[idof, i] = smatrix.T[idof, i] / div
                    smatrix[i, idof] = smatrix[i, idof] * div

        xyz2nm = np.reshape(smatrix, (nbmode, Natom*3))
        return nm2xyz, xyz2nm

    def xyz_to_nm(xyz2nm, geom_init, geomx):
        Natom = geom_init.shape[0]
        nbmode = xyz2nm.shape[0]
        result = np.zeros([geomx.shape[0], geomx.shape[1], nbmode])
        for i in range(geomx.shape[0]): # Loop over GWPs
            for j in range(geomx.shape[1]): # Loop over time steps
                # Compute the XYZ -> NM
                dispxyz=np.subtract(geomx[i,j],geom_init)
                result[i,j] = xyz2nm.dot(np.reshape(dispxyz, Natom*3))
        return result