import numpy as np
import copy

class MathUtils:
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