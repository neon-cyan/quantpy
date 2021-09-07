import numpy as np

# The link parsers are designed to be as simple as possible
# All processing logic ought to be handled elsewhere

class linkparsers():
    def L202(txt):
        sep_ctr = 0
        atoms = {}
        atoms_std = {}
        proton_nums = {}

        for line in txt:
            if '------' in line:
                sep_ctr +=1
                continue
            if sep_ctr==2:
                num, protonnum, atom_type, *coords = line.split()
                atoms[int(num)] = [int(atom_type), np.array(coords).astype(float)]
                proton_nums[int(num)] = int(protonnum)
            if sep_ctr==5:
                num, protonnum, atom_type, *coords = line.split()
                atoms_std[num] = [int(atom_type), np.array(coords).astype(float)]
        return {'geom_symm' : atoms_std, 'proton_nums' : proton_nums, 'geom' : atoms}

    def L716_hpmodes(txt, NAtoms, hpmode = True):
        #  Deal with frequency analysis HPMODE is govereneed by IOP(8)=1_

        #  Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering
        #  activities (A**4/AMU), depolarization ratios for plane and unpolarized
        #  incident light, reduced masses (AMU), force constants (mDyne/A),
        #  and normal coordinates Table
        result = {}
        if hpmode:
            tablestart = None

            for i, line in enumerate(txt):
                if 'and normal coordinates' in line : 
                    tablestart = i + 1
                    break
            if tablestart == None : raise Exception('''Unable to find the frequency table in L716!\n
            Perhaps #P was omitted?''')
            
            setctr = 0

            # Results
            symlabels = []
            freq = []
            redm = []
            forceconst = []
            irint = []
            first_vmat = True

            while True:
                startindex = setctr*(3*NAtoms + 7) + tablestart
                endindex = (setctr + 1)*(3*NAtoms + 7) + tablestart
                if 'Harmonic frequencies' in txt[startindex] : break
                _, symh, freqh, redmh, forceconsth, irinth, _, *atomicinfo = txt[startindex:endindex]

                symh = symh.split()
                symlabels.extend(symh)

                freqh = freqh.split('---')[1].split()
                freq.extend(freqh)

                redmh = redmh.split('---')[1].split()
                redm.extend(redmh)

                forceconsth = forceconsth.split('---')[1].split()
                forceconst.extend(forceconsth)

                irinth = irinth.split('---')[1].split()
                irint.extend(irinth)
                
                atomicinfo = [x.split() for x in atomicinfo]

                matrix = np.array([x[3:] for x in atomicinfo])
                matrix = matrix.astype(float)
                matrix = matrix.transpose()

                new_matrix = []
                for x in matrix:
                    reshaped =  x.reshape((NAtoms, 3))
                    new_matrix.append(reshaped)
                new_matrix = np.array(new_matrix)
                if first_vmat:
                    vibrationmatrix = new_matrix
                    first_vmat = False
                else:
                    vibrationmatrix = np.concatenate((vibrationmatrix, new_matrix))

                setctr += 1
            vibrationmatrix = np.array(vibrationmatrix)
            
            # quick sanity check
            # Make sure have correct number of frequencies
            
            assert(len(symlabels) == vibrationmatrix.shape[0])
            assert(len(freq) == vibrationmatrix.shape[0])
            assert(len(redm) == vibrationmatrix.shape[0])
            assert(len(forceconst) == vibrationmatrix.shape[0])
            assert(len(irint) == vibrationmatrix.shape[0])
            # Double check atoms
            # print(f'atom_numbers={atom_numbers}\atom_proton_num={atom_proton_num}\n')
            result.update({
                'vibdisps'    : vibrationmatrix,
                'vibfreqs'    : np.array(freq).astype(float),
                'vibirs'      : np.array(irint).astype(float),
                'vibsyms'     : symlabels,
                'vibredm'     : np.array(redm).astype(float),
            })
        #  Parse Thermo Table
        tablestart = None
        for i, line in enumerate(txt):
            if 'Thermochemistry' in line : 
                tablestart = i + 2
                break
        if tablestart == None : raise Exception('''Unable to find the thermochemistry table in L716!\n
        Perhaps #P was omitted?''')
        _, temp, _, _, pressure, _ = txt[tablestart].split()
        atom_weights = {}
        atom_proton_num = {}

        setctr = 0
        while True:
            if ' Molecular mass' in txt[tablestart + 1 + setctr] : break
            _, atom_num, _, _, _, atom_p_num, _, _, atom_mass = txt[tablestart + 1 + setctr].split()
            # print(f'an {atom_numbers[setctr]} {atom_num}')
            assert(1 + setctr == int(atom_num)) # make sure that atom numbers are well behaved i.e. 1->2->3
            atom_proton_num[int(atom_num)] = int(atom_p_num) # MAke sure the proton numbers line up
            atom_weights[int(atom_num)] = float(atom_mass)
            setctr+=1

        #  Parse Forces Tables
        tablestart = None
        for i, line in enumerate(txt):
            if 'Forces (Hartrees/Bohr)' in line : 
                tablestart = i
                break
        if tablestart == None : raise Exception('''Unable to find the forces table in L716!\n
        Perhaps #P was omitted?''')
        atom_forces = {}

        setctr = 0
        while True:
            if '---------' in txt[tablestart + 3 + setctr] : break
            atom_num, _, xx, yx, zx = txt[tablestart + 3 + setctr].split()
            atom_f = np.array([xx, yx, zx])
            atom_forces[int(atom_num)] = atom_f.astype(float)
            setctr+=1
        _, maxf, _, rmsf = txt[tablestart + 4 + setctr].split(':')[1].split()

        result.update({
            'atommasses'  : atom_weights,
            'atomnos'     : atom_proton_num,
            'temperature' : float(temp),
            'pressure'    : float(pressure),
            'forces'      : atom_forces,
            'maxforce'    : float(maxf),
            'rmsforce'    : float(rmsf)
        })
        return result # Immitiate the CCLIB format [https://cclib.github.io/data.html]
    
    def L601(txt, spin_dens = False, dipole=False):
        muliken = {}
        muliken_sum = {}

        if spin_dens:
            spin_den = {}
            spin_den_sum = {}

        for i, ln in enumerate(txt):
            if ('Mulliken charges and spin densities:' in ln) and spin_dens:
                # This is the unsummed table
                for subln in txt[i + 2:]:
                    if 'Sum' in subln : break
                    # print(subln)
                    atomnum, _, mq, sd = subln.split()
                    muliken[int(atomnum)] = float(mq)
                    spin_den[int(atomnum)] = float(sd)

            elif ('Mulliken charges:' in ln) and not spin_dens:
                # This is the unsummed table
                for subln in txt[i + 2:]:
                    if 'Sum' in subln : break
                    # print(subln)
                    atomnum, _, mq = subln.split()
                    muliken[int(atomnum)] = float(mq)

            elif (' Mulliken charges and spin densities with hydrogens summed into heavy atoms:' in ln) and spin_dens:
                # This is the summed H->Heavy Atom table
                for subln in txt[i + 2:]:
                    if 'Electronic' in subln : break
                    # print(subln)
                    atomnum, _, mq, sd = subln.split()
                    muliken_sum[int(atomnum)] = float(mq)
                    spin_den_sum[int(atomnum)] = float(sd)

            elif ('Mulliken charges with hydrogens summed into heavy atoms:' in ln) and not spin_dens:
                # This is the unsummed table
                for subln in txt[i + 2:]:
                    if 'APT' in subln : break
                    # print(subln)
                    atomnum, _, mq = subln.split()
                    muliken_sum[int(atomnum)] = float(mq)

            elif ('Dipole moment' in ln) and dipole:
                _, x, _, y, _, z, _, tot = txt[i+1].split()
                dm = [[float(x),float(y),float(z)],float(tot)]
            else: continue

        result = {
            'muliken'            : muliken,
            'mulliken_sum' : muliken_sum,
        }
        if spin_dens:
            result['spinden'] = spin_den
            result['spinden_sum'] = spin_den_sum

        if dipole:
            result['dipole'] = dm
        return result

    def L510_TD(txt):
        diabats = {}    # Called CSFs
        adiabats = {}   # Called CI states
        energy = None
        denergy = None
        for i, ln in enumerate(txt):
            if 'ITN=  2 MaxIt=***' in ln:
                en, de = ln.split(' E=')[1].split('DE=')
                de = de.split('Acc=')[0]
                denum, depow = de.split('D')
                energy = float(en)
                denergy= float(denum) * (10 ** int(depow))
            elif 'Current Time Dep wavefunction in basis of configuration state functions' in ln:
                for csf in txt[i+2:]:
                    try:
                        s, re, im = csf.split()
                        diabats[int(s)] = complex(float(re), float(im))
                    except Exception as e:
                        break


            elif 'Current Time Dep wavefunction in basis of ci state functions' in ln:
                for cis in txt[i+2:]:
                    try:
                        s, re, im = cis.split()
                        adiabats[int(s)] = complex(float(re), float(im))
                    except Exception as e : break
        return {'diabats' : diabats , 'adiabats' : adiabats, 'case' : energy, 'casde' : denergy}