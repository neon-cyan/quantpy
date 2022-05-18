import numpy as np

# The link parsers are designed to be as simple as possible
# All processing logic ought to be handled elsewhere

# TODO -> Move these into a link class with a .parse() method

class linkparsers():
    def L118(txt):
        # Step one - split the block into TRJ-TRJ-TRJ sections
        trjlines = []
        for i, ln in enumerate(txt):
            if 'TRJ-TRJ-TRJ' in ln:
                trjlines.append(i)
        assert(len(trjlines) % 2 == 0)
        # print(trjlines)
        if len(trjlines) == 2 : trjlines = [trjlines]
        else: trjlines = np.reshape(trjlines, (2, int(len(trjlines)/2)))
        # print(trjlines)
        # Step 2 : loop over the L118 TRJ blocks
        l118_steps = []
        for s, e in trjlines:
            # print(s,e)
            d = {}
            # Define 3 types of L118 block (Input / Startpoint / Iterative)
            if True in map(lambda x : 'INPUT DATA FOR L118' in x, txt[s:e]): d['type'] = 'Input'
            elif True in map(lambda x : 'Start new trajectory calculation' in x, txt[s:e]): d['type'] = 'Start'
            else: d['type'] = 'Iter'

            # Next step is to loop through and pull out the relevant information
            for i, ln in enumerate(txt[s:e]):
                if 'Trajectory Number' in ln:
                    # Parse lines of form
                    # Trajectory Number     1    Step Number     2
                    d['traj'], d['step'] = (int(x) for x in ln.replace('Trajectory Number', '').split('Step Number'))

                elif 'Time (fs)' in ln and 'Trajectory summary' not in txt[s:e][i-1]:
                    # Parse lines of form
                    # Time (fs)     0.139113
                    d['time'] = float(ln.replace('Time (fs)',''))

                elif 'EKin' in ln:
                    # Parse lines of form
                    # EKin =  0.0; EPot =    -1.7; ETot =    -1.7 A.U.
                    d['ekin'], d['epot'], d['etot'] = (float(k.split('=')[1]) for k in ln.rstrip('A.U.').split(';'))

                elif 'Total energy' in ln:
                    # Parse lines of form
                    # Total energy  -1.137304D+02 A.U.
                    # Total energy  -1.137304D+02  Delta-E   0.000000D+00 A.U.                    
                    if 'Delta-E' in ln:
                        te, de = ln.rstrip('A.U.').replace('Total energy','').split('Delta-E')
                        pnum, ppow = (float(j) for j in te.split('D'))
                        d['etot2'] = pnum * pow(10, ppow)
                        pnum, ppow = (float(j) for j in de.split('D'))
                        d['de2'] = pnum * pow(10, ppow)
                    else:
                        pnum, ppow = ln.rstrip('A.U.').replace('Total energy','').split('D')
                        d['etot2'] = float(pnum) * pow(10, int(ppow))


                elif 'Total angular momentum' in ln:
                    # Parse lines of form
                    # Total angular momentum   0.000000D+00 h-bar
                    # Total angular momentum   9.601915D-15  Delta-A  -6.606302D-18 h-bar
                    if 'Delta-A' in ln:
                        te, de = ln.rstrip('h-bar').replace('Total angular momentum','').split('Delta-A')
                        pnum, ppow = (float(j) for j in te.split('D'))
                        d['angmomtot'] = float(pnum) * pow(10, int(ppow))
                        pnum, ppow = (float(j) for j in de.split('D'))
                        d['dangmomtot'] = float(pnum) * pow(10, int(ppow))                        
                    else:
                        pnum, ppow = ln.rstrip('h-bar').replace('Total angular momentum','').split('D')
                        d['angmomtot'] = float(pnum) * pow(10, int(ppow))

                elif 'Cartesian coordinates: (bohr)' in ln:
                    # Parse block of form
                    # Cartesian coordinates: (bohr)
                    # I=    1 X=   2.050651524812D-04 Y=   1.252053407470D+00 Z=   0.000000000000D+00
                    # I=    2 X=   1.792736905953D+00 Y=   2.196753953306D+00 Z=   0.000000000000D+00
                    # I=    3 X=  -1.792014970837D+00 Y=   2.197343547857D+00 Z=   0.000000000000D+00
                    # I=    4 X=  -1.993362381758D-04 Y=  -1.216206700952D+00 Z=   0.000000000000D+00
                    atomdict = {}
                    for atomln in txt[s:e][i+1:]:
                        if 'I=' in atomln:
                            _, label, _, x, _, y, _, z = atomln.split()
                            pnum, ppow = x.split('D')
                            x = float(pnum) * pow(10, int(ppow))
                            
                            pnum, ppow = y.split('D')
                            y = float(pnum) * pow(10, int(ppow))

                            pnum, ppow = z.split('D')
                            z = float(pnum) * pow(10, int(ppow))
                            atomdict[int(label)] = np.array([x,y,z]) * 0.52918 # Convert Bohr -> Angstrom
                        else: break 
                    d['geom'] = atomdict
            l118_steps.append(d)
        return l118_steps

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
        return result # Immitate the CCLIB format [https://cclib.github.io/data.html]
    
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
                    try: atomnum, _, mq, sd = subln.split()
                    except: break
                    muliken[int(atomnum)] = float(mq)
                    spin_den[int(atomnum)] = float(sd)

            elif ('Mulliken charges:' in ln) and not spin_dens:
                # This is the unsummed table
                for subln in txt[i + 2:]:
                    try: atomnum, _, mq = subln.split()
                    except: break
                    muliken[int(atomnum)] = float(mq)

            elif (' Mulliken charges and spin densities with hydrogens summed into heavy atoms:' in ln) and spin_dens:
                # This is the summed H->Heavy Atom table
                for subln in txt[i + 2:]:
                    try: atomnum, _, mq, sd = subln.split()
                    except: break
                    muliken_sum[int(atomnum)] = float(mq)
                    spin_den_sum[int(atomnum)] = float(sd)

            elif ('Mulliken charges with hydrogens summed into heavy atoms:' in ln) and not spin_dens:
                # This is the unsummed table
                for subln in txt[i + 2:]:
                    try: atomnum, _, mq = subln.split()
                    except: break
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

    def L510_TD(txt, do_CI_States=False):
        d = {}
        for i, ln in enumerate(txt):
            if 'ITN=  2 MaxIt=***' in ln:
                en, de = ln.split(' E=')[1].split('DE=')
                de = de.split('Acc=')[0]
                denum, depow = de.split('D')
                d['case'] = float(en)
                d['casde']= float(denum) * (10 ** int(depow))

            elif 'Current Time Dep wavefunction in basis of configuration state functions' in ln:
                diabats = {}
                for csf in txt[i+2:]:
                    try:
                        s, re, im = csf.split()
                        diabats[int(s)] = complex(float(re), float(im))
                    except Exception as e:
                        break
                d['diabats'] = diabats

            elif 'EIGENVALUES AND  EIGENVECTORS OF CI MATRIX' in ln and do_CI_States:
                import re
                cies = {}
                ci_energies = {}
                for nci in txt[i+4:]:
                    if True in map(lambda x  : x in nci, ['iTDHX', '***', 'vector']): break
                    elif 'EIGENVALUE' in nci:
                        # print(nci.split('  '),nci.split('  ')[1],nci.split('  ')[4])
                        state = int(nci.split('  ')[1].replace('(', '').replace(')',''))
                        state_energy = float(nci.split('  ')[4])
                        ci_energies[state] = state_energy
                        cies[state] = {}
                    else:
                        for x in re.findall('\(\ *\d*\)-?\s?\d.\d*', nci):
                            # print(x.replace('(','').split(')'))
                            n,c = x.replace('(','').split(')')
                            cies[state][int(n)] = float(c)

                d['cic'] = cies
                d['cie'] = ci_energies


            elif 'Current Time Dep wavefunction in basis of ci state functions' in ln:
                adiabats = {}
                for cis in txt[i+2:]:
                    try:
                        s, re, im = cis.split()
                        adiabats[int(s)] = complex(float(re), float(im))
                    except Exception as e : 
                        
                        break
                d['adiabats'] = adiabats    

            elif  'Real A(t) energy is' in ln:
                d['raenergy'] = float(ln.replace('Real A(t) energy is', ''))

            elif  'A(t) energy is' in ln:
                d['aenergy'] = float(ln.replace('A(t) energy is', ''))

        # print(d)
        return d