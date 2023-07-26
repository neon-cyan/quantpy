from glogpy.jobio import gaussian_jobio
import numpy as np
from glogpy.linkparser import linkparsers


# This is a definition for a single traj Eh L118 job
# TODO test with GradOnly fixed timestep

class l118_job(gaussian_jobio):
    def __init__(self, txtf, allow_partial=False):
        super().__init__(txtf, allow_partial=allow_partial)

    def parse(self, print_info=True, allow_truncate=False):

        # Check the links in l118 loop - L405 and 510 MUST be there, 601 MAY be there
        loopidx, olay = list(filter(lambda x : 118 in x[1][0] and x[1][2]<0, enumerate(self.route)))[0]
        loop_n = olay[2]
        # print(loopidx, loop_n)
        linkinloop = set()
        for x in self.route[loopidx+loop_n: loopidx]:
            linkinloop = linkinloop.union(x[0])
        # print(linkinloop)
        assert(510 in linkinloop)
        assert(405 in linkinloop)
        do_spin_ana = True if 601 in linkinloop else False # Do not do spin if not in L118 loop
        do_force_ana = True if 716 in linkinloop else False # Do not do force if not in L118 loop
        if print_info: print('Found suitable L118 loop : L510 OK L405 OK L601 {} L716 {}'.
            format('OK' if do_spin_ana else 'ABSENT', 'OK' if do_force_ana else 'ABSENT'))

        # Next, check Number of L510s & 405s
        l405s = list(filter(lambda x: x.number==405, self.link_list))
        l510s = list(filter(lambda x: x.number==510, self.link_list))
        if not allow_truncate: assert(len(l510s) == len(l405s)) # MUST have a L4105 per L510
        assert(len(l510s) > 0) # MUST have at least one step!
        if print_info: print(f'Found {len(l510s)} x L510s. Consistent with L405==L510')

        # Next check L202s - need just the first one for init_geom & atomic labels
        l202s = list(filter(lambda x: x.number==202, self.link_list))
        assert(len(l202s) >= 1)
        l202_init = l202s[0].text_parse(self.txtfile, linkparsers.L202)
        if print_info: print(f'Found a suitable L202')

        l405 = l405s[0].text_parse(self.txtfile, linkparsers.L405)     
        if print_info: print(f'Parsed the first L405 OK')

        # Need to decide if spindens is available
        spin_dens = False not in [72 in x.iops for x in l510s]
        if print_info: print('{}'.format('Spin denisty will be extracted' if spin_dens else '72=1 not set in L510 => No Spin Density'))

        # Check how many steps to parse
        l118s = list(filter(lambda x: x.number==118, self.link_list))

        # Need a special case for L118 IOP(5) => Initial momenta NM specified
        # L118 provides 2 steps in first call otherwise...
        double_118 = False if 5 in l118s[0].iops else True
        # Extra case if L118 pulse logic is used
        if 9 in l118s[0].iops:
            if l118s[0].iops[9] in [10,11] : double_118 = False

        if print_info: 
            print(f'Found {len(l118s)} x L118. Will discard the first step with integration info ...')
        l118s = l118s[1:] # Drop the first L118 datapoint with integration info

        tot = len(l118s)
        # if double_118 : tot +=1 
        # When double L118 it will double print in the beginning but have an extra point at the end

        # Generate the 610s list if needed
        if do_spin_ana: l601s = list(filter(lambda x: x.number==601, self.link_list))
        if do_force_ana: l716s = list(filter(lambda x: x.number==716, self.link_list))

        if allow_truncate:
            # Use for incomplete jobs
            if print_info: print('Checking truncation ...')
            totnew = min(tot, len(l510s))
            if do_spin_ana: totnew = min(totnew, len(l601s))
            if do_force_ana: totnew = min(totnew, len(l716s))
            if tot != totnew:
                if print_info:
                    print('About to truncate from {tot} -> {totnew} steps')
                    print('MAKE SURE THIS BEHAVIOUR IS EXPECTED!')
                tot = totnew

        if print_info: print(f'Expect to parse out {tot} steps')

        printmod = int(np.ceil(tot/10))
        xn = []
        print(double_118)
        for n in range(tot):
            if print_info and (n+1)%printmod == 0: print(f'=> Parsing {n+1}/{tot}')
            
            if double_118:
                # print(l118s[0].text_parse(self.txtfile, linkparsers.L118))
                if n == 0: info118 = l118s[0].text_parse(self.txtfile, linkparsers.L118)[0]
                elif n == 1: info118 = l118s[0].text_parse(self.txtfile, linkparsers.L118)[1]
                else: info118 = l118s[n-1].text_parse(self.txtfile, linkparsers.L118)[0]
            else: 
                info118 = l118s[n].text_parse(self.txtfile, linkparsers.L118)[0]

            info510 = l510s[n].text_parse(self.txtfile, linkparsers.L510_TD, do_CI_States=True)
            # print(info510)

            info601 = l601s[n].text_parse(self.txtfile, linkparsers.L601, spin_dens=spin_dens) if do_spin_ana else {}
            info716 = l716s[n].text_parse(self.txtfile, linkparsers.L716_hpmodes, NAtoms=len(l202_init['proton_nums']), hpmode=False) if do_force_ana else {}
            xn.append([info510, info118, info601, info716])
            n += 1
        if print_info: print('Parser finished OK')
        return (l202_init, xn, l405)