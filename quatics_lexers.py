import re

class QuanticsParsers():
    def parse_input(data):
        # Wil return a dict containing important quantics inputs:
        # name = QUANTICS output dir
        # ngwp = Number of GWPs
        # data = GAUSSIAN output dir
        # freqf = Name of the frequency file
        # tfinal = simulation length (fs)
        # tpsi = step size

        d = {}
        for line in data.splitlines():
            try: 
                if line[0] == '#': continue
            except IndexError: continue # Skip empty lines

            if 'name' in line:
                group = re.search('name\s?=\s?\S{1,}', line).group()
                d['name'] = group.split('=')[1].strip() # Clean up whitespace

            if 'ngwp' in line:
                group = re.search('ngwp\s?=\s?\d{1,}', line).group()
                d['ngwp'] = int(group.split('=')[1])
            
            if ('data' in line) and ('subcmd' not in line):
                group = re.search('data\s?=\s?\S{1,}', line).group()
                d['data'] = group.split('=')[1].strip() # Clean up whitespace

            if 'transfile' in line:
                group = re.search('transfile\s?=\s?\S{1,}', line).group()
                d['freqf'] = group.split('=')[1].strip() # Clean up whitespace

            if 'tfinal' in line:
                group = re.search('tfinal\s?=\s?\d{1,}', line).group()
                d['tfinal'] = float(group.split('=')[1])

            if 'tpsi' in line:
                group = re.search('tpsi\s?=\s?(\d|\.){1,}', line).group()
                d['tpsi'] = float(group.split('=')[1])

        return d
    def parse_output(data):
        # Wil return a dict containing quantics output lists:
        # time  = sim timestamp
        # Etot  = total energy
        # DelE  = Delta E
        # DiagD = Diagonal Densities * 1000
        # GGP   = Gross Gaussian Populations * 10
        # 
        # MEQ   = expectation []
        # MEdQ  = Variance in expectation []
        # MEP   = momentum []
        # MEdP  = Variance in momentum []

        def parse_dd_ggp(inp):
            ans = []
            for line in inp.splitlines()[1:]:
                if ':' in line:
                    nums = line.split(':')[1]
                elif '>' in line:
                    nums = line.split('>')[1]
                nums = [float(i) for i in nums.split()]
                ans.extend(nums)
            return ans

        DIVISOR = '{}\n\n'.format('-'*78)
        
        datapoints = []

        for g in data.split(DIVISOR)[1:-1]:
            d = {}

            group = re.search('Time\ {1,}=\ {1,}-?(\.|\d){1,}', g).group()
            d['time'] = float(group.split('=')[1]) # Clean up whitespace

            group = re.search('E-tot\ {1,}=\ {1,}-?(\.|\d){1,}', g).group()
            d['Etot'] = float(group.split('=')[1])

            group = re.search('Delta-E\ {1,}=\ {1,}-?(\.|\d){1,}', g).group()
            d['DelE'] = float(group.split('=')[1])

            _, diagden, grosspop, modeexp = g.split('\n\n')

            d['DiagD'] = parse_dd_ggp(diagden)
            d['GGP'] = parse_dd_ggp(grosspop)

            meq_arr = []
            medq_arr = []
            mep_arr = []
            medp_arr = []

            for line in modeexp.splitlines()[1:]:
                nums = line.split(':')[1].strip().split()[1::2] 
                meq, medq, mep, medp = [float(i) for i in nums]
                meq_arr.append(meq)
                medq_arr.append(medq)
                mep_arr.append(mep)
                medp_arr.append(medp)
            d['MEQ'] = meq_arr
            d['MEdQ'] = medq_arr
            d['MEP'] = mep_arr
            d['MEdP'] = medp_arr
            datapoints.append(d)
        return datapoints