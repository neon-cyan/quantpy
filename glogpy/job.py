import re
import sys
from glogpy.linkparser import linkparsers

# This is a parser for gaussian log files.
# The goal is to offer a linkwise parsing of jobs
# Support of iops in each link is theoretically possible - but, as of yet, not implemented


class link(): # The link class contains both the definition and textual output
    def __init__(self, number, iops, text):
        self.number = number
        self.iops = iops
        self.text = text

class gaussian_job():
    def __init__(self, txt):
        lines = txt.splitlines()
        if 'Normal' not in lines[-1]: raise Exception("[GLEX] Gaussian job ended in error - can not parse")
        
        start = 0
        links = []
        for linenum, line in enumerate(lines):
            # Check if start of link
            if re.search('Enter.{1,}l\d{1,4}.exe', line) != None : 
                #print(line)
                start = linenum

            # Check is end of link
            end = False
            if re.search('Leave Link {1,8}\d{1,3}', line) != None : end = True
            if 'Normal termination of Gaussian' in line : end = True
            if end:
                links.append(lines[start : linenum])
                continue

        # Isolate the route from link1
        x = 0
        for linenum, line in enumerate(links[0]):
            if '-------' in line: x = linenum
        route = links[0][x+1:]
        self.overlays = []

        # Build a list of overlays consisting of [[links], iops , jump]
        for line in route:
            overlay, iops, lks = line.rstrip(';').split('/')
            iopd = {}
            if iops != '':
                iops = iops.split(',')
                for iop in iops:
                    k, v = iop.split('=')
                    iopd[int(k)] = v

            overlay = int(overlay)
            jmp = 0
            if '(' in lks:
                jmp = lks[lks.find('('):].strip('(').rstrip(')')
                jmp = int(jmp)
                lks = lks[:lks.find('(')]
            lks = [int(i) for i in lks.split(',')]
            lks = [overlay*100 + i for i in lks]
            if len(lks) > 1 and jmp != 0:
                raise Exception(f'[GLEX] Found {len(lks)} links for a {jmp} jump in a single line - this is ambigous')
            self.overlays.append([lks, iopd, jmp])
        
        # Fianlly loop through the output links, match them to overlays and add link object to self.link_list
        i, j = 0, 1 # Overlay index // Seen links index [Start at 1 to skip L1]
        
        def extract_link_txt(txt):
                group = re.search('l\d{1,4}.exe', txt[0]).group()
                txt_lk = group.split('.')[0].strip('l')
                return int(txt_lk)
        # print([i[0] for i in links])
        # print(route)
        self.link_list = []
        while True: # Need a looping structure -> accomodate -ve jumps
            for lk in self.overlays[i][0]:
                txt = links[j]
                #print(extract_link_txt(txt), lk)
                assert(extract_link_txt(txt) == lk) # VIBE CHECK - make sure loading link we expect!
                self.link_list.append(link(lk, self.overlays[i][1], txt))
                j += 1
            if lk == 9999: # Last link => Break
                break
            # print(f'jmp={self.overlays[i][2]}')
            if self.overlays[i][2] == 0:         # Spooky shiddy happens here
                i += 1                      
            elif self.overlays[i][2] > 0:        # Positive jumps [skips some steps if possible // Spooky AF]
                if self.overlays[i][2] + i > len(self.overlays)-1:
                    i += 1
                else:
                    i += 1 + self.overlays[i][2]
            else:                           # Negative jumps -> loops
                k = i + self.overlays[i][2]      # Compute the potential jump point

                #print(f'LOOP-DECIDE [{i}] ?-[{extract_link_txt(links[j])}]> L{self.overlays[k][0][0]} @ [{k}]')
                jumplink = extract_link_txt(links[j]) == self.overlays[k][0][0]
                next_link_after_jump = self.overlays[k][0][1] if len(self.overlays[k][0])>1 else self.overlays[k+1][0][0]
                jumplinkp = extract_link_txt(links[j+1]) ==  next_link_after_jump
                if jumplinkp and jumplink: i = k
                else: i += 1

    def prettyprint(self):
        print('ROUTE SECTION:')
        for i, olay in enumerate(self.overlays):
            jmpstring = '🠗' if olay[2] >= 0 else '🠕'
            for j, link in enumerate(olay[0]):
                if len(olay[0]) - 1 == j: print(f'L{link}\t OVERLAY {i} {jmpstring}[{olay[2]}]')
                elif j == 0 : print(f'L{link}\t ∨')
                else : print(f'L{link}\t ∨')

    def parse():
        raise Exception("Attempted to parse a generic job - need to define a inherited class!")
        