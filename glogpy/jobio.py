import re
import sys
from glogpy.linkparser import linkparsers

# This is a parser for gaussian log files.
# The goal is to offer a linkwise parsing of jobs
# This version is IO heavy but light on RAM.
# If parsing small logs (<1gb or so) please use job class instead that laods file into RAM
# Support of iops in each link is theoretically possible - but, as of yet, not implemented


class linkio(): # The link class contains both the definition and textual output
    def __init__(self, number, iops, start, end):
        self.number = number
        self.iops = iops
        self.start = start
        self.end = end

    def text_parse(self, txtfile, fn, *a, **ka):
        with open(txtfile, 'r') as f:
            f.seek(self.start)
            data = f.read(self.end-self.start)
        data = data.splitlines()
        r = fn(data, *a, **ka)
        del(data) # Clean up
        return r

class gaussian_jobio():
    def __init__(self, txtfile, pos=0, allow_partial=False):
        self.txtfile = txtfile
        self.pos = pos
        link1=False

        with open(txtfile,'r') as f:
            f.seek(pos)
            ln = f.readline()
            while ln:
                # print(ln)
                if 'Entering Link 1' in ln: # Special case for Link1 (route)
                    self.route = []
                    # If you fed in a multi-job log we reset to make sure last job is properly parsed
                    self.links= [] 
                    start = f.tell()
                    link1=True

                elif re.search('Enter.{1,}l\d{1,}.exe', ln) != None : 
                    start = f.tell()

                elif re.search('Leave Link {1,}\d{1,}', ln) != None : 
                    s = re.search('Leave Link {1,}\d{1,}', ln)[0]
                    lk = int(s.split()[-1])
                    end = f.tell()
                    link1=False
                    self.links.append([lk, start,end])

                elif link1: # ROUTE BUILDER
                    if '-----------------' in ln:
                        self.route = [] # Clear the route in case of a NonStd conflict
                        
                    s = re.search('\d+\/(\d+\=-?\d+,?)*\/(\d,?)*(\(-?\d+\))?', ln)
                    if s != None:
                        # Split lines like 1/18=10,22=1/23(-9);
                        i, j, k = s[0].split('/')

                        overlay = int(i)
                        
                        if j.split(',') == ['']: iops = {}
                        else: iops={int(k):int(v) for k,v in [x.split('=') for x in j.split(',')]}
                        
                        k = k.split('(')
                        lks = [int(x) for x in k[0].split(',')]
                        jmp = 0 if len(k) == 1 else int(k[1].replace(')','')) # #Deals with 1,2,3(1)

                        if len(lks) > 1 and jmp != 0:
                            raise Exception(f'[GLEX] Found {len(lks)} links for a {jmp} jump in a single line - this is ambigous!')

                        self.route.append([[x+overlay*100 for x in lks], iops, jmp])

                elif 'Normal termination of Gaussian' in ln: # Special case for L9999
                    end = f.tell()
                    lk = int(9999)
                    self.links.append([lk, start,end])

                ln = f.readline()

        # print(self.links, self.route)

        # Finally loop through the output links, match them to overlays
        i, j = 0, 1
        self.link_list = []
        while True: # Need a looping structure -> accomodate -ve jumps
        
            for lk in self.route[i][0]:
                # print(self.links[j][0] , lk)
                # print( allow_partial, j, len(self.links))
                if allow_partial and len(self.links) == j: 
                    self.partial=True
                    return
                assert(self.links[j][0] == lk) # VIBE CHECK - make sure loading link we expect!
                self.link_list.append(linkio(lk, self.route[i][1], self.links[j][1], self.links[j][2]))
                j += 1
                
            if lk == 9999: # Last link => Break
                self.partial=False
                return

            if self.route[i][2] == 0:
                i += 1

            elif self.route[i][2] > 0: # Positive jumps [skips some steps if possible // Spooky AF]
                if self.route[i][2] + i > len(self.route)-1:
                    i += 1
                else:
                    i += 1 + self.route[i][2]
            else:                           
                # Negative jumps -> loops
                k = i + self.route[i][2] # Compute the potential jump point
                jumplink = self.links[j][0] == self.route[k][0][0]
                next_link_after_jump = self.route[k][0][1] if len(self.route[k][0])>1 else self.route[k+1][0][0]
                jumplinkp = self.links[j+1][0] == next_link_after_jump
                if jumplinkp and jumplink: i = k
                else: i += 1

    def parse():
        raise Exception("Attempted to parse a generic job - need to define a inherited class!")
        