"""
sparkylist read sparky list files 

The header of the file can be different according to export from 
SPARKY or Analysis.

The header can exist, not exist, and have an empty line space-between header and peak list.
Also the column order of HN, and N, is different.

Analysis export: Resn - HN - N - no line skip between header and files. Have Height/Volumne column.
SPARKY export: Resn - N - HN - one line skip between header and files
Personal peak: Probably don't have header line.

------------------------------
Ex 1: Analysis export.
------------------------------
      Assignment         w1         w2       Height       Volume
              ?-?      9.061    132.333   1.14E+05   7.66E+05 bx
         V2HN-V2N      8.856    131.178   2.33E+03   9.26E+03 bx
         F3HN-F3N      8.856    131.178   2.33E+03   9.26E+03 bx
------------------------------
Ex 2: Sparky export
------------------------------
      Assignment         w1         w2

           G4N-HN    108.175      8.386
           E7N-HN    129.130      8.131
------------------------------

"""
import numpy as np
import csv
import re
import numpy.lib.recfunctions as rfn

one2three ={'V':'VAL', 'I':'ILE', 'L':'LEU', 'E':'GLU', 'Q':'GLN', 'D':'ASP', 'N':'ASN', 'H':'HIS', 'W':'TRP', \
'F':'PHE', 'Y':'TYR', 'R':'ARG', 'K':'LYS', 'S':'SER', 'T':'THR', 'M':'MET', 'A':'ALA', \
'G':'GLY', 'P':'PRO', 'C':'CYS'}

def read_HN_N(peakfile):
    # Test the header order of the file
    with open(peakfile, 'r') as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
        lines = []
        for line in reader:
            lines.append(line)
    header = lines[:3]
    emptyline = -1 
    hasheader = -1
    haspeak = -1
    peakorder = ['resn','N','HN']
    headerorder = peakorder
    for i in range(len(header)):
        line = header[i]
        # Test emptylines
        try:
            empty = line[1]
        except IndexError:
            emptyline = i
            continue
        # Test for header/peak
        try: #If peak
            isfloat1 = float(line[1])
            isfloat2 = float(line[2])
            haspeak = i
            if isfloat1 > 50.0 and isfloat1 < 50.0:
                peakorder = ['resn','N','HN']
            elif isfloat1 < 50.0 and isfloat1 > 50.0:
                peakorder = ['resn','HN','N']
        except ValueError: #If header
            hasheader = i
            headerorder = header[i]
    headerorder[:3] = peakorder
    skiplines = max(hasheader,emptyline)+1
    # Read the peak file as structured array
    pf = np.recfromtxt(peakfile, names=headerorder, skip_header=skiplines)
    regint = re.compile('\d+')
    resn_one_list = []
    resn_three_list = []
    resi_list = []
    for i in range(len(pf['resn'])):
        resnstr = pf['resn'][i]
        resi = regint.findall(resnstr)
        resi = resi[0]
        resi_list.append(resi)
        resimatch = re.search(resi,resnstr)
        resis,resie = resimatch.span()
        resn_one = resnstr[:resis]
        resn_one_list.append(resn_one)
        resn_three = one2three[resn_one]
        resn_three_list.append(resn_three)
        pf['resn'][i] = resnstr[:resie]
    out = zip(resn_one_list,resn_three_list,resi_list)
    out = np.array(out, dtype=[('oneAA','|S1'),('AA','|S3'),('resi',int)])
    #out = out.view(np.recarray)
    out = rfn.merge_arrays((pf,out),flatten=True)
    out = out.view(np.recarray)
    return(out)

           


