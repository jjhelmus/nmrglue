#!/sbinlab2/software/python-enthought-dis/epd-7.3-2-rh5-x86_64/bin/python2.7
import sys, os
import nmrglue as ng
import matplotlib.pyplot as plt
import matplotlib.cm
import numpy as np
import numpy.lib.recfunctions as rfn
import math

# Integration options
dx = 1
dy = 1
f = open("CPMG.NESSY",'w')

if len(sys.argv) < 2:
    sys.exit('Usage: %s proc_list.txt 1' % sys.argv[0])
# Open CPMG settings
filelist = np.recfromtxt(sys.argv[1], names="ftfile, peakfile, timeT2, NI, nu")

####################
try:
    nrplots = int(sys.argv[2])
except (IndexError,ValueError):
    nrplots = 1
try:
    if int(sys.argv[3]):plot_ppm = True
    else: plot_ppm = False
except (IndexError,ValueError):
    plot_ppm = True
print plot_ppm
####################

# plot parameters
cmap = matplotlib.cm.Blues_r    # contour map (colors to use for contours)
contour_start = 30000           # contour level start value
contour_num = 20                # number of contour levels
contour_factor = 1.20           # scaling factor between contour levels
## calculate contour levels
cl = [contour_start * contour_factor ** x for x in xrange(contour_num)]

## read in the data from a NMRPipe file
dic_list = []
data_list = []
pp_list = []
nu_list = []
for line in filelist:
    ftfile=line['ftfile']
    peakfile=line['peakfile']
    nu=line['nu']
    nu_list.append(nu)
    timeT2=line['timeT2']
    dic, data = ng.pipe.read(ftfile)
    FDF1CAR=dic['FDF1CAR']
    dic_list.append(dic); data_list.append(data)
    peakfiles = ng.fileio.sparkylist.read_HN_N(peakfile)
    peakfilespoints = ng.fileio.pipepoints.calc(dic,peakfiles)
    # Integrate peaks
    vol_list = []
    for peak in peakfilespoints:
        ptsX=peak['ptsX']
        ptsY=peak['ptsY']
        x0 = ptsX-dx
        x1 = ptsX+dx
        y0 = ptsY-dy
        y1 = ptsY+dy
        vol = data[int(round(y0)):int(round(y1)) + 1, int(round(x0)):int(round(x1)) + 1].sum()
        vol_list.append(vol)
    out = np.array(vol_list, dtype=[('Vol',np.float32)])
    out = rfn.merge_arrays((peakfilespoints,out),flatten=True)
    out = out.view(np.recarray)
    pp_list.append(out)
exp_list = zip(dic_list,data_list,pp_list,nu_list) 
# Create nessy database
nr_exp = len(exp_list)
f.write("Project folder:<>nessy_files\n")
reksp = int(round(nr_exp/30))+1
SEQ = ["'%s',"%timeT2 for x in range(reksp)]
f.write("CPMG delay:<>["+' '.join(SEQ)+"]\n")
SEQ = ["'%s'"%FDF1CAR for x in range(reksp)]
f.write("Spec freq:<>["+' ,'.join(SEQ)+"]\n")
f.write("Number of experiments:<>%s\n"%reksp)
max_resi = max(exp_list[0][2]['resi'])
residual_exp = nr_exp
for ds in range(0,reksp):
    residual_exp-=ds*30 
    SEQ = ["''," for x in range(max_resi)]
    for i in range(len(exp_list[0][2])):
        line = exp_list[0][2][i]
        resi = line['resi']
        resn = line['resn']
        SEQ[resi-1]="'%sN-H',"%resn
    f.write("Sequence:<>%s<>["%ds+' '.join(SEQ)+" '']"+"\n")
    if residual_exp > 30: dse = 30
    else: dse = residual_exp
    for i in range(ds*30,ds*30+dse):
        dic, data, peaks,nu = exp_list[i]
        SEQ = ["''," for x in range(max_resi)]
        for j in range(len(peaks)):
            line = peaks[j]
            resi = line['resi']
            vol = line['Vol']
            SEQ[resi-1]="'%s',"%vol
        f.write("Datasets:<>%s<>%s<>["%(ds,i+1-ds*30)+' '.join(SEQ)+" '']"+"\n")
    f.write("CPMG relaxation delay:<>%s<>%s\n"%(ds,timeT2))
    f.write("Experiment:<>%s<>cpmg\n"%ds)
    SEQ = []
    nuN = 0
    for i in range(ds*30,ds*30+dse):
        x = nu_list[i]
        if float(x) == 0:
            x = int(x)
        SEQ.append("'%s',"%x)
        nuN+=1
    for i in range(nuN,29):
        SEQ.append("'',")
    f.write("CPMG frequencies:<>%s<>["%ds+' '.join(SEQ)+" '']"+"\n")
print "Made nessy project."
print "Open with: nessy CPMG.NESSY"

## create the figure
for exp in exp_list[:nrplots]:
    dic, data, peaks,nu = exp
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if plot_ppm: ### make ppm scales
        uc_1_H = ng.pipe.make_uc(dic, data, dim=1)
        ppm_1_H_S, ppm_1_H_E = uc_1_H.ppm_limits()
        uc_15_N = ng.pipe.make_uc(dic, data, dim=0)
        ppm_15_N_S, ppm_15_N_E = uc_15_N.ppm_limits()
    ## plot the contours
    if plot_ppm:
        CS = ax.contour(data, cl, cmap=cmap, extent=(ppm_1_H_S, ppm_1_H_E, ppm_15_N_S, ppm_15_N_E))
    else:
        CS = ax.contour(data, cl, cmap=cmap, extent=(0, data.shape[1]-1, 0, data.shape[0]-1))
    # loop over the peaks
    for peak in peaks:
        ptsX=peak['ptsX']
        ptsY=peak['ptsY']
        ptsX_or=ptsX
        ptsY_or=ptsY
        x0 = ptsX-dx
        x1 = ptsX+dx
        y0 = ptsY-dy
        y1 = ptsY+dy
        resn=peak['resn']
        textp = x1+0.5
        valp = y0-0.5
        vol = data[int(round(y0)):int(round(y1)) + 1, int(round(x0)):int(round(x1)) + 1].sum()
        if plot_ppm:
            ptsX=uc_1_H.ppm(ptsX)
            ptsY=uc_15_N.ppm(ptsY)
            x0=uc_1_H.ppm(x0)
            x1=uc_1_H.ppm(x1)
            y0=uc_15_N.ppm(y0)
            y1=uc_15_N.ppm(y1)
            textp=uc_1_H.ppm(textp)
            valp=uc_15_N.ppm(valp)
        ## plot a box around each peak and label - for points
        ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], 'k')
        ax.plot(ptsX, ptsY, "x", color='y')
        ax.text(textp, y0, resn, size=10, color='r')
        ax.text(x1, valp, "Vol=%1.6e"%vol, size=10, color='r')
        if resn in ['E35']:
            #print x0,x1,y0,y1
            xslice = data[ptsY_or, :]
            yslice = data[:, ptsX_or]
            X = xrange(data.shape[1]) # np.ones(data.shape[1])
            X_V = xslice/contour_start
            Y =  xrange(data.shape[0]) # np.ones(data.shape[0])
            Y_V = -yslice/contour_start
            if plot_ppm:
                #X_H = X_V + ptsY_or
                #Y_H = Y_V + ptsX_or
                #X_T= [(uc_15_N.ppm(x)) for x in X_H]
                #Y_T= [(uc_1_H.ppm(x)) for x in Y_H]
                #ax.plot(X,X_T)
                #ax.plot(Y_T,Y)
                pass
            else:
                X_H = X_V + ptsY_or
                Y_H = Y_V + ptsX_or
                ax.plot(X,X_H)
                ax.plot(Y_H,Y)
    ## decorate the axes
    ax.set_title("Protein H_NH HSQC spectrum")
    plt.clabel(CS, inline=1, fontsize=10)
    labels = ['label-%s'%x for x in range(len(exp_list))]
    if plot_ppm:
        ax.set_ylabel("15N (ppm)")
        ax.set_xlabel("1H (ppm)")
        ax.set_xlim(10.4,6.4)
        ax.set_ylim(132.2, 100.0)
    else:
        ax.set_ylabel("15N (pts)")
        ax.set_xlabel("1H (pts)")
        ax.set_xlim(0, data.shape[1])
        ax.set_ylim(0, data.shape[0])
plt.show()
f.close()
