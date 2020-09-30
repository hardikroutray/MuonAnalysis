import ROOT
import pandas as pd
import csv
import os
from array import array

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import glob
from matplotlib.colors import LogNorm

def grid(x, y, z, resX=50j, resY=50j):
    from scipy.interpolate import griddata
    xi, yi = np.mgrid[min(x):max(x):resX, min(y):max(y):resY]
    Z = griddata((x, y), z, (xi[None,:], yi[None,:]), method="linear") #, interp="linear")                                                                                          
    #     Z = griddata((x, y), z, (xi[None,:], yi[None,:]), method="cubic", rescale=True) #, interp="linear")                                                                       
    Z = Z[0]
    return xi, yi, Z

x = []
y = []
zexp50 = []
accbr = []
acc = []
acc_lxybin0 = []
acc_lxybin1 = []
acc_lxybin2 = []
acc_lxybin3 = []
acc_lxybin4 = []
acc_lxybin5 = []


# masses = [0.5, 0.6, 0.75, 1.25, 1.5, 2, 2.5, 4, 5, 6, 8, 10, 12, 15, 18, 21, 25]
# masses = [0.5,0.525,0.55,0.575,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.775,0.8,0.825,0.85,0.875,0.9,0.925,0.95,1.25,1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5,3.75,4,4.25,4.5,4.75,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]

masses = [15]

# masses = [2,4]

# ctaus = [15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100]
ctaus = [1]

# ctaus = [1,5,10]

dct = {}
dct1 = {}

for k in range(len(ctaus)):
    dct['ctau{}_exp50'.format(ctaus[k])] = []

for j in range(len(masses)):

    os.chdir("./mass_{}".format(masses[j]))

    dct['mass{}_exp50'.format(masses[j])] = []

    for k in range(len(ctaus)):

        print "looking at mass ", masses[j], "and ctau ", ctaus[k]


        from acc_fit import acceptance
        acceptancearray = acceptance(masses[j],ctaus[k],sample = "hzd")
        # print acceptancearray
        
        acc_lxyall = acceptancearray[0]
        # br_zdtomumu = br
        print "The acceptance of this mass,ctau point is ", acc_lxyall
        
        g = []
        signal_rates = []    
        lxybins = np.array([[0.0,0.2], [0.2,1.0], [1.0,2.4], [2.4,3.1], [3.1,7.0], [7.0,11.0]])

        for l in range(len(lxybins)):

            if  acceptancearray[l+1] > 0.000000000001:
                signal_rates.append(acceptancearray[l+1])
            else:
                signal_rates.append(0.00000001)

        print "The signal rates for lxy bins", signal_rates

        acc_lxybin0.append(signal_rates[0])
        acc_lxybin1.append(signal_rates[1])
        acc_lxybin2.append(signal_rates[2])
        acc_lxybin3.append(signal_rates[3])
        acc_lxybin4.append(signal_rates[4])
        acc_lxybin5.append(signal_rates[5])

        print "The sum of signal rates or total acceptance is", sum(signal_rates)

        if sum(signal_rates) >= 0.01:
            signal_rates = [sr * 10000.0 for sr in signal_rates]
        elif sum(signal_rates) < 0.01 and sum(signal_rates) >= 0.001:
            signal_rates = [sr * 100000.0 for sr in signal_rates]
        elif sum(signal_rates) < 0.001 and sum(signal_rates) >= 0.0001:
            signal_rates = [sr * 1000000.0 for sr in signal_rates]
        elif sum(signal_rates) < 0.0001:
            signal_rates = [sr * 10000000.0 for sr in signal_rates]

        print "The scaled signal rates for lxy bins", signal_rates
        totalsignalrate = sum(signal_rates)

        print "The sum of scaled signal rates", sum(signal_rates)

        for file in glob.glob("simple*.txt"):
            # print(file)
            if "Lxy0.0_0.2_cheb" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy0.0_0.2.txt".format(signal_rates[0],file,masses[j],ctaus[k]))
            elif "Lxy0.0_0.2_counting" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy0.0_0.2.txt".format(signal_rates[0]*0.9244,file,masses[j],ctaus[k]))
            elif "Lxy0.2_1.0_cheb" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy0.2_1.0.txt".format(signal_rates[1],file,masses[j],ctaus[k]))
            elif "Lxy0.2_1.0_counting" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy0.2_1.0.txt".format(signal_rates[1]*0.9244,file,masses[j],ctaus[k]))
            elif "Lxy1.0_2.4_cheb" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy1.0_2.4.txt".format(signal_rates[2],file,masses[j],ctaus[k]))
            elif "Lxy1.0_2.4_counting" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy1.0_2.4.txt".format(signal_rates[2]*0.9244,file,masses[j],ctaus[k]))
            elif "Lxy2.4_3.1_cheb" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy2.4_3.1.txt".format(signal_rates[3],file,masses[j],ctaus[k]))
            elif "Lxy2.4_3.1_counting" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy2.4_3.1.txt".format(signal_rates[3]*0.9244,file,masses[j],ctaus[k]))
            elif "Lxy3.1_7.0_cheb" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy3.1_7.0.txt".format(signal_rates[4],file,masses[j],ctaus[k]))
            elif "Lxy3.1_7.0_counting" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy3.1_7.0.txt".format(signal_rates[4]*0.9244,file,masses[j],ctaus[k]))
            elif "Lxy7.0_11.0_cheb" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy7.0_11.0.txt".format(signal_rates[5],file,masses[j],ctaus[k]))
            elif "Lxy7.0_11.0_counting" in file:
                os.system("sed 's/rate\s100/rate {}/g' {} >  datacard_mass{}_ctau{}_Lxy7.0_11.0.txt".format(signal_rates[5]*0.9244,file,masses[j],ctaus[k]))

        os.system('combineCards.py -S datacard_mass{}_ctau{}_Lxy*.txt > datacard_mass{}_ctau{}_alllxy.txt'.format(masses[j],ctaus[k],masses[j],ctaus[k]))
        os.system('combine -M  AsymptoticLimits -m {} --rAbsAcc=0.0001 --rRelAcc=0.001 datacard_mass{}_ctau{}_alllxy.txt > com.out'.format(masses[j],masses[j],ctaus[k]))
        os.system('cat com.out')                                                                                                                                                    
        os.system('rm datacard_mass{}_ctau{}_Lxy*.txt'.format(masses[j],ctaus[k]))

        com_out = open('com.out','r')                                                                                                                                               
        
        for line in com_out:                                                                                                                                                        
        
            if line[:15] == 'Observed Limit:':                                                                                                                                  
                coml_obs = float(line[19:])                                                                                                                                 

            elif line[:15] == 'Expected  2.5%:':                                                                                                                                
                coml_2sd = float(line[19:])                                                                                                                                 

            elif line[:15] == 'Expected 16.0%:':                                                                                                                                
                coml_1sd = float(line[19:])                                                                                                                                 

            elif line[:15] == 'Expected 50.0%:':                                                                                                                                
                coml_exp = float(line[19:])                                                                                                                                 

            elif line[:15] == 'Expected 84.0%:':                                                                                                                                
                coml_1su = float(line[19:])                                                                                                                                 

            elif line[:15] == 'Expected 97.5%:':                                                                                                                                    
                coml_2su = float(line[19:])
        
        exp_xsec = (coml_exp*totalsignalrate)/(10.1*acc_lxyall)
        print "The expected 50% UL xsec is ", exp_xsec

        x.append(masses[j])
        y.append(ctaus[k])
        zexp50.append(exp_xsec)
        acc.append(acc_lxyall)
        # accbr.append(acc_lxyall*br)

        dct['mass{}_exp50'.format(masses[j])].append(exp_xsec)
        dct['ctau{}_exp50'.format(ctaus[k])].append(exp_xsec)
        


    print "the expected 50% UL xsec for this mass for different ctaus is", dct['mass{}_exp50'.format(masses[j])]
    os.chdir("./..")
 
# print "the expected 50% UL xsec for ctau 20 for different masses is", dct['ctau20_exp50']

print "The ordered mass of all points", x
print "The ordered ctau of all points", y
print "The ordered acceptance of all points",acc
print "The ordered exp 50% UL of all points", zexp50

print "The ordered acceptance of all points in lxybin0",acc_lxybin0
print "The ordered acceptance of all points in lxybin1",acc_lxybin1
print "The ordered acceptance of all points in lxybin2",acc_lxybin2
print "The ordered acceptance of all points in lxybin3",acc_lxybin3
print "The ordered acceptance of all points in lxybin4",acc_lxybin4
print "The ordered acceptance of all points in lxybin5",acc_lxybin5
'''
print "The exp 50% UL of all mass points and ctau 15", dct['ctau15_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau20_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau25_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau30_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau35_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau40_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau45_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau50_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau55_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau60_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau65_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau70_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau75_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau80_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau85_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau90_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau95_exp50']
print "The exp 50% UL of all mass points and ctau 15", dct['ctau100_exp50']
'''

#####
#####
#zobs = df_lims["obs"]                                                                                                                                                               
'''
X, Y, Zexp = grid(x,y,zexp50)
#_, _, Zobs = grid(x,y,zobs)                                                                                                                                                         
fig, ax = plt.subplots()
levels = [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001][::-1]                                                                                               
# levels = [500, 400, 300, 175, 125, 100, 75, 50, 25, 15, 10, 5, 3, 1, 0.5][::-1]
# levels = [,75, 65, 55, 45, 35, 25, 15, 10, 5, 3, 1, 0.1][::-1]
# levels = [20, 18, 15, 10, 5, 1, 0.5, 0.3, 0.1, 0.05, 0.03][::-1]
# levels = np.logspace(np.log10(zexp.min()), np.log10(zexp.max()), 30)                                                                                                               
# levels = 10                                                                                                                                                                        
csexp = ax.contour(X, Y, Zexp, levels=levels, cmap="rainbow")
# csexp = ax.contour(X, Y, Zobs, levels=levels, cmap="rainbow")                                                                                                                      
ax.clabel(csexp, inline=1, fontsize=10, fmt="%2g pb")
ax.scatter(x, y, marker="x", color="k", label="sample")

ax.set_xlabel("mass (GeV)")
ax.set_ylabel(r"c$\tau$ (mm)")

ax.set_xlim([0,25])
# ax.set_ylim([0.,130])                                                                                                                                                              
ax.set_ylim([0.05,130])
# ax.set_xscale("log")                                                                                                                                                               
ax.set_yscale("log")
ax.grid(alpha=0.2)
ax.legend()
ax.text(0.99, 1.01,"10.1 fb${}^\mathregular{-1}$", horizontalalignment='right', verticalalignment='bottom', transform = ax.transAxes, size="large")

fig.savefig("limit_mass_1.png")
'''

hzdlimits = open("hzdlimits_v2.txt", "w")

hzdlimits.write(" masses = {}\n".format(x))
hzdlimits.write(" ctaus = {}\n".format(y))
hzdlimits.write(" xsec = {}\n".format(zexp50))

hzdlimits.write("ctau1_exp50 = {}\n".format(dct['ctau1_exp50']))
hzdlimits.write("ctau5_exp50 = {}\n".format(dct['ctau5_exp50']))
hzdlimits.write("ctau10_exp50 = {}\n".format(dct['ctau10_exp50']))

'''
hzdlimits.write("ctau15_exp50 = {}\n".format(dct['ctau15_exp50']))
hzdlimits.write("ctau20_exp50 = {}\n".format(dct['ctau20_exp50']))
hzdlimits.write("ctau25_exp50 = {}\n".format(dct['ctau25_exp50']))
hzdlimits.write("ctau30_exp50 = {}\n".format(dct['ctau30_exp50']))
hzdlimits.write("ctau35_exp50 = {}\n".format(dct['ctau35_exp50']))
hzdlimits.write("ctau40_exp50 = {}\n".format(dct['ctau40_exp50']))
hzdlimits.write("ctau45_exp50 = {}\n".format(dct['ctau45_exp50']))
hzdlimits.write("ctau50_exp50 = {}\n".format(dct['ctau50_exp50']))
hzdlimits.write("ctau55_exp50 = {}\n".format(dct['ctau55_exp50']))
hzdlimits.write("ctau60_exp50 = {}\n".format(dct['ctau60_exp50']))
hzdlimits.write("ctau65_exp50 = {}\n".format(dct['ctau65_exp50']))
hzdlimits.write("ctau70_exp50 = {}\n".format(dct['ctau70_exp50']))
hzdlimits.write("ctau75_exp50 = {}\n".format(dct['ctau75_exp50']))
hzdlimits.write("ctau80_exp50 = {}\n".format(dct['ctau80_exp50']))
hzdlimits.write("ctau85_exp50 = {}\n".format(dct['ctau85_exp50']))
hzdlimits.write("ctau90_exp50 = {}\n".format(dct['ctau90_exp50']))
hzdlimits.write("ctau95_exp50 = {}\n".format(dct['ctau95_exp50']))
hzdlimits.write("ctau100_exp50 = {}\n".format(dct['ctau100_exp50']))
'''

hzdlimits.close()

'''
matplotlib.pyplot.ioff()

from matplotlib import pyplot as plt1
# plt1.plot(masses,dct['ctau0.5_exp50'],linestyle='--', marker='o', label='ctau = 0.5')
# plt1.plot(masses,dct['ctau1_exp50'],linestyle='--', marker='o', label='ctau = 1')
# plt1.plot(masses,dct['ctau2_exp50'],linestyle='--', marker='o', label='ctau = 2')
# plt1.plot(masses,dct['ctau5_exp50'],linestyle='--', marker='o', label='ctau = 5')
# plt1.plot(masses,dct['ctau1_exp50'],linestyle='--', marker='o', label= 'ctau = 1')
# plt1.plot(masses,dct['ctau4_exp50'],linestyle='--', marker='o', label= 'ctau = 4')
# plt1.plot(masses,dct['ctau8_exp50'],linestyle='--', marker='o', label= 'ctau = 8')
plt1.plot(masses,dct['ctau15_exp50'],linestyle='--', marker='o', label= 'ctau = 15')
plt1.plot(masses,dct['ctau20_exp50'],linestyle='--', marker='o', label= 'ctau = 20')
plt1.plot(masses,dct['ctau30_exp50'],linestyle='--', marker='o', label= 'ctau = 30')
plt1.plot(masses,dct['ctau40_exp50'],linestyle='--', marker='o', label= 'ctau = 40')
plt1.plot(masses,dct['ctau50_exp50'],linestyle='--', marker='o', label= 'ctau = 50')
plt1.plot(masses,dct['ctau60_exp50'],linestyle='--', marker='o', label= 'ctau = 60')
plt1.plot(masses,dct['ctau70_exp50'],linestyle='--', marker='o', label= 'ctau = 70')
plt1.plot(masses,dct['ctau80_exp50'],linestyle='--', marker='o', label= 'ctau = 80')
plt1.plot(masses,dct['ctau90_exp50'],linestyle='--', marker='o', label= 'ctau = 90')
plt1.plot(masses,dct['ctau100_exp50'],linestyle='--', marker='o', label= 'ctau = 100')
# plt1.plot(masses,dct['ctau120_exp50'],linestyle='--', marker='o', label='ctau = 120')
plt1.axvspan(2.9, 3.8, alpha=0.8, color='cyan', label = "J-Psi/Psi'")
plt1.axvspan(0.95, 1.05, alpha=0.5, color='magenta', label = "Phi")
plt1.axvspan(0.75, 0.8, alpha=0.5, color='burlywood', label = "Rho/Omega")

plt1.xlabel('mass[GeV]')
plt1.ylabel('x-sec * BR [fb]')
plt1.yscale("log")
# plt1.xscale("log")
# plt1.show()
plt1.legend(ncol=2,fancybox=True,framealpha=0.2)
plt1.title("Expected 50% UL(H -> ZdZd, 10.1 fb-1)")
plt1.savefig("limits_mass_ctaugt15_v2.png")
'''
