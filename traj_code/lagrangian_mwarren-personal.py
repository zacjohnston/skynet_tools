import numpy as np
import scipy
import matplotlib.pyplot as plt
import simfiles
import hdf5file
import singlerun
import mathRecipes
from scipy.interpolate import interp1d

## Get progenitor data
file2= '/Users/mackenziewarren/Research/Woosley_Progenitors/Sukhbold/s12.0_presn'
mass, radius, velx, dens, temp, entr, ye = np.loadtxt(file2, usecols=(1,2,3,4,5,8,11),unpack=True,skiprows=2)
print(mass[-1]/1.989e33)


## Define mass grid
n=1000
mass_grid = [(mass[-1]/1.989e33-2.0)*i/n+2.0 for i in range(n)]


##  Get FLASH simulation data
path = "./Data_12.0/stir2_14may19_s12.0_alpha1.25"

myFiles = simfiles.FlashOutputFiles(path)
print (myFiles.dirname)
myChkfiles = myFiles.chkFilePaths()
files = [chkfiles for chkfiles in myChkfiles]
print ("first checkpoint:", files[0])
print ("last checkpoint:", files[-1])
print ("datfile:", myFiles.datfile)

defVar = hdf5file.FLASH1dFile(files[-1]).loadDefindedVariables()
simTime = (lambda simFile: simFile.simulationTime)
time = singlerun.CCSN1dRun(path).chkfileloop(simTime)
table = [None]*(len(defVar))


output = []
out_time = []
out_rad = []
out_mass = []
out_dens = []
out_ye = []
out_temp = []
out_entr = []
out_velx = []
for f in files:
    ### Read in data
    for ind in defVar:
        data2 = hdf5file.FLASH1dFile(f).getData(str(ind))
        table[defVar.index(ind)] = data2
    time1 = time[files.index(f)]
    temp1 = table[defVar.index('temp')]
    dens1 = table[defVar.index('dens')]
    rad1 = table[defVar.index('radi')]
    ye1 = table[defVar.index('ye  ')]
    entr1 = table[defVar.index('entr')]
    mass1 = table[defVar.index('mass')]
    velx1 = table[defVar.index('velx')]
    mass2 = []
    rad2 = []
    temp2 = []
    dens2 = []
    ye2 = []
    entr2 = []
    velx2 = []
    for i in range(len(mass1)):
        #output.append([mass1[i],rad1[i],dens1[i],temp1[i],entr1[i],ye1[i]])
        mass2.append(mass1[i])
        rad2.append(rad1[i])
        temp2.append(temp1[i])
        dens2.append(dens1[i])
        ye2.append(ye1[i])
        entr2.append(entr1[i])
        velx2.append(velx1[i])
    for i in range(len(mass)):
        if mass[i]/1.989e33 > mass1[-1]:
            mass2.append(mass[i]/1.989e33)
            rad2.append(radius[i])
            temp2.append(temp[i])
            dens2.append(dens[i])
            ye2.append(ye[i])
            entr2.append(entr[i])
            velx2.append(velx[i])

    ## Map to new mass grid
    f_rad = interp1d(mass2,rad2,kind='next')
    radhold = f_rad(mass_grid)

    f_temp = interp1d(mass2,temp2)
    temphold = f_temp(mass_grid)

    f_dens = interp1d(mass2,dens2)
    denshold = f_dens(mass_grid)

    f_ye = interp1d(mass2,ye2)
    yehold = f_ye(mass_grid)

    f_entr = interp1d(mass2,entr2)
    entrhold = f_entr(mass_grid)

    f_velx = interp1d(mass2,velx2)
    velxhold = f_velx(mass_grid)

    out_time.append(time1)
    out_mass.append(masshold)
    out_rad.append(radhold)
    out_temp.append(temphold)
    out_dens.append(denshold)
    out_ye.append(yehold)
    out_entr.append(entrhold)
    out_velx.append(velxhold)
tstart = time1

##  Output
for i in range(len(mass_grid)):
    print mass_grid[i]
    with open("./NDWtraj/stir2_oct8_s12.0_alpha1.25_tracer"+str(i)+".dat","w") as file:
        file.write('mass element = '+str(mass_grid[i])+" M_sun\n")
        file.write("   t [s]      T [K]      rho [g/cm^3]      r [cm]      Ye   Entropy: S =    "+str(out_entr[0][i])+"     k  \n")
        file.write("--------------------------------------------------------------------------------------\n")
        for j in range(len(out_time)):
            file.write("%18.10E" %(out_time[j])+ "\t"+"%18.10E" %(out_temp[j][i])+"\t"+"%18.10E" %(out_dens[j][i])+"\t"+"%18.10E" %(out_rad[j][i])+"\t"+"%18.10E"%(out_ye[j][i])+"\n")
        file.close()
