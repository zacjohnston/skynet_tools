import numpy as np
import simfiles
import hdf5file
import singlerun
from scipy.interpolate import interp1d

# Get progenitor data
file2 = '/home/zac/projects/data/progenitors/sukhbold_2016/s12.0_presn'
mass, radius, velx, dens, temp, entr, ye = np.loadtxt(file2, usecols=(1, 2, 3, 4, 5, 8, 11), unpack=True, skiprows=2)
# mass, radius, velx, dens, temp, entr, ye = np.loadtxt(file2, usecols=(1, 0, 7, 2, 3, 6, 9),unpack=True,skiprows=11)
print(mass[-1] / 1.989e33)

# Define mass grid
n = 42
start_M = 1.49
end_M = 1.9
mass_grid = np.linspace(start_M, end_M, n)
print (mass_grid[-1])

#  Get FLASH simulation data
path = "/home/zac/projects/data/stir/run_12.0/output/stir2_14may19_s12.0_alpha1.25"

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
table = [None] * (len(defVar))

output = []
out_time = []
out_rad = []
out_mass = []
out_dens = []
out_ye = []
out_temp = []
out_entr = []
out_velx = []
out_enue = []
out_enua = []
out_fnue = []
out_fnua = []

for f in files:
    print (f)
    # Read in data
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
    enue1 = table[defVar.index('enue')]
    enua1 = table[defVar.index('enua')]
    fnue1 = table[defVar.index('fnue')]
    fnua1 = table[defVar.index('fnua')]
    mass2 = []
    rad2 = []
    temp2 = []
    dens2 = []
    ye2 = []
    entr2 = []
    velx2 = []
    enue2 = []
    enua2 = []
    fnue2 = []
    fnua2 = []
    for i in range(len(mass1)):
        mass2.append(mass1[i])
        rad2.append(rad1[i])
        temp2.append(temp1[i])
        dens2.append(dens1[i])
        ye2.append(ye1[i])
        entr2.append(entr1[i])
        velx2.append(velx1[i])
        enue2.append(enue1[i])
        enua2.append(enua1[i])
        fnue2.append(fnue1[i])
        fnua2.append(fnua1[i])

    for i in range(len(mass)):
        if mass[i] / 1.989e33 > mass1[-1]:
            mass2.append(mass[i] / 1.989e33)
            rad2.append(radius[i])
            temp2.append(temp[i])
            dens2.append(dens[i])
            ye2.append(ye[i])
            entr2.append(entr[i])
            velx2.append(velx[i])
            enue2.append(0)
            enua2.append(0)
            fnue2.append(0)
            fnua2.append(0)

    # Map to new mass grid
    f_rad = interp1d(mass2, rad2, kind='next')
    radhold = f_rad(mass_grid)

    f_temp = interp1d(mass2, temp2)
    temphold = f_temp(mass_grid)

    f_dens = interp1d(mass2, dens2)
    denshold = f_dens(mass_grid)

    f_ye = interp1d(mass2, ye2)
    yehold = f_ye(mass_grid)

    f_entr = interp1d(mass2, entr2)
    entrhold = f_entr(mass_grid)

    f_velx = interp1d(mass2, velx2)
    velxhold = f_velx(mass_grid)

    f_enue = interp1d(mass2, enue2)
    enuehold = f_enue(mass_grid)

    f_enua = interp1d(mass2, enua2)
    enuahold = f_enua(mass_grid)

    f_fnue = interp1d(mass2, fnue2)
    fnuehold = f_fnue(mass_grid)

    f_fnua = interp1d(mass2, fnua2)
    fnuahold = f_fnua(mass_grid)

    out_time.append(time1)
    out_rad.append(radhold)
    out_temp.append(temphold)
    out_dens.append(denshold)
    out_ye.append(yehold)
    out_velx.append(velxhold)
    out_enue.append(enuehold)
    out_enua.append(enuahold)
    out_fnue.append(fnuehold)
    out_fnua.append(fnuahold)
tstart = time1

#  Output
for i in range(len(mass_grid)):
    print mass_grid[i]
    with open("../input/traj_s12.0/stir_s12.0_tracer_" + str(i) + ".dat", "w") as f:
        # with open("./data/traj_s11.0/stir2_s11.0_alpha1.25_tracer"+str(i)+".dat","w") as file:
        f.write('mass element = ' + str(mass_grid[i]) + " M_sun\n")
        f.write(
            "   t [s]      T [K]      rho [g/cm^3]      r [cm]      Ye    enue [MeV]   enua [MeV]   fnue [erg/cm^2/s/1e-51]    fnua [erg/cm^2/s/1e-51]\n")
        for j in range(len(out_time)):
            f.write("%18.10E" % (out_time[j]) + "\t" + "%18.10E" % (out_temp[j][i]) + "\t"
                    + "%18.10E" % (out_dens[j][i]) + "\t" + "%18.10E" % (out_rad[j][i]) + "\t"
                    + "%18.10E" % (out_ye[j][i])
                    + "%18.10E" % (out_enue[j][i]) + "%18.10E" % (out_enua[j][i])
                    + "%18.10E" % (out_fnue[j][i]) + "%18.10E" % (out_fnua[j][i])
                    + "\n")
        # f.close()
