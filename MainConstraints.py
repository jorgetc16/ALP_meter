import Modules.StatsFunctions as stat
import Modules.FileReading as file
import Modules.Plots as plot
import Modules.AxionCalculations as calc
import numpy as np


#First we fill the list with the name of the data files (the extension asumed is .dat, if is anjother it should be changed in all the modules)
#and the dark matter density for the earth (first entry) and tçfor each of the sources (in the same order as the filelist)
filelist = ["J00001","J00002","J00003","J00004","J00005","J00006","J00007","J00008","J00009","J00010"]
rho = np.array([0.35, 0.35, 0.29, 0.35, 0.33, 0.32, 0.36, 0.53, 0.56, 0.45, 0.44])

file.ReadAndCompileMultiple(filelist,1) 
#We are putting constraints so we are interested in running the phi95.f file, therefore second argument is 1

file.RunMultiple(filelist,1)
#We are putting constraints so we are interested in running the phi95.f file, therefore second argument is 1


nu, sigma = file.GetScale(filelist)
#With the phi95 bounds on a generic harmonic signal we get the scale for the Rayleigh distribution to which we fit the phi


m_a = calc.NuToMa(nu)
n_sources = len(filelist)
n_universes = 100000
n_nu = len(nu)

g_total = np.array(calc.loop(n_sources, n_universes, m_a, sigma, rho, n_nu))

g_limit=[]

for i in range(n_nu):
        g_limit.append(np.percentile(g_total[:,i],95))

file.SaveResult(m_a, g_limit)

y_range = [1.01*10**(-14), 2.99*10**(-8)]
x_range = [.07001, 1100]

plot.plot_bound(m_a,g_limit,x_range, y_range, "Result","New analysis",'forestgreen')