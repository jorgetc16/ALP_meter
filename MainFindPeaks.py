import Modules.StatsFunctions as stat
import Modules.FileReading as file
import Modules.Plots as plot
import Modules.AxionCalculations as calc
import numpy as np


#First we fill the list with the name of the data files (the extension asumed is .dat, if is anjother it should be changed in all the modules)
#and the dark matter density for the earth (first entry) and tçfor each of the sources (in the same order as the filelist)
filelist = ["J00001","J00002","J00003","J00004","J00005","J00006","J00007","J00008","J00009","J00010"]

file.ReadAndCompileMultiple(filelist,0) 
#We are searching for peaks so we are interested in running the peak_finder.f file, therefore second argument is 0

file.RunMultiple(filelist,0)
#We are searching for peaks so we are interested in running the peak_finder.f file, therefore second argument is 0

for i in range(len(filelist)):
    plot.plot_periodogram(filelist[i],"periodogram"+filelist[i],filelist[i])


for i in range(len(filelist)):
    file.AreTherePeaks(filelist[i])