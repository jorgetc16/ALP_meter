import Module.StatsFunctions as stat
import Module.FileReading as file
import Module.Plots as plot
import Module.AxionCalculations as calc
import numpy as np


#First we fill the list with the name of the data files (the extension asumed is .dat, if is anjother it should be changed in all the modules)
#and the dark matter density for the earth (first entry) and tçfor each of the sources (in the same order as the filelist)
filelist = ["J00001","J00007","J00008"]

file.ReadAndCompileMultiple(filelist,0) 
#We are searching for peaks so we are interested in running the peak_finder.f file, therefore second argument is 0

file.RunMultiple(filelist,0)
#We are searching for peaks so we are interested in running the peak_finder.f file, therefore second argument is 0

for i in range(len(filelist)):
    plot.plot_periodogram(filelist[i],"Periodogram"+filelist[i],filelist[i])


for i in range(len(filelist)):
    file.AreTherePeaks(filelist[i])