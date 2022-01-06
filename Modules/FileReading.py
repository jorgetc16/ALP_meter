import os
import pandas as pd
from Modules.StatsFunctions import sigray


def ReadAndCompileMultiple(filelist,program):
    """
        This function takes a list of file names as input and compiles the Fortran code
        peak_finder.f or phi95.f on each of them. The output of peak_finder.f are: the Lomb-Scargle
        periodogram (pLS_XXXXX.dat), the FAP values for 32%, 5% and 1% (pLS_FAP_XXXXX.dat)
        and a list of the peaks above 1% (pLS_peaks_XXXXX.dat); the output of phi95.f is
        the 95% CL amplitude for the harmonic signal for each frequency (phi95_XXXXX.dat)
        where XXXXX is the name of each file.

        Parameters
        ----------
        filelist : list
            A list of file names (strings).

        program: integer 
            The number of the program to be run.
                
        Notes
        -----
        The second argument must be 0 for peak_finder.f or 1 for phi_95.f

    """
    if program == 0:
        fortran ="peak_finder"
    elif program == 1:
        fortran ="phi95"
    else:
        print("The second argument must be 0 for peak_finder.f or 1 for phi_95.f")
        return 
        
    df  = pd.read_csv("InputData/"+filelist[0]+".dat",delim_whitespace=True, header=None)
    Size=[]
    Size.append(str(len(df[0])))

    os.system("sed -i 's/npuls1=100/npuls1="+Size[0]+"/g' Modules/"+fortran+".f")
    os.system("sed -i 's/XXXXX.dat/"+filelist[0]+".dat/g' Modules/"+fortran+".f")
    os.system("gfortran -o FortranBinaries/"+filelist[0]+"_"+fortran+".out Modules/"+fortran+".f `cernlib -safe mathlib`")

    for i in range(len(filelist)-1):
        df  = pd.read_csv("InputData/"+filelist[i+1]+".dat",delim_whitespace=True, header=None)
        Size.append(str(len(df[0])))

        os.system("sed -i 's/npuls1="+Size[i]+"/npuls1="+Size[i+1]+"/g' Modules/"+fortran+".f")    
        os.system("sed -i 's/"+filelist[i]+".dat/"+filelist[i+1]+".dat/g' Modules/"+fortran+".f")
        os.system("gfortran -o FortranBinaries/"+filelist[i+1]+"_"+fortran+".out Modules/"+fortran+".f `cernlib -safe mathlib`")
    
    os.system("sed -i 's/"+filelist[len(filelist)]+".dat/XXXXX.dat/g' Modules/"+fortran+".f")
    os.system("sed -i 's/npuls1="+Size[len(filelist)]+"/npuls1=100/g' Modules/"+fortran+".f") 

def ReadAndCompileSingle(filename,program):
    """
        This function takes a file name as input and compiles the Fortran code
        peak_finder.f or phi95.f on it. 

        Parameters
        ----------
        filename : string
            The file name

        program: integer 
            The number of the program to be run.
                
        Notes
        -----
        The second argument must be 0 for peak_finder.f or 1 for phi_95.f

    """
    if program == 0:
        fortran ="peak_finder"
    elif program == 1:
        fortran ="phi95"
    else:
        print("The second argument must be 0 for peak_finder.f or 1 for phi_95.f")
        return 

    df  = pd.read_csv("InputData/"+filename+".dat",delim_whitespace=True, header=None)
    Size=str(len(df[0]))

    os.system("sed -i 's/npuls1=100/npuls1="+Size+"/g' Modules/"+fortran+".f")
    os.system("sed -i 's/XXXXX.dat/"+filename+".dat/g' CompactGithub/"+fortran+".f")
    os.system("gfortran -o FortranBinaries/"+filename+"_"+fortran+".out Modules/"+fortran+".f `cernlib -safe mathlib`")
    os.system("sed -i 's/"+filename+".dat/XXXXX.dat/g' CompactGithub/"+fortran+".f")
    os.system("sed -i 's/npuls1="+Size+"/npuls1=100/g' Modules/"+fortran+".f")

def RunMultiple(filelist,program):
    """
        This function takes a list of file names as input and runs the Fortran code
        peak_finder.f or phi95.f previously compiled. The output of peak_finder.f are: the Lomb-Scargle
        periodogram (pLS_XXXXX.dat), the FAP values for 32%, 5% and 1% (pLS_FAP_XXXXX.dat)
        and a list of the peaks above 1% (pLS_peaks_XXXXX.dat); the output of phi95.f is
        the 95% CL amplitude for the harmonic signal for each frequency (phi95_XXXXX.dat)
        where XXXXX is the name of each file.

        Parameters
        ----------
        filelist : list
            A list of file names (strings).

        program: integer 
            The number of the program to be run.
                
        Notes
        -----
        The second argument must be 0 for peak_finder.f or 1 for phi_95.f

    """
    if program == 0:
        fortran ="peak_finder"
    elif program == 1:
        fortran ="phi95"
    else:
        print("The second argument must be 0 for peak_finder.f or 1 for phi_95.f")
        return 
        
    for i in range(len(filelist)):
        os.system("./FortranBinaries/"+filelist[i]+"_"+fortran+".out")


def RunSingle(filename,program):
    """
        This function takes a file name as input and runs the Fortran code
        peak_finder.f or phi95.f previously compiled. The output of peak_finder.f are: the Lomb-Scargle
        periodogram (pLS_XXXXX.dat), the FAP values for 32%, 5% and 1% (pLS_FAP_XXXXX.dat)
        and a list of the peaks above 1% (pLS_peaks_XXXXX.dat); the output of phi95.f is
        the 95% CL amplitude for the harmonic signal for each frequency (phi95_XXXXX.dat)
        where XXXXX is the name of each file.

        Parameters
        ----------
        filename : string
            The file name

        program: integer 
            The number of the program to be run.
                
        Notes
        -----
        The second argument must be 0 for peak_finder.f or 1 for phi_95.f

    """
    if program == 0:
        fortran ="peak_finder"
    elif program == 1:
        fortran ="phi95"
    else:
        print("The second argument must be 0 for peak_finder.f or 1 for phi_95.f")
        return 

    os.system("./FortranBinaries/"+filename+"_"+fortran+".out")


def GetScale(filelist):
    df_si = pd.DataFrame()
    sk=0
    W=13
    MP=2
    wt=None

    for i in range(len(filelist)):
        df_Pul = pd.read_csv("Output/Files/phi95_"+filelist[i]+".dat",delim_whitespace=True, header=None, skiprows=sk,error_bad_lines=False)
        
        nu_Pul = df_Pul.rolling(W, min_periods=MP,win_type=wt).mean()[0]
        phi_95 = df_Pul.rolling(W, min_periods=MP,win_type=wt).mean()[1]       

        if i == 0:
            nu = pd.DataFrame({"nu": nu_Pul})
            df_si.append(nu)

        Phi = pd.DataFrame({"Pulsar"+filelist[i]: sigray(phi_95)})
        df_si.append(Phi)

    nu = df_si.iloc[0:,0].to_numpy()
    sigma_array = df_si.iloc[0:,1:].to_numpy()
    return nu, sigma_array

def SaveResult(m_a, g_a):

    g_95 = pd.DataFrame({'m_a': m_a,  'g': g_a})
    g_95.to_csv('Output/Files/Result.dat', header=None, index=None, sep ='\t')

def AreTherePeaks(file_name):
    df_Peaks  = pd.read_csv("Output/Files/pLS_peaks_"+source+".dat",delim_whitespace=True, header=None)
    nu = ' '.join(str(e) for e in df_Peaks[1])
    if nu:
        print("In "+file_name+" I have found peaks at "+nu+" days^-1")
    else:
        print("In "+file_name+" I did not find any peak")
