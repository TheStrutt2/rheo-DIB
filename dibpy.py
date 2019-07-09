## Import plots, and numpy
import matplotlib.pyplot as plt
import numpy as np
import scipy
import glob, os
from scipy.optimize import curve_fit
from scipy.integrate import odeint
import pandas as pd

# Represents Ficks first law as interpreted in terms of RheoDib permeability
def PermeabilityRheoDIB (VA,A,CEQ,CAt,dCdt):
    "Returns a value for P"
    return (VA / (2*A*(CEQ-CAt))) * (dCdt)

#Determines the area of a circle 
def circlearea (r):
    "Return the area of a circle based on the radius"
    return np.pi * r **2

def Propagation_of_P_error (P,dVA,VA,dA,A,dcdtse,diffterm,dC,CEQ,CAt):
    "Returns a value for the propagation of error of P"
    return P*(sqrt(((dVA/VA)**2) + ((dA/A)**2) + ((dcdtse/diffterm)**2)+((dC/(CEQ-CAt))**2)))


### File: previously OpenAndFitData1.py
## Extract data from file
def extractSpectraFromFile(file_name):
    # open file, get data, close file
    with open(file_name) as myFile:
        text = myFile.read()
        delims_str = '>>>>>Begin Spectral Data<<<<<'
        sample_info, spectra = text.split(delims_str) # split data into sample info and raw spectra
        myFile.close()
    return spectra, sample_info

## Convert spectra from str to float
def convertSpectraToFloat(spectr_str):
    # split string into list
    spectra = spectr_str.split()
    
    wavelength = []
    intensity = []

    for i, val in enumerate(spectra):
        if i % 2 == 0:
            wavelength.append(float(val))
        else:
            intensity.append(float(val))
    
    return wavelength, intensity


## Normalzie and reduce data
def normalizeAndReduce(wavelength, intensity, lower_bin, upper_bin):
    count1 = 0
    count2 = 0
    sumwave1 = 0
    sumwave2 = 0
    
    # enumerate data - lower bin first 
    for i, val in enumerate(wavelength):
        if val >= lower_bin[0] and val <= lower_bin[1]:
            sumwave1 +=  intensity[i]
            count1 += 1
    # enumerate data - upper bin second 
    for i, val in enumerate(wavelength):
        if val >= upper_bin[0] and val <= upper_bin[1]:
            sumwave2 = sumwave2 + intensity[i]
            count2 += 1

    return ((sumwave1 / count1) / (sumwave2 / count2)) * 100

def fitfunc(t, k, MTV, CaO):
        #'Function that returns Ca computed from an ODE for a rate k, moles MTV, and IC CaO'
        #Ca0 = Ca_data[0]
        # the lambda function syntax is sometimes more compact and readable
        Casol = odeint(lambda Ca, t: -k * (2 * Ca - MTV), CaO, t)
        #Casol = odeint(myode, CaO, t) ## or call myode as a nested function above
        return Casol[:,0]

def performCalibration(rootdir):
    base_files = os.listdir(rootdir)

    Spectraldata = []  #extract the result from each loop iteration / directory
    Timedata = []

    all_data_proxy = []

    for bf in base_files:

        if not bf.startswith('.'): #and os.path.isfile(os.path.join(rootdir, bf)):
            
            level1_files = os.listdir(rootdir+'/'+bf)  #allows file to be read

            fitdata1 = []
            timestamp = []
            print(bf)  # print title

            for l1f in level1_files:



            # get spectra and sample info and break specra into wavelength and intensity as float
                spectra, sample_info = extractSpectraFromFile(rootdir+'/'+bf+'/'+l1f)
                wavelength, intensity = convertSpectraToFloat(spectra)

                # sample name, get time stamp
                #sample = "HR4C15671_09-25-48-792.txt" ## use this to override sample name
                sample = l1f
                #print ('sample label:', sample)
                hrs, mins, secs = map(int, sample.split("_")[1].split("-")[:3])
                timestamp.append((hrs * 60 * 60) + (mins * 60) + (secs))

                normdata1 = normalizeAndReduce(wavelength, intensity, lower_bin=(240, 270), upper_bin=(500, 550))

                #print(normdata1)
                fitdata1.append(normdata1)

            ## Perform fit to 1st order ODE
            # given data we want to fit and shift time to zero
            timestamp[:]=[x-timestamp[0] for x in timestamp]
            tspan = timestamp
            Ca_data = fitdata1

            #fit and covariance of fit with initial guesses p0 
            #k_fit, kcov = curve_fit(fitfunc, tspan, Ca_data, p0=[0.00005,0.1,62]) #p0=[k, MTV, CaO]
            #perr = np.sqrt(np.diag(kcov)) # 1 standard deviation
            #print ('k value is ', k_fit)
            #print ('first point in data is ', Ca_data[0])
            #print ('covariance value is ', kcov)
            #print ('the variance on k and CaO (diagonals) are ', perr)
            #print(tspan, Ca_data)

            #tfit = np.linspace(0,tspan[-1],num = 1000)
            #fit = fitfunc(tfit, k_fit[0],k_fit[1],k_fit[2])
            #k_fit[0] = k,  k_fit[1] = MTV,  k_fit[2] = CaO
            
            plt.plot(tspan, Ca_data, 'ro', label='data')
            #plt.plot(tfit, fit, 'b-', label='fit')
            #plt.legend(loc='best')
            plt.ylabel('Abs a.u.')
            plt.xlabel('Time, seconds')
            plt.show()
            Spectraldata.append(Ca_data) # update list for both time and data
            Timedata.append(tspan)

            all_data_proxy.append({'File name': bf, 
                            'Timestamps': tspan,
                            'Absorbance': Ca_data, 
                            'Average of run': np.mean(Ca_data),
                            'Average of run error': (np.std(Ca_data, ddof=1))/len(Ca_data)})  

            all_calibdata = pd.DataFrame(all_data_proxy) # this is where the dataframe is completed.
    return all_calibdata


def applyCalibration(source, sink, all_data, raw_intensity, calibintercept, calibslope):
    eq = (sink + source)/2

    solved_data = []

    for i in range (0,1):   # the second number is how many folders there are in the directory
        raw_intensity = np.array(all_data['Absorbance'][i])
        time = np.array(all_data['Timestamps'][i])
        name = all_data['File name'][i]
        conc_solved = (raw_intensity - calibintercept) / calibslope
        conc_normalised = (conc_solved) / eq
        solved_data.append({'File name': name, 
                          'Timestamps': time,
                          'Concentration': conc_solved, 
                          'Normalised concentration': conc_normalised}) 
        plt.scatter(time, conc_solved)
        plt.ylabel('[Acceptor Paracetamol] / \u03BCM')
        plt.xlabel('Time / seconds')
        plt.grid(color='black', which='major', axis='y', linestyle='solid')

    all_solveddata = pd.DataFrame(solved_data) # this is where the dataframe is completed.
    return all_solveddata  


def filterdatframe(dataframe, count, upper, lower):
    conc = np.array(dataframe['Concentration'][count])
    time = np.array(dataframe['Timestamps'][count])
    x = np.array(list(zip(conc, time)))
    conc2 = x[:,0]
    fp = x[conc2<upper]   #enter y upper limit here
    fpy = fp[:,0]
    filtered_perm = fp[fpy>lower]   #enter y lower limit here
    filtered_pconc = filtered_perm[:,0]
    filtered_ptime = filtered_perm[:,1]
    linearP = polyfit(filtered_ptime, filtered_pconc, 1, cov=True) # 1 indicates linear fit
    Pparams=linearP[0]
    pm=Pparams[0]
    pc=Pparams[1]
    return linearP, pm, pc, filtered_pconc, filtered_ptime  #returns the m and c


def performPermCalc(volume, volume_se, average_DIBarea, average_DIBarease, donor_0conc, Eq, conc_error):
    # this is where the dataframe is built
    dataframe = []  
    for i in range (0, 1):            # no. of folders we are looking in
        values = filterdatframe(all_solveddata,i,500,0)  # this will vary the timestamps in each file as the point at which each files reaches 0.6mM
        conc = values[3]
        filt_time = values[4]
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(filt_time, conc) # perform scipy statistics on fit
        permvalues = []
        permerrorvalues = []
        timevalues = []
        for x in conc, filt_time:
            concA = conc[range(0,8)]    # adjust this value for number of points within the conc range
            concD = donor_0conc - concA
            timefor = filt_time[range(0,8)]
            Pvalue = PermeabilityRheoDIB(volume,average_DIBarea,Eq,concA,slope)
            Perror = Propagation_of_P_error(Pvalue,volume_se,volume,average_DIBarease,
                                            average_DIBarea,std_err,slope,conc_error,Eq,concA)
        permvalues.append(Pvalue)
        permerrorvalues.append(Perror)
        timevalues.append(timefor)
        dataframe.append({'File name': all_solveddata['File name'][i], 
                          'Permeability': np.mean(permvalues),
                          'Standard error of P': (np.std(permvalues, ddof=1))/len(timevalues),
                          'Permeability propagated error': np.mean(permerrorvalues),
                          'Standard error of propagated error': (np.std(permerrorvalues,ddof=1))/len(timevalues)})   # this extracts all the data into the dataframe

        #Plotting for each iteration

        plt.errorbar(timevalues, permvalues, yerr=permerrorvalues, linestyle = 'none')
        plt.scatter(timevalues, permvalues, )
        #plt.yscale('log')
        #plt.ylim(0.0001,0.002)
        plt.ylabel('Permeability (cm/s)')
        plt.xlabel('Time (s)')
        plt.title(all_solveddata['File name'][i])
        plt.show()
        print ('The values for permeability are:', (permvalues))
        print('The error for each permeability is:', permerrorvalues)
        print ('The R2 value for the linear fit between the upper and lower concentration bounds is:', r_value)
        print ('The SE value for the linear fit between the upper and lower concentration bounds is:', std_err)

    all_data = pd.DataFrame(dataframe) # this is where the dataframe is completed.
    return all_data