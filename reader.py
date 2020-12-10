import os
import numpy as np
from constants import *
import h5py

# def read_temp(file):
#     '''
#     Read in data for initial temperatures

#     :param file: name of the file which holds the temperature data

#     :return input_temp: temperature vector from a specified csv file
#     :return input_year_temp: corresponding time vector (in years)
#     '''
#     # spot = os.path.dirname(sys.argv[0])
#     spot = os.getcwd()

#     FID_temp        = os.path.join(spot, file)
#     data_temp       = np.genfromtxt(FID_temp, delimiter=',')
#     input_year_temp = data_temp[0, :]
#     input_temp      = data_temp[1, :]
#     if input_temp[0] < 0.0:
#         input_temp = input_temp + K_TO_C

#     return input_temp, input_year_temp

# def read_bdot(file):
#     '''
#     Read in data for initial accumulation rates

#     :param file: name of the file which holds the accumulation rate data

#     :return input_bdot: accumulation rate vector from a specified csv file
#     :return input_year_bdot: corresponding time vector (in years)
#     '''

#     spot = os.getcwd()

#     FID_bdot        = os.path.join(spot, file)
#     data_bdot       = np.genfromtxt(FID_bdot, delimiter=',')
#     input_year_bdot = data_bdot[0, :]
#     input_bdot      = data_bdot[1, :]

#     return input_bdot, input_year_bdot

def read_input(filename):
    '''
    Read in data from csv input files

    :param file: name of the file which holds the accumulation rate data

    :return input_bdot: accumulation rate vector from a specified csv file
    :return input_year_bdot: corresponding time vector (in years)
    '''

    spot = os.getcwd()

    FID        = os.path.join(spot, filename)
    data       = np.loadtxt(FID, delimiter=',') #changed 3/6/17 to loadtxt from genfromtxt; much faster
    input_year = data[0, :]
    input_data = data[1, :]

    return input_data, input_year


# def read_input2(file):
#     '''
#     Read in data for initial accumulation rates

#     :param file: name of the file which holds the accumulation rate data

#     :return input_bdot: accumulation rate vector from a specified csv file
#     :return input_year_bdot: corresponding time vector (in years)
#     '''

#     spot = os.getcwd()

#     FID        = os.path.join(spot, file)
#     data       = np.genfromtxt(FID, delimiter=',')
#     input_year = data[0, :]
#     input_data = data[1, :]

#     return input_data, input_year

def read_snowmelt(file):
    '''
    Read in data for initial melt rates

    :param file: name of the file which holds the accumulation rate data

    :return input_bdot: accumulation rate vector from a specified csv file
    :return input_year_bdot: corresponding time vector (in years)
    '''

    spot = os.getcwd()

    FID_melt        = os.path.join(spot, file)
    data_melt       = np.genfromtxt(FID_melt, delimiter=',')
    input_year_melt = data_melt[0, :]
    input_melt      = data_melt[1, :]

    return input_snowmelt, input_year_snowmelt

def read_snowmelt(file):
    '''
    Read in data for initial melt rates

    :param file: name of the file which holds the accumulation rate data

    :return input_bdot: accumulation rate vector from a specified csv file
    :return input_year_bdot: corresponding time vector (in years)
    '''

    spot = os.getcwd()

    FID_melt        = os.path.join(spot, file)
    data_melt       = np.genfromtxt(FID_melt, delimiter=',')
    input_year_melt = data_melt[0, :]
    input_melt      = data_melt[1, :]

    return input_snowmelt, input_year_snowmelt

# def read_init(folder):
#     '''
#     Read in data for initial depth, age, density, and temperature to run the model without spin

#     :param folder: the folder containing the files holding depth, age, density, and temperature

#     :return initDepth: initial depth vector from a specified csv file
#     :return initAge: initial age vector from a specified csv file
#     :return initDensity: initial density vector from a specified csv file
#     :return initTemp: initial temperature vector from a specified csv file
#     '''

#     densityPath = os.path.join(folder, 'densitySpin.csv')
#     tempPath    = os.path.join(folder, 'tempSpin.csv')
#     agePath     = os.path.join(folder, 'ageSpin.csv')
#     depthPath   = os.path.join(folder, 'depthSpin.csv')

#     initDepth   = np.genfromtxt(depthPath, delimiter = ',')
#     initAge     = np.genfromtxt(agePath, delimiter = ',' )
#     initDensity = np.genfromtxt(densityPath, delimiter = ',')
#     initTemp    = np.genfromtxt(tempPath, delimiter = ',')

#     return initDepth, initAge, initDensity, initTemp

# def read_init(folder, resultsFileName):
#     '''
#     Read in data for initial depth, age, density, and temperature to run the model without spin

#     :param folder: the folder containing the files holding depth, age, density, and temperature

#     :return initDepth: initial depth vector from a specified csv file
#     :return initAge: initial age vector from a specified csv file
#     :return initDensity: initial density vector from a specified csv file
#     :return initTemp: initial temperature vector from a specified csv file
#     '''

#     f5 = h5py.File(os.path.join(folder, resultsFileName),'r')

#     initDensity = f5['densitySpin'][:]
#     initAge = f5['ageSpin'][:]
#     initDepth = f5['depthSpin'][:]
#     initTemp = f5['tempSpin'][:]

#     # densityPath = os.path.join(folder, 'densitySpin.csv')
#     # tempPath    = os.path.join(folder, 'tempSpin.csv')
#     # agePath     = os.path.join(folder, 'ageSpin.csv')
#     # depthPath   = os.path.join(folder, 'depthSpin.csv')

#     # initDepth   = np.genfromtxt(depthPath, delimiter = ',')
#     # initAge     = np.genfromtxt(agePath, delimiter = ',' )
#     # initDensity = np.genfromtxt(densityPath, delimiter = ',')
#     # initTemp    = np.genfromtxt(tempPath, delimiter = ',')

#     return initDepth, initAge, initDensity, initTemp

def read_init(folder, resultsFileName, varname):
    '''
    Read in data for initial depth, age, density, and temperature to run the model without spin

    :param folder: the folder containing the files holding depth, age, density, and temperature

    :return initDepth: initial depth vector from a specified csv file
    :return initAge: initial age vector from a specified csv file
    :return initDensity: initial density vector from a specified csv file
    :return initTemp: initial temperature vector from a specified csv file
    '''

    f5 = h5py.File(os.path.join(folder, resultsFileName),'r')
    init_value = f5[varname][:]


    # initAge = f5['ageSpin'][:]
    # initDepth = f5['depthSpin'][:]
    # initTemp = f5['tempSpin'][:]

    # initDensity = f5['densitySpin'][:]
    # initAge = f5['ageSpin'][:]
    # initDepth = f5['depthSpin'][:]
    # initTemp = f5['tempSpin'][:]

    # densityPath = os.path.join(folder, 'densitySpin.csv')
    # tempPath    = os.path.join(folder, 'tempSpin.csv')
    # agePath     = os.path.join(folder, 'ageSpin.csv')
    # depthPath   = os.path.join(folder, 'depthSpin.csv')

    # initDepth   = np.genfromtxt(depthPath, delimiter = ',')
    # initAge     = np.genfromtxt(agePath, delimiter = ',' )
    # initDensity = np.genfromtxt(densityPath, delimiter = ',')
    # initTemp    = np.genfromtxt(tempPath, delimiter = ',')

    return init_value

class SigmacoReader():

    def __init__(self):
        """

        """
        self.rho_co = 804.3
        return

    def _function_to_map(self, z_coarse, rho_coarse, sigma_coarse, z_co, dz):
        """

        """
        z_fine = np.arange(z_co-5, z_co+5, dz)
        rho_fine = np.interp(z_fine, z_coarse, rho_coarse)
        sigma_fine = np.interp(z_fine, z_coarse, sigma_coarse)
        sigma_co = sigma_fine[rho_fine>=self.rho_co][0]

        return sigma_co


    def read_sigma_co(self, filepath, interp_z = 1, save_ascii = False):
        """

        """

        f = h5py.File(filepath,'r')

        z = f["depth"][1:]
        rho = f["density"][1:]
        temperature = f["temperature"][1:]
        age = f["age"][1:]
        climate = f["Modelclimate"][1:]
        iso_sigmaD = f["iso_sigmaD"][1:]
        iso_sigma18 = f["iso_sigma18"][1:]
        iso_sigma17 = f["iso_sigma17"][1:]

        model_time = np.array(([a[0] for a in z[:]]))
        temp_forcing = climate[:,2]
        accum_forcing = climate[:,1]

        depth_co = np.array(([z[i][1:][rho[i][1:]>=self.rho_co][0] for i in range(len(z))]))
        depth_co_fine = np.array((list(map(self._function_to_map, z, rho, z, depth_co, np.zeros(np.size(depth_co))+interp_z))))
        age_co_fine = np.array((list(map(self._function_to_map, z, rho, age, depth_co, np.zeros(np.size(depth_co))+interp_z))))
        sigma_co_17_fine = np.array((list(map(self._function_to_map, z, rho, iso_sigma17, depth_co, np.zeros(np.size(depth_co))+interp_z))))
        sigma_co_18_fine = np.array((list(map(self._function_to_map, z, rho, iso_sigma18, depth_co, np.zeros(np.size(depth_co))+interp_z))))
        sigma_co_D_fine = np.array((list(map(self._function_to_map, z, rho, iso_sigmaD, depth_co, np.zeros(np.size(depth_co))+interp_z))))

        return {"model_time": model_time, "depth_co": depth_co_fine,
        "age_co": age_co_fine, "sigmaD_co": sigma_co_D_fine, "sigma18_co": sigma_co_18_fine,
        "sigma17_co": sigma_co_17_fine}
