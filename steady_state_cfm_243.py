from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import os
import os.path
import sys
import json
import string
import h5py
import shutil
import time
import diffusivity_vas
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# plt.rc('font', size = 12)
plt.rc('xtick', direction = 'in')
plt.rc('ytick', direction = 'in')
mylabelsize = 14
plt.rc('ytick', labelsize = mylabelsize)
plt.rc('xtick', labelsize = mylabelsize)
plt.rc('axes', labelsize = mylabelsize)

def run_st_st(configName = "./vas_simple_tests_st_st.json", temp = 243.15, \
accum = 0.20, years = 1000, physRho = "Goujon2003", **kwargs):



    #Build steady state forcing files
    f = open("./cfm_mytests/cfm_my_forcings/temp_forcing_st_st.csv", "w")
    f.write("%i,%0.2f\n" %(years*(-1), 0.0))
    f.write("%0.5f, %0.5f" %(temp, temp))
    f.close()

    f = open("./cfm_mytests/cfm_my_forcings/accum_forcing_st_st.csv", "w")
    f.write("%i,%0.2f\n" %(years*(-1), 0.0))
    f.write("%0.5f, %0.5f" %(accum, accum))
    f.close()



    with open(configName, "r") as f:
        jsonString = f.read()
        c = json.loads(jsonString)

    c["InputFileNameTemp"] = "cfm_mytests/cfm_my_forcings/temp_forcing_st_st.csv"
    c["InputFileNamebdot"] = "cfm_mytests/cfm_my_forcings/accum_forcing_st_st.csv"
    c["resultsFolder"] = "cfm_mytests/steady_state_test_results_temp/"
    c["physRho"] = physRho
    c["resultsFileName"] = "steady_state_test_" + str(np.round(temp, 2)).replace(".", "_") + "_" + str(np.round(accum, 4)).replace(".", "_")\
     + "_" + physRho + ".hdf5"

    if len(kwargs.keys()) != 0:
        for i in kwargs.keys():
            if i in c.keys():
                save_json = True
                c[i] = kwargs[i]
                print("changing json parameter %s, new value %s" %(i, c[i]))
            else: save_json = False

    with open(configName, "w") as f:
        json.dump(c, f, indent = 4, sort_keys = True)
    sys_command = "python main.py vas_simple_tests_st_st.json -n"
    os.system(sys_command)
    shutil.move(os.path.join(c["resultsFolder"], c["resultsFileName"]), "cfm_mytests/steady_state_test_results/")

    return c



def read_st_st_test():
    """

    """
    plt.close("all")
    plt.ion()
    results_folder = "./cfm_mytests/steady_state_test_results"
    list_of_files = os.listdir(results_folder)
    plt.figure(12)
    for fil in list_of_files:
        list_fil = string.split(fil, "_")
        # print list_fil

        #list of strings not allowed in file name - continues with next iteration if True
        str_list_excl = ["heat"]
        if np.any([k in fil for k in str_list_excl]):
            continue

        str_list_incl = ["hdf5", "Goujon"]
        if np.all([k in fil for k in str_list_incl]):
            print(fil)
            f = h5py.File(os.path.join(results_folder, fil))
            # print f.keys()
            z = f["depth"][:]
            rho = f["density"][:]
            temperature = f["temperature"][:]
            age = f["age"][:]
            climate = f["Modelclimate"][:]
            iso_sigmaD = f["iso_sigmaD"][:]
            iso_sigma18 = f["iso_sigma18"][:]
            iso_sigma17 = f["iso_sigma17"][:]
            plt.subplot(121)
            plt.plot(rho[-1][1:], z[-1][1:], label = list_fil[-2] + list_fil[-1])
            plt.ylim(300, 0)
            plt.subplot(122)
            plt.plot(iso_sigma18[-1][1:], z[-1][1:], label = list_fil[-2] + list_fil[-1])
            plt.ylim(300, 0)
            plt.title(str_list_incl[1])
            f.close()

    plt.show()
    return

def export_st_st_sigmas(rho_phys = "HLdynamic"):
    """

    """
    results_folder = "./cfm_mytests/steady_state_test_results"
    list_of_files = np.sort(os.listdir(results_folder))

    model_list = []
    temp_array = np.array(())
    depth_co_array = np.array(())
    accum_array = np.array(())
    sigma17_co = np.array(())
    sigma18_co = np.array(())
    sigma_D_co = np.array(())
    age_co_array = np.array(())
    rho_co = 804.3


    for fil in list_of_files:
        list_fil = string.split(fil, "_")

        #list of strings not allowed in file name - continues with next iteration if True
        str_list_excl = ["heat", "off"]
        if np.any([k in fil for k in str_list_excl]):
            continue

        #list of strings to be included in file name
        str_list_incl = ["hdf5", rho_phys, "Lamb"]
        if np.all([k in fil for k in str_list_incl]):
            print(fil)
            f = h5py.File(os.path.join(results_folder, fil))
            # print f.keys()
            z = f["depth"][-1]
            rho = f["density"][-1]
            temperature = f["temperature"][-1]
            age = f["age"][-1]
            climate = f["Modelclimate"][:]
            iso_sigmaD = f["iso_sigmaD"][-1]
            iso_sigma18 = f["iso_sigma18"][-1]
            iso_sigma17 = f["iso_sigma17"][-1]



            model_list.append(rho_phys)
            temp_array = np.append(temp_array, temperature[rho>rho_co][0])
            depth_co_array = np.append(depth_co_array, z[rho>rho_co][0])
            accum_array = np.append(accum_array, climate[-1,1])
            age_co_array = np.append(age_co_array, age[rho>rho_co][0])
            sigma17_co = np.append(sigma17_co, iso_sigma17[rho>rho_co][0])
            sigma18_co = np.append(sigma18_co, iso_sigma18[rho>rho_co][0])
            sigma_D_co = np.append(sigma_D_co, iso_sigmaD[rho>rho_co][0])
    print(temp_array)
    print(accum_array)
    print(age_co_array)

    dataout = np.transpose(np.vstack((temp_array, accum_array, depth_co_array, age_co_array, \
        sigma17_co, sigma18_co, sigma_D_co)))
    f = open("./cfm_mytests/export_st_st_sigmas.out", "w")
    f.write(model_list[-1] + "\n")
    f.write("temp\taccum\tdepth_co\tage_co\tsigma17_co\tsigma18_co\tsigmaD_co\n")
    np.savetxt(f, dataout, fmt = "%0.4f\t%0.6f\t%0.4f\t%0.4f\t%0.5f\t%0.5f\t%0.5f")
    f.close()

    return

def plot_profiles_eq_terms():
    f = h5py.File("./cfm_mytests/steady_state_test_results/steady_state_test_243_15_0_2_HLdynamic.hdf5")
    print(f.keys())
    z = f["depth"][-1][1:]
    rho = f["density"][-1][1:]
    temperature = f["temperature"][-1]
    age = f["age"][-1]
    climate = f["Modelclimate"][:]
    iso_sigmaD = f["iso_sigmaD"][-1]
    iso_sigma18 = f["iso_sigma18"][-1][1:]
    iso_sigma17 = f["iso_sigma17"][-1]
    drho_dt = f["drho_dt"][-1][1:]
    iso_dsigma2_dt_18 = f["iso_dsigma2_dt_18"][-1][1:]
    print(rho)

    plt.ion()
    plt.close("all")
    f1, axess = plt.subplots(nrows = 1, ncols = 1 , num = 567, figsize = (10,6), tight_layout = True)
    axess.plot(z, iso_dsigma2_dt_18, linewidth = 1, color = "r")
    axess.plot(z, drho_dt*iso_sigma18**2/rho, linewidth = 1, color = "b")
    axess.plot(z, iso_dsigma2_dt_18/2 + drho_dt*iso_sigma18**2/rho, linewidth = 1, color = "k")
    # axess.plot(z, iso_sigma18, linewidth = 1, color = "r")
    run_st_st(physRho = "HLdynamic", tortuosity = "johnsen2000")

    return


#
if __name__ == "__main__":
    run_st_st(physRho = "HLdynamic", tortuosity = "johnsen2000")
    plot_profiles_eq_terms()
    # os.remove("./steady_state_test_results_temp/steady_state_test_243_15_0_2_HLdynamic.hdf5")
#     t1 = time.time()
#     b = np.logspace(np.log(0.015),np.log(0.3), 30, base = np.e)
#     a = np.linspace(213.0, 253, 30)
#     print(a)
#     print(b)
#     fout = open("./f_factors.dat", "w")
#     fout.write("Temperature\talpha_D_merlivat\talpha_D_ellehoj\talpha_D_lamb\talpha_18_majoube\talpha_18_ellehoj\talpha_17_majoube\talpha_17_ellehoj\n")
#     for i in range(len(a)):
#         dif_inst = diffusivity_vas.FirnDiffusivity(rho = 550.0, P = 0.7, T = a[i])
#         dataout = np.transpose(np.vstack((a[i], dif_inst.deuterium(f_factor_version = "Merlivat"),\
#         dif_inst.deuterium(f_factor_version = "Ellehoj"), dif_inst.deuterium(f_factor_version = "Lamb"),\
#         dif_inst.o18(f_factor_version = "Majoube"), dif_inst.o18(f_factor_version = "Ellehoj"),\
#         dif_inst.o17(f_factor_version = "Majoube"), dif_inst.o17(f_factor_version = "Ellehoj"))))
#         print(dataout)
#         np.savetxt(fout, dataout, fmt = "%0.6e", delimiter = "\t")
#     #    run_st_st(temp = a[i], accum = b[i], years = 500, physRho = "HLdynamic")
#     # export_st_st_sigmas(rho_phys = "HLdynamic")
# #     # read_st_st_test()
#     # run_st_st(temp = 250.0, HbaseSpin = 2800)
# #     print("\n\nScript time: %0.2f minutes" %((time.time() - t1)/60.))
#     fout.close()
