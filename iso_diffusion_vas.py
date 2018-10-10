import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, FixedFormatter
from scipy import interpolate
from scipy import integrate
from scipy import stats
import copy
import diffusivity_vas
import constants



class Sigma():
    """
  Class Sigma
  Performs calculations of diffusion lengths
    """
    def __init__(self, P = 1., rho_o = 330.1, rho_co = 804.3, fo = 1, f1 = 1, f_factor_o18 = "Majoube", \
        f_factor_deuterium = "Merlivat", f_factor_o17 = "Majoube"):
        """

        Sigma class instance

        Arguments:
        ------------------
            P = 1.: Pressure in atm

            rho_o = 330: density at z = 0 [kgrm-3]

            rho_i = 917., Ice density [kgrm-3]

            rho_c = 550., Densification transition density [kgrm-3]

            rho_co = 804., Close off density [kgrm-3]

            fo = 1: HL scaling factor for ko

            f1 = 1: HL scaling factor for k1

            f_factor_o18: Version of fract. factor for o18. "Majoube" or "Ellehoj"

            f_factor_o17: Version of fract. factor for o17. "Majoube" or "Ellehoj"

            f_factor_deuterium: Version of fract. factor for deuterium. "Merlivat" or "Ellehoj"


        Examples
        --------

        """
        self.P = P
        self.rho_o = rho_o/1000 #Convert to Mgrm-3
        self.rho_i = constants.RHO_I/1000 #Convert to Mgrm-3
        self.rho_c = constants.RHO_1/1000 #Convert to Mgrm-3
        self.rho_co = rho_co/1000 #Convert to Mgrm-3
        self.fo = fo
        self.f1 = f1
        self.f_factor_o18 = f_factor_o18
        self.f_factor_o17 = f_factor_o17
        self.f_factor_deuterium = f_factor_deuterium
        self.R = constants.R #Gas constant JK-1mol-1

        return



    def analytical_HL(self, rho_array = np.array((804.3,)), T = 218.15, accum_ice = 0.025):
        """
        Calculation of  diffusion lengths for a certain density\n
        Uses an analytical solution of the diffusion length equation (Johnsen2000 eq.1)\n
        assuming simple strain rate and applying a Herron-Langway firn model\n

        Arguments:\n
        ----------

        rho_array: Array of densitites to be evaluated [kgrm-3]\n
        T: Temperature [K]\n
        accum_ice: Accumulation **ice equivalent** [myr-1]\n

        Returns:\n
        --------
        A dictionary d where\n
        {"D": sigma_D, "18": sigma_o18, "17": sigma_o17}\n
        """
        accum = accum_ice*constants.RHO_I_MGM
        rho_array = np.reshape(rho_array, -1)/1000 #Convert to Mgrm-3
        rho_array_max = self.rho_i/(np.sqrt(1.3))
        rho_array_clipped = np.clip(rho_array, self.rho_o, rho_array[rho_array<rho_array_max][-1])

        sigma_sq_D = np.zeros(np.size(rho_array))
        sigma_sq_o18 = np.zeros(np.size(rho_array))
        sigma_sq_o17 = np.zeros(np.size(rho_array))


        try:
            i_up_part = np.where(rho_array_clipped <= self.rho_c)[0]
        except:
            i_up_part = np.array(())
        try:
            i_bot_part = np.where(rho_array_clipped > self.rho_c)[0]
        except:
            i_bot_part = np.array(())

        m = 18.e-6 #molar weight in Mgr
        p_sat = diffusivity_vas.P_Ice().clausius_clapeyron_Lt(T) # p_sat in Pa
        ko = self.fo*11.*np.exp(-10160./(self.R*T))
        k1 = self.f1*575.*np.exp(-21400./(self.R*T))
        alpha = 1.#HL accumulation depandancy
        beta = 0.5 ##HL accumulation depandancy


        #### Deuterium Block  #####
        alpha_D = diffusivity_vas.FractionationFactor(T).deuterium()[self.f_factor_deuterium]
        air_diffusivity_D = diffusivity_vas.AirDiffusivity(T, self.P).deuterium()*3600*24*365.25 # Convert to m2yr-1
        Z_D = m*p_sat*air_diffusivity_D/(self.R*T*alpha_D)

        if np.size(i_up_part) != 0:
            sigma_sq_D[i_up_part] = Z_D/(self.rho_i*ko*accum**alpha*rho_array_clipped[i_up_part]**2)*\
                (rho_array_clipped[i_up_part]**2- self.rho_o**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_up_part]**4 - self.rho_o**4))
        else:
            pass

        if np.size(i_bot_part) != 0:
            sigma_sq_D[i_bot_part] = Z_D/(k1*accum**beta*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (rho_array_clipped[i_bot_part]**2 - self.rho_c**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_bot_part]**4\
                - self.rho_c**4)) + Z_D/(ko*accum**alpha*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (self.rho_c**2 - self.rho_o**2 - 1.3/(2*self.rho_i**2)*(self.rho_c**4 - self.rho_o**4))
        else:
            pass


        #### Oxygen 18 Block  ####

        alpha_o18 = diffusivity_vas.FractionationFactor(T).o18()[self.f_factor_o18]
        air_diffusivity_o18 = diffusivity_vas.AirDiffusivity(T, self.P).o18()*3600*24*365.25
        Z_o18 = m*p_sat*air_diffusivity_o18/(self.R*T*alpha_o18)

        if np.size(i_up_part) != 0:
            sigma_sq_o18[i_up_part] = Z_o18/(self.rho_i*ko*accum**alpha*rho_array_clipped[i_up_part]**2)*\
                (rho_array_clipped[i_up_part]**2- self.rho_o**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_up_part]**4 - self.rho_o**4))
        else:
            pass

        if np.size(i_bot_part) != 0:
            sigma_sq_o18[i_bot_part] = Z_o18/(k1*accum**beta*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (rho_array_clipped[i_bot_part]**2 - self.rho_c**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_bot_part]**4\
                - self.rho_c**4)) + Z_o18/(ko*accum**alpha*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (self.rho_c**2 - self.rho_o**2 - 1.3/(2*self.rho_i**2)*(self.rho_c**4 - self.rho_o**4))
        else:
            pass

        #### Oxygen 17 Block  ####

        alpha_o17 = diffusivity_vas.FractionationFactor(T).o17()[self.f_factor_o17]
        air_diffusivity_o17 = diffusivity_vas.AirDiffusivity(T, self.P).o17()*3600*24*365.25
        Z_o17 = m*p_sat*air_diffusivity_o17/(self.R*T*alpha_o17)

        if np.size(i_up_part) != 0:
            sigma_sq_o17[i_up_part] = Z_o17/(self.rho_i*ko*accum**alpha*rho_array_clipped[i_up_part]**2)*\
                (rho_array_clipped[i_up_part]**2- self.rho_o**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_up_part]**4 - self.rho_o**4))
        else:
            pass

        if np.size(i_bot_part) != 0:
            sigma_sq_o17[i_bot_part] = Z_o17/(k1*accum**beta*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (rho_array_clipped[i_bot_part]**2 - self.rho_c**2 - 1.3/(2*self.rho_i**2)*(rho_array_clipped[i_bot_part]**4\
                - self.rho_c**4)) + Z_o17/(ko*accum**alpha*self.rho_i*rho_array_clipped[i_bot_part]**2)*\
                (self.rho_c**2 - self.rho_o**2 - 1.3/(2*self.rho_i**2)*(self.rho_c**4 - self.rho_o**4))
        else:
            pass

        sigma_D = np.sqrt(sigma_sq_D)*rho_array_clipped/rho_array
        sigma_o18 = np.sqrt(sigma_sq_o18)*rho_array_clipped/rho_array
        sigma_o17 = np.sqrt(sigma_sq_o17)*rho_array_clipped/rho_array
        return {"D": sigma_D, "18": sigma_o18, "17": sigma_o17}


class Sigma2Prime():
    """
    Class SigmaPrime: performs calculations of dsigma2/dt
    """

    def __init__(self, params_dict, physical_param_dict):
        """

        Sigma2Prime class instance

        Arguments:
        ------------------
            params_dict: the dictionary of parameters from the json file of CFM
            physical_param_dict: the PhysParams dictionary as defined in the time_evolve method in
                                 firn_density_spin.py and firn_density_nospin.py


        Examples
        --------

        """

        self.cd = params_dict
        for k,v in physical_param_dict.items():
            setattr(self,k,v)

        return

    def dsigma2_dt(self, drho_dt, iso_sigma_dict):
        """
        Calculates of the sigma2 derivative with time (dsigma2/dt) for O17, O18, D\n
        It uses dsigma2/dt = 2*(D_firn - 1/rho*drho/dt*sigma2) where D_firn the firn diffusivity\n
        drho/dt the densification rate and sigma2 the diffusion length\n

        Arguments:\n
        ----------\n
        drho_dt: the densification rate array (calculated with CFM within the time_evolve method)\n
        iso_sigma_dict: the diffusion length dictionary {"D": sigma_D, "18": sigma_o18, "17": sigma_o17}\n

        Returns:\n
        --------
        A dictionary d where\n
        {"D": dsigmaD_2_dt, "18": dsigma18_2_dt, "17": dsigma17_2_dt}
        """
        sigmaD_2 = iso_sigma_dict["D"]**2
        sigma18_2 = iso_sigma_dict["18"]**2
        sigma17_2 = iso_sigma_dict["17"]**2

        diffusivity_inst = diffusivity_vas.FirnDiffusivityFast(rho = self.rho, rho_co = self.cd["rho_co_iso"], \
            T = self.Tz, P = self.cd["P_atm"])

        dsigmaD_2_dt = 2*(diffusivity_inst.deuterium(f_factor_version = "Merlivat") - drho_dt*sigmaD_2/self.rho)
        dsigma18_2_dt = 2*(diffusivity_inst.o18(f_factor_version = "Majoube") - drho_dt*sigma18_2/self.rho)
        dsigma17_2_dt = 2*(diffusivity_inst.o17(f_factor_version = "Majoube") - drho_dt*sigma17_2/self.rho)

        # dsigmaD_2_dt = diffusivity_inst.deuterium(f_factor_version = "Merlivat")/iso_sigma_dict['D'] - iso_sigma_dict['D']*drho_dt/self.rho
        # dsigma18_2_dt = diffusivity_inst.o18(f_factor_version = "Majoube")/iso_sigma_dict['18'] - iso_sigma_dict['18']*drho_dt/self.rho
        # dsigma17_2_dt = diffusivity_inst.o17(f_factor_version = "Majoube")/iso_sigma_dict['17'] - iso_sigma_dict['17']*drho_dt/self.rho

        return {"D": dsigmaD_2_dt, "18": dsigma18_2_dt, "17": dsigma17_2_dt}
