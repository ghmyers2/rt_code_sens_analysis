############
# CORRECT VERSION
############

# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 23:42:57 2021

@author: ghmyers2

This code will read the input generated by inputGenRandom or inputGenGrid
then drive the radiative transfer code to do the calculation

last modified: 08/20/2022
"""

from interpScatMatr import interpScatProp as interp
import os, time, sys, re, subprocess
import pandas as pd
import numpy as np
from shutil import copyfile
from netCDF4 import Dataset

sys.path.append("../lib")
# from altToPres import altToPressure
import createInput as info
import scatteringFiles as scat
import rtcOutput as rtc
import scipy
import matplotlib.pyplot as plt
import math


REFRAC_INDEX_DICT = {
    411.40854 : (1.3184161, 2.79E-11),
    469.2736  : (1.3145687, 1.99E-10),
    554.93747 : (1.3108153, 2.59E-09),
    670.02308 : (1.3075744, 1.92E-08),
    861.89361 : (1.3038436, 2.29E-07),
    863.51    : (1.30382, 2.35000e-07),
    959.35478 : (1.3021846, 7.53E-07),
    1588.5468 : (1.2896261, 0.000317418),
    1884.3153 : (1.279063, 0.000308622),
    2265.5691 : (1.2566663, 0.000256603),
    2111.9325 : (1.2689038, 0.00073104962)
}

def prepare_snow_iop(ar, d, reff, f_col, scatt_dir):
    ## interpolate and mix to get the snow IOP
    # ar: aspect ratio for columns
    # d: roughness parameter
    # reff: effective radius in microns
    # f_col: the fraction of columns in the mixture

    # interpolate to get the IOP for columns
    interpolant_1 = [ar, d, reff]
    # print("AR, D, REFF", interpolant_1)
    c = interp(scatt_dir, wl, interpolant_1)
    ssc_prop_interp_arr_col, scattering_matrix_interp_arr_col = c.interp_NDLinear()

    # interpolate to get the IOP for plates
    interpolant_2 = [1. / ar, d, reff]
    c = interp(scatt_dir, wl, interpolant_2)
    ssc_prop_interp_arr_pla, scattering_matrix_interp_arr_pla = c.interp_NDLinear()

    # mix the IOPs with given f_col
    ssc_prop_arr = f_col * np.array(ssc_prop_interp_arr_col) + \
                   (1. - f_col) * np.array(ssc_prop_interp_arr_pla)
    ssc_prop_arr = list(ssc_prop_arr)

    scattering_matrix_arr = []
    for i in range(len(scattering_matrix_interp_arr_col)):
        scattering_matrix_element = f_col * scattering_matrix_interp_arr_col[i] + \
                                    (1. - f_col) * scattering_matrix_interp_arr_pla[i]
        scattering_matrix_arr.append(scattering_matrix_element)

    scatt_angles = c.scatt_angles

    return ssc_prop_arr, scatt_angles, scattering_matrix_arr


if __name__ == "__main__":

    #DESCRIPTIVE VARIABLES OF THE SNOW PACK AND SOOT IMPURITIES
    mu = 0.046  # mu of impurity
    sigma = 1.5  # sigma of impurity
    mu_new = 0.1100
    sigma_new = 0.380
    # f_imp = 1.0  # ppmw
    rho_snow = 0.1  # g/cm^3
    rho_imp = 2.0  # g/cm^3
    nr_imp = 1.80  #
    ni_imp = 0.6  #

    dep_2 = 0.98  # meters
    dep_1 = 0.02

    reference_wl = 554.93747
    reference_wl_str = '555'

    wl_list = [554.93747, 469.2736, 670.02308, 863.51, 1588.5468, 2111.9325, 2265.5691]  # 1884.3153,
    wl_list_str = ['555', '470', '670', '865', '1589', '2112', '2266']  # '1884',
    scatt_dir_list = ['/home/accurt/Downloads/ice_reff_555/', '/home/accurt/Downloads/ice_reff_469/', '/home/accurt/Downloads/ice_reff_670/',
                      '/home/accurt/Downloads/ice_reff_865/', '/home/accurt/Downloads/ice_reff_1589/',
                      '/home/accurt/Downloads/ice_reff_2112/',
                      '/home/accurt/Downloads/ice_reff_2266/']  # '/home/accurt/Downloads/ice_reff_1884/',

    ALBEDO = [0.00] * len(wl_list)  # same number for each wavelength [-1.30]
    SRFFILELIST = ["oceanl2264test"] * len(wl_list)  # surface_filename rename ["oceanl2250test"]
    NTYPE = 3
    NPERS, NLAM, NLAYER, IPRT, IREMV = 2, len(wl_list), 8, 0, 1  # CHECK IF JUST 2 LAYERS OR INCLUDE ATMOSPHERE
    A = [mu_new, 0.0, 0.0]  # add
    B = [sigma_new, 0.0, 0.0]  # add
    NR = [[nr_imp], [1.0], [1.0]] * len(wl_list)  # [[1.450000, 1.410000]]
    NI = [[ni_imp], [0.0], [0.0]] * len(wl_list)  # [[0.00E-2, 0.00E-3]]
    R1 = 0.0
    R2 = 10.0
    NSD = [3, 6, 6]
    DELP = [0.11772121669935132, 0.0024023581127039506, 218.1123876207905, 178.4078979709592, 144.42710392617516,
            207.17752947768372, 265.00435142957934, 0.000606]
    TAUREP = 1

    NGAS = 0
    TAUABS = []

    parameter_file_name = '/home/accurt/Shared/SOBOL/param_values_SOBOL.txt'
    input_data = pd.read_csv(parameter_file_name, header=None, delimiter='\s+')
#     dict_param = {0: "sza", 1: "ar_1", 2: "d_1",
#                   3: "f_1", 4: "reff_1", 5: "ar_2", 6: "d_2",
#                   7: "f_2", 8: "reff_2", 9: 'soot_1', 10: "soot_2", 11: 'aod'}
    
    dict_param = {0:'ar_1', 1:'d_1', 2:'f_1', 3:'Reff_1', 4:'ar_2', 5:'d_2', 6:'f_2', 7:'Reff_2',
                  8:'thick_top', 9:'dens_top', 10:'dens_bottom', 11:'aod', 12:'soot_1', 13:'soot_2'}

    input_data.rename(columns=dict_param, inplace=True)

    num_cases = len(input_data)
    # num_cases = 5

    for icase in range(num_cases):
        RSPFILELIST = []
        REFRACFILELIST = []

        # GRABBING THE FREE PARAMETER VALUES FROM THE .txt FILE GENERATED BY inputGenRandom.py

        df = input_data.iloc[[icase]]
        raz = 0.0
        sza = 45.0
        
        ar_1 = df.ar_1.values[0]
        d_1 = df.d_1.values[0]
        f_1 = df.f_1.values[0]
        reff_1 = df.Reff_1.values[0]
        
        ar_2 = df.ar_2.values[0]
        d_2 = df.d_2.values[0]
        f_2 = df.f_2.values[0]
        reff_2 = df.Reff_2.values[0]
        
        aod = df.aod.values[0]
        soot_1 = df.soot_1.values[0]
        soot_2 = df.soot_2.values[0]


#         print("__________________________________")
#         print(f"___CURRENT CASE NUMBER: {icase}___")
#         # sza, raz, ar_1, d_1, f_1, reff_1, soot_1, ar_2, d_2, f_2, reff_2, soot_2, aod
#         print(f"sza:{sza}, raz:{raz}, ar_1:{ar_1}, d_1:{d_1}, f_1:{f_1}, reff_1:{reff_1}")
#         print(f"soot_1:{soot_1}, ar_2:{ar_2}, d_2:{d_2}, f_2:{f_2}")
#         print(f"reff_2:{reff_2}, soot_2:{soot_2}, aod:{aod}")
#         print("__________________________________")

        # GENERATING THE INTERPOLATED SCATTERING MATRICES FOR EACH WL
        for i in range(len(wl_list)):
            wl = wl_list[i]
            wl_str = wl_list_str[i]
            scatt_dir = scatt_dir_list[i]
            output_file = 'SOBOL_rtc_output_mult' + str(icase) + f"_{wl_str}"
            RSPFILELIST.append(output_file)

            # find and write snow IOP for the first layer
            ssc_prop_arr, scatt_angles, scattering_matrix_arr = prepare_snow_iop(ar_1, d_1, reff_1, f_1,
                                                                                 scatt_dir=scatt_dir)
            output_dir = '/home/accurt/rt_code/rt_code/info/'

            output_name_1 = f"a_snow_iop_1_mult_{icase}_{wl_str}.txt"
            snow_iop_1_name = output_dir + output_name_1

            print("WAVELENGHT : ", wl_str, "Eff Radius :", reff_1, reff_2)

            interp.write_snow_ssc_prop(snow_iop_1_name, ssc_prop_arr, scatt_angles, scattering_matrix_arr, Reff=reff_1,
                                       NR=REFRAC_INDEX_DICT[wl][0], NI=REFRAC_INDEX_DICT[wl][1], Lambda=wl / 1000)

            # find and write snow IOP for the second layer
            ssc_prop_arr, scatt_angles, scattering_matrix_arr = prepare_snow_iop(ar_2, d_2, reff_2, f_2,
                                                                                 scatt_dir=scatt_dir)

            output_name_2 = f"a_snow_iop_2_mult_{icase}_{wl_str}.txt"
            snow_iop_2_name = output_dir + output_name_2

            interp.write_snow_ssc_prop(snow_iop_2_name, ssc_prop_arr, scatt_angles, scattering_matrix_arr, Reff=reff_2,
                                       NR=REFRAC_INDEX_DICT[wl][0], NI=REFRAC_INDEX_DICT[wl][1], Lambda=wl / 1000)

            REFRACFILELIST.append("scatMatrPlaceholder")
            REFRACFILELIST.append(snow_iop_2_name)
            REFRACFILELIST.append(snow_iop_1_name)

        # CALCULATING REFERENCE MEASUREMENTS FOR BOTH IMPURITIES AND SNOW

        # mu has a unit (units of length), must have same units as wavelength

        imp = scat.getIOPs.sizedis_spher(mu_new, sigma_new)
        tau_soot_2 = scat.getIOPs.calcTau_imp(imp, soot_2, rho_snow, rho_imp, reference_wl / 1000., nr_imp, ni_imp,
                                              dep_2)
        tau_soot_1 = scat.getIOPs.calcTau_imp(imp, soot_1, rho_snow, rho_imp, reference_wl / 1000., nr_imp, ni_imp,
                                              dep_1)

        ALAM = [wl / 1000. for wl in wl_list]
        mu0 = np.cos(np.radians(sza))
        phi = raz

        snow_iop_ref_2_name = f'/home/accurt/rt_code/rt_code/info/a_snow_iop_2_mult_{icase}_{reference_wl_str}.txt'
        snow_iop_ref_1_name = f'/home/accurt/rt_code/rt_code/info/a_snow_iop_1_mult_{icase}_{reference_wl_str}.txt'
        tau_snow_bottom = scat.calcTau(snow_iop_ref_2_name, rho_snow, dep_2)
        tau_snow_top = scat.calcTau(snow_iop_ref_1_name, rho_snow, dep_1)

        # PACKAGING REFERENCE VARIABLES INTO ARRAY OF SPECIES AMOUNTS FOR EACH PRESSURE LAYER

        NZITEMS = [[tau_soot_2, tau_soot_1, aod, 0.00E-02, 0.00E-01, 0.00E-01, 0.00E-02, 0.00E-01],
                   [tau_snow_bottom, 0.00E-01, 0.00E-01, 0.00E-02, 0.00E-01, 0.00E-01, 0.00E-02, 0.00E-01],
                   [0.00E-02, tau_snow_top, 0.00E-01, 0.00E-02, 0.00E-01, 0.00E-01, 0.00E-02,
                    0.00E-01]]  # species_amount rename

        file_name = 'lut_rewrite_' + str(icase) + '_a2_multiple_wl.info'
        NTAU = 24
        NTAU2 = 24
        nPhi = 256
        nGauss = 24

        info.createInfo(fileName=file_name, A=A, B=B, R1=R1, R2=R2, NSD=NSD,
                        mu0=float(mu0), phi=float(phi), NTYPE=NTYPE, NLAM=NLAM, ALAM=ALAM, NLAYER=NLAYER,
                        NGAS=NGAS, IPRT=IPRT, IREMV=IREMV,
                        ALBEDO=ALBEDO, SRFFILELIST=SRFFILELIST, RSPFILELIST=RSPFILELIST,
                        NR=NR, NI=NI, REFRACFILELIST=REFRACFILELIST,
                        DELP=DELP, NZITEMS=NZITEMS, TAUABS=TAUABS, NTAU=NTAU, NTAU2=NTAU2, nPhi=nPhi, nGauss=nGauss)

        print("______________________________")
        print("AODs of soot immersed in snow", tau_soot_2, tau_soot_1)
        print("AODs of snow layers", tau_snow_bottom, tau_snow_top)
        # quit()

        os.chdir("/home/accurt/rt_code/rt_code/")
        subprocess.Popen("make all", shell=True, cwd="/home/accurt/rt_code/rt_code/").wait()
        subprocess.Popen(f"./vec_generate_obs info/{file_name} 0 1", shell=True,
                         cwd="/home/accurt/rt_code/rt_code/").wait()

        #### NOW CHECK THAT THE RTC RAN PROPERLY AND GENERATD AN OUTPUT FILE
        output_test = RSPFILELIST[0]
        currentOutput = rtc.rtcOutput(output_test + '.rsp')
        Ival = currentOutput.RV11 + currentOutput.RZ11
        Qval = currentOutput.RV21 + currentOutput.RZ21
        Uval = currentOutput.RV31 + currentOutput.RZ31
        polRef = np.sqrt(Qval ** 2 + Uval ** 2)
        DolP = polRef / Ival

        if math.isnan(Ival[0]) or math.isnan(Qval[0]) or math.isnan(Uval[0]) or math.isnan(polRef[0]) or math.isnan(
                DolP[0]):
            print("NAN VALUE FROM RTC OUTPUT")
            print(f"I: {Ival}, Q:{Qval}, U:{Uval}")
            print(f"Rp:{polRef}, DoLP:{DolP}")

        # Now we delete the unnecessary files (scattering matrices and rsp files without .rsp suffix)

        for wl_val in wl_list_str:
            # Deleting snow scattering matrices
            print(f"REMOVING FOR CURRENT WL{wl_val}")
            output_dir = '/home/accurt/rt_code/rt_code/info/'
            output_name_1 = f"a_snow_iop_1_mult_{icase}_{wl_val}.txt"
            snow_iop_1_name = output_dir + output_name_1
            output_name_2 = f"a_snow_iop_2_mult_{icase}_{wl_val}.txt"
            snow_iop_2_name = output_dir + output_name_2
            os.remove(snow_iop_1_name)
            os.remove(snow_iop_2_name)

            # Deleting the info and output files (saving the first .info file)
            output_file = 'SOBOL_rtc_output_mult' + str(icase) + f"_{wl_val}"
            copyfile("/home/accurt/rt_code/rt_code/" + output_file + ".rsp",
                     "/home/accurt/Shared/SOBOL/" + output_file + ".rsp")

            os.remove("/home/accurt/rt_code/rt_code/" + output_file)
            os.remove("/home/accurt/rt_code/rt_code/" + output_file + ".rsp")

        # Deleting the .info file
        if icase > 0:
            os.remove("/home/accurt/rt_code/rt_code/info/" + file_name)
