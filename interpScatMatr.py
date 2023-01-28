import re
import sys
# sys.path.append("../")
# import lib.rtcOutput as rt
import pandas as pd
import numpy as np
import os
import time
import scipy.interpolate as interp
from os import path



class interpScatProp:
    
    ## The parameter table for interpolation, interpolant value must not exceed the range of each parameter 
    
    ## This table is for non-sperical snow particles 
    ## AR -- aspect ratio, D -- roughness parameter, Reff -- effective radius
    AR = [0.037, 0.044, 0.051, 0.06, 0.07, 0.082, 0.096, 0.112, 0.131, 0.153, 0.179, 0.209, 0.245, 0.286, 0.334, 0.391,
            0.457, 0.535, 0.625, 0.731, 0.855, 1.0, 1.169, 1.367, 1.599, 1.87, 2.187, 2.557, 2.99, 3.497, 4.089, 4.782,
            5.592, 6.539, 7.646, 8.942, 10.456, 12.227, 14.299, 16.721, 19.553, 22.865, 26.738]
    D = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
    Reff = [5.0, 7.0, 10.0, 14.0, 20.0, 28.0, 40.0, 56.0, 80.0, 113.0, 160.0, 226.0, 320.0, 640.0, 1280.0, 2560.0] 

    table_interpolation = [AR, D, Reff]
    
    def __init__(self, scatt_dir, wl, interpolant):
        
        self.scatt_dir  = scatt_dir
        ## the wavelength is an input parameter and not interpolated
        self.wl = wl     
        
        # file_name = 'FWDinput_0863_155_000.46_005_phase.dat'
        # full_name = self.scatt_dir + file_name
        # _ = self.read_scattering_properties(full_name)
        self.interpolant = interpolant
        
        index_interpolation = self.get_interp_indices(interpolant)
        # print("___INDEX INTERPOLATION____", index_interpolation)
        self.index_interpolation = index_interpolation
        # print(index_interpolation)
       
        par_arr, ssc_prop_arr, scattering_matrix_arr = self.read_scattering_file_array_snow()
        
        self.par_arr = par_arr
        self.ssc_prop_arr = ssc_prop_arr
        self.scattering_matrix_arr = scattering_matrix_arr    

    def interp_NDLinear(self):
        
        ## This routine will interpolate the single scattering properties
        ## and the scattering matrix
        
        interpolant = self.interpolant
        
        ssc_prop_interp_arr = []
        scattering_matrix_interp_arr = []
        
        # interpolate cext, csca, ssa, asf 
        for icol in range(self.ssc_prop_arr.shape[1]):
            x_arr = self.par_arr
            y_arr = self.ssc_prop_arr[:,icol]
            x = interpolant
            result = self.lin_interpl(x_arr,y_arr, x)
            ssc_prop_interp_arr.append(float(result))
        
        # interpolate the scattering matrix
        for icol in range(self.scattering_matrix_arr.shape[1]):
            x_arr = self.par_arr
            y_arr = self.scattering_matrix_arr[:,icol,:]
            x = interpolant
            result = self.lin_interpl(x_arr,y_arr, x)
            scattering_matrix_interp_arr.append(np.squeeze(result))
        
        return ssc_prop_interp_arr, scattering_matrix_interp_arr

    def read_scattering_file_array_snow(self):
        
        wl = self.wl
     
        par_arr = []
        ssc_prop_arr = []
        scattering_matrix_arr = []
        
        for ar_ind in self.index_interpolation[0]:
            for d_ind in self.index_interpolation[1]:
                for r_ind in self.index_interpolation[2]:
                    ar = self.AR[ar_ind]
                    d  = self.D[d_ind]
                    r  = self.Reff[r_ind]
                    par_arr.append([ar, d, r])
                    file_name = self.get_file_name_snow(wl, r, ar, d)
                    full_name = self.scatt_dir + file_name
                    scatt_angles, ssc_prop, scattering_matrix = self.read_scattering_properties(full_name)
                    ssc_prop_arr.append(ssc_prop)
                    scattering_matrix_arr.append(scattering_matrix)
        
        self.scatt_angles = scatt_angles
        
        return par_arr, np.array(ssc_prop_arr), np.array(scattering_matrix_arr)
    
    @staticmethod
    def read_scattering_properties(file_name, header_length=6):
        
        ## To read the snow single scattering property file 
                
        f = open(file_name, "r")
        
        lines = f.readlines()
        file_length = len(lines)
        num_scattering_angles = file_length - header_length
        
        title = lines[5]
        
        x = re.split('; |, |\*|\n|=|    ', lines[3])
        reff = float(x[1]) # effective radius
        veff = float(x[3]) # effective variance
        area = float(x[5]) # average area
        wl   = float(x[7]) # wavelength
        nr   = float(x[9]) # real refractive index 
        ni   = float(x[11]) # imaginary refractive index 
        
        x = re.split('; |, |\*|\n|=|    ', lines[4])
       
        cext = float(x[1]) # extinction cross-section 
        csca = float(x[3]) # scattering cross-section 
        ssa  = float(x[5]) # single scattering albedo
        asf = float(x[7])  # asymmetry parameter      
        ssc_prop = [cext, csca, ssa, asf] # pack the single scattering properties
        
        scatt_angles = []
        P11 = []
        P22_over_P11 = []
        P33_over_P11 = []
        P44_over_P11 = []
        P12_over_P11 = []
        P34_over_P11 = []
        
        for i in range(num_scattering_angles):
            line = lines[i + header_length]
            scatt_angles_loop = float(line.split()[0])
            scatt_angles.append(scatt_angles_loop)
            P11_loop =  float(line.split()[1])
            P11.append(P11_loop)
            P22_over_P11_loop =  float(line.split()[2])
            P22_over_P11.append(P22_over_P11_loop)
            P33_over_P11_loop =  float(line.split()[3])
            P33_over_P11.append(P33_over_P11_loop)
            P44_over_P11_loop =  float(line.split()[4])
            P44_over_P11.append(P44_over_P11_loop)
            P12_over_P11_loop =  float(line.split()[5])
            P12_over_P11.append(P12_over_P11_loop)
            P34_over_P11_loop =  float(line.split()[6])
            P34_over_P11.append(P34_over_P11_loop)
        
        scatt_angles = np.array(scatt_angles)
        P11 = np.array(P11)
        P22_over_P11 = np.array(P22_over_P11)
        P33_over_P11 = np.array(P33_over_P11)
        P44_over_P11 = np.array(P44_over_P11)
        P12_over_P11 = np.array(P12_over_P11)
        P34_over_P11 = np.array(P34_over_P11)
        
        # pack the scattering matrix
        scattering_matrix = [P11, P22_over_P11, P33_over_P11, P44_over_P11, P12_over_P11, P34_over_P11]
        
        f.close()
        
        return scatt_angles, ssc_prop, scattering_matrix
    
    @staticmethod
    def lin_interpl(x_arr, y_arr, x):
        
        lin_interp = interp.LinearNDInterpolator(x_arr, y_arr)      
        y_interp = lin_interp(x)
        
        return y_interp
    
    @staticmethod
    def get_file_name_snow(wl, r, ar, d):
        
        ## To compile the input file name for snow IOP 
        # wl: wavelength in nm
        # r: effective radius in micron
        # ar: aspect ratio
        # d: roughness
        # print("WAVLENEHTN OF IMPORTANCE", wl) jadhad
        dict_of_file_sig = {
            0.41027:"AS" + ("%06.3f" % ar) + "_R" + ("%06.2f" % r) + "_WL" + str(int(round(wl, 0))) + "_D" + "{:01.2f}".format(d) + ".dat",
            0.46913:"AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL0469"  + "_D" + "{:01.2f}".format(d) +"_K0.75" + ".dat",
            0.4692736: "AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL0469" + "_D" + "{:01.2f}".format(d) + "_K0.75" + ".dat",
            0.5549374699999999:  "AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL0555" + "_D" + "{:01.2f}".format(d) + "_K0.75" + ".dat",
            0.55496:"AS" + ("%06.3f" % ar) + "_R" + ("%06.2f" % r) + "_WL" + str(int(round(wl, 0))) + "_D" + "{:01.2f}".format(d) + ".dat",
            0.67001:"AS" + ("%06.3f" % ar) + "_R" + ("%06.2f" % r) + "_WL" + str(int(round(wl, 0))) + "_D" + "{:01.2f}".format(d) + ".dat",
            0.67002308:"AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL0670" + "_D" + "{:01.2f}".format(d) +"_K0.75" + ".dat",
            0.6700230800000001:"AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL0670" + "_D" + "{:01.2f}".format(d) +"_K0.75" + ".dat",
            0.86351:"AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL0862"  + "_D" + "{:01.2f}".format(d) +"_K0.75" + ".dat",
            0.96:"AS" + ("%06.3f" % ar) + "_R" + ("%06.2f" % r) + "_WL" + str(int(round(wl, 0))) + "_D" + "{:01.2f}".format(d) + ".dat",
            1.58886:"AS" + ("%06.3f" % ar) + "_R" + ("%06.2f" % r) + "_WL1589"  + "_D" + "{:01.2f}".format(d) + ".dat",
            1.5885468: "AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL1589"  + "_D" + "{:01.2f}".format(d) +"_K0.75" + ".dat",
            1.8843153: "AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL1884" + "_D" + "{:01.2f}".format(d) + "_K0.75" + ".dat",
            2.2655691: "AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL2266" + "_D" + "{:01.2f}".format(d) + "_K0.75" + ".dat",
            2.1119325:"AS" + ("%06.3f" % ar) + "_REFF" + ("%08.2f" % r) + "_WL2112"  + "_D" + "{:01.2f}".format(d) +"_K0.75" + ".dat",
            1.88:"AS" + ("%06.3f" % ar) + "_R" + ("%06.2f" % r) + "_WL" + str(int(round(wl, 0))) + "_D" + "{:01.2f}".format(d) + ".dat",
            2.26438:"AS" + ("%06.3f" % ar) + "_R" + ("%06.2f" % r) + "_WL" + str(int(round(wl, 0))) + "_D" + "{:01.2f}".format(d) + ".dat"
        }
        file_name = "FWDinput_" + "{:04d}".format(int(wl)) + "_" + "{:03d}".format(
                    int(r)) + "_" + "{:06.2f}".format(round(ar, 2)) + "_" + "{:03d}".format(int(100 * d)) + "_phase.dat"
        file_name = dict_of_file_sig[wl/1000]
        
        return file_name


    @classmethod
    def get_interp_indices(self, interpolant): 
        
        index_interpolation = []
        
        for i in range(len(interpolant)):
            index = self.interpolationSearch(self.table_interpolation[i], interpolant[i])
            index_interpolation.append(index)
        
        return index_interpolation
      
    @staticmethod
    def interpolationSearch(arr, x):
        
        ## search for the indices of the closest neighbors to the interpolant
        ## the interpolant needs to be within the range of the table in each dimension
        
        if x > arr[len(arr)-1] or x < arr[0]:
            print('interpolant {0} out of range, please check the table range'.format(x))
            sys.exit()
            
        if x in arr:
            
            if arr.index(x) == len(arr) - 1:
                index = [arr.index(x) - 1, arr.index(x)]
            else:
                index = [arr.index(x), arr.index(x) + 1]
            
            return index
        
        else:
            
            for i in range(len(arr)):
                if arr[i] < x and arr[(i + 1)] > x:
                    index = [i, i + 1]
                    
            return index

    @staticmethod
    def write_snow_ssc_prop(output_path, ssc_prop_arr, scatt_angles, scattering_matrix_arr, Lambda=0.46913,
                            NR=1.3038E0, NI=0.0E0, Reff=4.7E2, Veff=1.0, Area=4.48902E6):
        
        ## Write out snow single scattering properties
        
        # the following variables are fixed dummy variables in the header for now
        # Reff = 4.7e2
        # Veff = 1.0
        # Area = 4.48902e6
        # Lambda = 8.6351e-1
        # NR=1.3038E0
        # NI=0.0E0
        
        cext, csca, ssa, asf = ssc_prop_arr
        angles = scatt_angles
        P11 = scattering_matrix_arr[0]
        P22overP11 = scattering_matrix_arr[1]
        P33overP11 = scattering_matrix_arr[2]
        P44overP11 = scattering_matrix_arr[3]
        P12overP11 = scattering_matrix_arr[4]
        P34overP11 = scattering_matrix_arr[5]        
        
        file_name = os.path.basename(output_path)
        
        # write out 
        f = open(output_path, "w")
        
        f.write("Ice optical properties for given values of\n")
        f.write("effective radius (Reff), aspect ratio (AR) and distortion parameter (delta)\n")
        f.write("Source File: " + file_name + "\n")
        numSpaceVeff = 1 - len(str(Veff).split(".")[0])
        f.write("Reff= " + format(Reff, ".5E") + "    Veff=" + " " * numSpaceVeff + ("{:.4f}").format(
            Veff) + "    Area= " + format(Area, ".5E") + "    Lambda= " + format(Lambda, '.5E'))
        f.write("    NR= " + format(NR, '.5E') + "    NI= " + format(NI, '.5E') + "\n")                
        f.write("Cext= " + format(cext, ".5E") + "    Csca= " + format(csca, ".5E") + "    Pizero=" + ("{:.6f}").format(
            ssa) + "    Cosbar=" + ("{:.6f}").format(asf) + "\n")
        f.write(
            "Angle (degs)          P11          P22/P11         P33/P11         P44/P11         P12/P11         P34/P11")
        f.write("\n")
        for i in range(len(angles)):
            f.write(" " + format(angles[i], ".5E"))
            ###         Note that we must adjust the number of whitespaces if the number is negative
            ###         Thus, nspace is 5 if the number is positive, 4 if negative (to leave space for negative sign)
            nSpace = (5 if P11[i]   >= 0.0 else 4)
            f.write(" " * nSpace + format(P11[i] , ".5E"))
            nSpace = (5 if P22overP11[i] >= 0.0 else 4)
            f.write(" " * nSpace + format(P22overP11[i] , ".5E"))
            nSpace = (5 if P33overP11[i] >= 0.0 else 4)
            f.write(" " * nSpace + format(P33overP11[i] , ".5E"))
            nSpace = (5 if P44overP11[i] >= 0.0 else 4)
            f.write(" " * nSpace + format(P44overP11[i] , ".5E"))
            nSpace = (5 if P12overP11[i]  >= 0.0 else 4)
            f.write(" " * nSpace + format(P12overP11[i] , ".5E"))
            nSpace = (5 if P34overP11[i]  >= 0.0 else 4)
            f.write(" " * nSpace + format(P34overP11[i] , ".5E"))
            f.write("\n")        
        f.close()
        
        return        

if __name__ == '__main__':

    st = time.time()
    
    scatt_dir = '/Users/cyberbass/Downloads/FWD_INPUT_SNOW/'    
    interpolant_test = [0.02, 0.000, 150.0]
    wl = 863.51
    
    c = interpScatProp(scatt_dir, wl, interpolant_test)
    ssc_prop_interp_arr, scattering_matrix_interp_arr = c.interp_NDLinear()
    
    # output_dir = './'
    # output_name = 'test.txt'
    # output_full_name = output_dir + output_name
    # c.write_snow_ssc_prop(output_full_name, ssc_prop_interp_arr, c.scatt_angles, scattering_matrix_interp_arr)
    
    
    print("time of file write", time.time() - st)
 


