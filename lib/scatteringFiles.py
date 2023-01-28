'''
Python file containing functions that read and write ice crystal scattering files.
'''
from os import path
from netCDF4 import Dataset
import re
import numpy as np
import matplotlib.pyplot as plt
import sys, time
import glob

def twoColRead(fileName):
    """
    Function that reads a two column file in the desired format
    :param fileName: Name of file to be read (str)
    :return: LH, RH, the left and right columns of data from the file (lists of floats)
    """
    f = open(fileName, "r")
    contents = f.readlines()
    contents = [i.split() for i in contents]
    LH = [float(i[0]) for i in contents]
    RH = [float(i[1]) for i in contents]
    RH_return = []
    for r in RH:
        if r == -0.0:
            RH_return.append(0.0)
        else:
            RH_return.append(r)
    return LH, RH_return

def twoColWrite(fileName,LH,RH):
    """
    Function that writes a two column file of the standard format.
    :param fileName: Name of file to be saved to (str)
    :param LH: Left hand column of file (list of floats)
    :param RH: Right hand column of file (list of floats)
    :return: None
    """
    f = open(fileName,"w")
    if (len(LH) != len(RH)):
        return "Input arrays are of different length"
    for i in range(len(LH)):
        nSpace = 3 - len(str(int(LH[i])))
        f.write(" " * nSpace + "{:.4f}".format(LH[i]) + " " + format(RH[i], '.3E') +"\n")
    return


def createScatMatr(sourceFile, Lambda = 8.6351E-01, NR = 1.3038E0, NI = 0.0E0, Reff = 4.7E2, Area = 4.48902E6, AR= 1.0, Delta = 0.7
              , Cext = 5.20041E5, Csca = 5.19059E5, ssa=0.998112, g = 0.713641, Veff = 1.0, Pizero = 0.998112, cosb = 0.713641):
    """
    Function that accepts snow parameters (aspect ratio, effective radius, etc.) and creates a scattering matrix .dat
    file in the format necessary for use in the radiative transfer code. Does not return anything, but saves file
    :param sourceFile: Name of file to create scattering matrix from, of the form ASxx.xxx_Ryyy.yy_WLzzzz_D0.00 (string)
    :param Lambda: Wavelength (float)
    :param NR: Real component of index of refraction (float)
    :param NI: Imaginary component of index of refraction (float)
    :param Reff: Effective radius of crystal (float)
    :param Area: Area of crystal (float)
    :param AR: Aspect ratio of crystal (float)
    :param Delta: Roughness of crystal (float)
    :param Cext: Extinction crosssection (float)
    :param Csca: Scattering crosssection (float)
    :param ssa: Single scattering albedo (float)
    :param g: Assymetry parameter (float)
    :param Veff: Effective Variance (float)
    :return: None
    """
    f = open(sourceFile + ".dat", "w")
    ###         Condition to check whether the sourceFile name is valid (i.e. whether .p11, .p12, etc. files exist with that name)
    # print("first Conditiom", path.exists(sourceFile + ".p11"))
    # print("second condition", path.exists(sourceFile + ".p12"))
    # print("third condition", path.exists(sourceFile + ".p22"))
    # print("fourth condition", path.exists(sourceFile + ".p33"))
    # print("fifth condition", path.exists(sourceFile + ".p34"))
    # print("sixth", path.exists(sourceFile + ".p44"))
    condition = path.exists(sourceFile + ".p11") and path.exists(sourceFile + ".p12") and path.exists(sourceFile + ".p22") \
    and path.exists(sourceFile + ".p33") and path.exists(sourceFile + ".p34") and path.exists(sourceFile + ".p44")

    ###         Checking if condition is not met. If it is not met, error message is written in output file, function
    ###         is terminated
    if not condition:
        f.write("No raytracing files (.Pnm files) are present in directory for given filename. Check filename.")
        f.close()
        return

    ###         Creating file of the name given with a .dat file type
    f.write("Ice optical properties for given values of\n")
    f.write("effective radius (Reff), aspect ratio (AR) and distortion parameter (delta)\n")
    f.write("Source File: " + sourceFile + "\n")
    # f.write("Lambda= " + format(Lambda, '.5E') + "    " + "NR= " + format(NR, '.5E') + "    " + "NI= " + format(NI, '.5E') + "\n")
    # f.write("Reff= "+format(Reff, ".5E")+"    Area= "+format(Area,".5E")+"    AR="+("{:.6f}").format(AR)+"    Delta="+("{:.2f}").format(Delta)+"\n")
    # f.write("Cext= " + format(Cext,".5E") + "    Csca= "+format(Csca,".5E")+"    ssa="+("{:.6f}").format(ssa)+"     g ="+("{:.6f}").format(g)+"\n")
    numSpaceVeff = 1 - len(str(Veff).split(".")[0])
    f.write("Reff= "+format(Reff, ".5E") + "    Veff=" + " "*numSpaceVeff + ("{:.4f}").format(Veff) + "    Area= "+format(Area,".5E") + "    Lambda= " + format(Lambda, '.5E'))
    f.write("    NR= " + format(NR, '.5E') + "    NI= " + format(NI, '.5E') + "\n")
    f.write("Cext= " + format(Cext, ".5E") + "    Csca= " + format(Csca, ".5E") + "    Pizero=" + ("{:.6f}").format(Pizero) + "    Cosbar=" + ("{:.6f}").format(cosb) + "\n")
    f.write("Angle (degs)          P11          P22/P11         P33/P11         P44/P11         P12/P11         P34/P11")
    f.write("\n")
    angles, P11 = twoColRead(sourceFile + ".p11")
    angles, P12 = twoColRead(sourceFile + ".p12")
    angles, P22 = twoColRead(sourceFile + ".p22")
    angles, P33 = twoColRead(sourceFile + ".p33")
    angles, P34 = twoColRead(sourceFile + ".p34")
    angles, P44 = twoColRead(sourceFile + ".p44")
    for i in range(len(angles)):
        f.write(" "+format(angles[i],".5E"))
        ###         Note that we must adjust the number of whitespaces if the number is negative
        ###         Thus, nspace is 5 if the number is positive, 4 if negative (to leave space for negative sign)
        nSpace = (5 if P11[i]>=0.0 else 4)
        f.write(" "*nSpace + format(P11[i],".5E"))
        nSpace = (5 if P22[i] >= 0.0 else 4)
        f.write(" "*nSpace + format(P22[i] , ".5E"))
        nSpace = (5 if P33[i] >= 0.0 else 4)
        f.write(" "*nSpace + format(P33[i] , ".5E"))
        nSpace = (5 if P44[i] >= 0.0 else 4)
        f.write(" "*nSpace + format(P44[i] , ".5E"))
        nSpace = (5 if P12[i] >= 0.0 else 4)
        f.write(" "*nSpace + format(P12[i] , ".5E"))
        nSpace = (5 if P34[i] >= 0.0 else 4)
        f.write(" "*nSpace + format(P34[i] , ".5E"))
        f.write("\n")
    f.close()
    return

def mixCrystals(filename, nRough, nAr, nReff, f_col, outFileName):
    """
    Function that "mixes" column and plate ice crystals of equivalent aspect ratios. Uses functions defined in this class
    to save appropriate files.
    :param filename: Name of file .cdf file
    :param nRough: Index of the roughness array to use (int)
    :param nAr: Index of the aspect ratio array to use, must be less than 25 (length of the array divided by 2) (int)
    :param nReff: Index of the effective radius to use (int)
    :param f_col: Relative area fraction of columns in mixture (float)
    :return: None
    """
    data = Dataset(filename, mode='r')

    ###         Defining column variables from data file
    ###         Note that in order to satisfy that nAr_columns = 1/ nAr, we use 50 - nAr for the
    C_ext_col = data["EXTINCTION_CROSSSECTION"][:][nRough][50 - nAr][nReff]
    area_col = data["AREA"][:][50 - nAr][nReff]
    C_sca_col = data["EXTINCTION_CROSSSECTION"][:][nRough][50 - nAr][nReff]
    g_col = data["ASYMMETRY_PARAMETER"][:][nRough][50 - nAr][nReff]
    P11_col = [j[nRough][50 - nAr][nReff] for j in data["P11"][:]]
    P12_col = [j[nRough][50 - nAr][nReff] for j in data["P12"][:]]
    P22_col = [j[nRough][50 - nAr][nReff] for j in data["P22"][:]]
    P33_col = [j[nRough][50 - nAr][nReff] for j in data["P33"][:]]
    P34_col = [j[nRough][50 - nAr][nReff] for j in data["P34"][:]]
    P44_col = [j[nRough][50 - nAr][nReff] for j in data["P44"][:]]

    ###         Defining plate variables from data file
    P11_pla = [j[nRough][nAr][nReff] for j in data["P11"][:]]
    P12_pla = [j[nRough][nAr][nReff] for j in data["P12"][:]]
    P22_pla = [j[nRough][nAr][nReff] for j in data["P22"][:]]
    P33_pla = [j[nRough][nAr][nReff] for j in data["P33"][:]]
    P34_pla = [j[nRough][nAr][nReff] for j in data["P34"][:]]
    P44_pla = [j[nRough][nAr][nReff] for j in data["P44"][:]]
    area_pla = data["AREA"][:][nAr][nReff]
    g_pla = data["ASYMMETRY_PARAMETER"][:][nRough][nAr][nReff]

    ###     Defining intermediate variables
    Q_ext_col = C_ext_col / area_col
    Q_sca_col = C_sca_col / area_col
    Q_ext_pla = C_ext_col / area_pla
    Q_sca_pla = C_sca_col / area_pla
    Q_ext_mix = f_col * Q_ext_col + (1 - f_col) * Q_ext_pla
    Q_sca_mix = f_col * Q_sca_col + (1 - f_col) * Q_sca_pla

    ###     Defining mixed variables
    P11_mix = [(f_col * P11_col[i] * Q_sca_col + (1 - f_col) * P11_pla[i] * Q_sca_pla) / Q_sca_mix for i in
               range(len(P11_pla))]
    P12_mix = [(f_col * P12_col[i] * Q_sca_col + (1 - f_col) * P12_pla[i] * Q_sca_pla) / Q_sca_mix for i in
               range(len(P12_pla))]
    P22_mix = [(f_col * P22_col[i] * Q_sca_col + (1 - f_col) * P22_pla[i] * Q_sca_pla) / Q_sca_mix for i in
               range(len(P22_pla))]
    P33_mix = [(f_col * P33_col[i] * Q_sca_col + (1 - f_col) * P33_pla[i] * Q_sca_pla) / Q_sca_mix for i in
               range(len(P33_pla))]
    P34_mix = [(f_col * P34_col[i] * Q_sca_col + (1 - f_col) * P34_pla[i] * Q_sca_pla) / Q_sca_mix for i in
               range(len(P34_pla))]
    P44_mix = [(f_col * P44_col[i] * Q_sca_col + (1 - f_col) * P44_pla[i] * Q_sca_pla) / Q_sca_mix for i in
               range(len(P44_pla))]
    g_mix = (f_col * g_col * Q_sca_col + (1 - f_col) * g_pla * Q_sca_pla) / Q_sca_mix
    Area_mix = area_col
    C_ext_mix = Q_ext_mix * Area_mix
    C_sca_mix = Q_sca_mix * Area_mix
    Ssa_mix = Q_sca_mix / Q_ext_mix
    ar_mix = data["ASPECT_RATIO"][:][nAr]
    scat_mix = data["SCATTERING_ANGLES"][:]
    lambda_mix = data["WAVELENGTH"][:][0] / 1000
    NR_mix = data["REAL_REFRACTIVE_INDEX"][:][0]
    NI_mix = data["IMAG_REFRACTIVE_INDEX"][:][0]
    Reff_mix = data["REAL_EFFECTIVE_RADIUS"][:][nAr][nReff]
    rough_mix = data["DISTORTION"][:][nRough]

    twoColWrite(outFileName + ".p11", scat_mix, P11_mix)
    twoColWrite(outFileName + ".p12", scat_mix, P12_mix)
    twoColWrite(outFileName + ".p22", scat_mix, P22_mix)
    twoColWrite(outFileName + ".p33", scat_mix, P33_mix)
    twoColWrite(outFileName + ".p34", scat_mix, P34_mix)
    twoColWrite(outFileName + ".p44", scat_mix, P44_mix)
    # print("FILE PASSED TO CREATESCATMATR",outFileName)
    createScatMatr(outFileName , lambda_mix, NR_mix, NI_mix, Reff_mix, Area_mix, ar_mix, rough_mix, C_ext_mix, C_sca_mix, Ssa_mix, g_mix)
    return

def calcTau(filename, rho_snow, thickness):
    """
    Function that calculates the optical depth for a given scattering matrix file and a given layer thickness
    :param filename: Filename of scattering matrix file
    :param thickness: Thickness of snow layer (m)
    :return: Optical Depth of given layer
    """
    # rho_snow = 0.100 #snow density(g / cm3)
    thk_snow = thickness * 100 #thickness of the snow layer converted to cm (cm)
    rho_ice = 0.917 #ice density(g / cm3)
    f = open(filename,"r")
    # print(re.split("AS|_", filename))
    # AR = float(re.split("AS|_", filename)[3])
    lines = f.readlines()
    nums = re.split("Reff=|Veff=|Area=|Lambda=|NR=|NI=|\n",lines[3])
    nums2 = re.split("Cext=|Csca=|Pizero=|Cosbar=|\n", lines[4])
    vars = [float(nums[i]) for i in range(len(nums)) if (nums[i]!='' and nums[i]!=" ")]
    vars2 = [float(nums2[i]) for i in range(len(nums2)) if (nums2[i]!='' and nums2[i]!=" ")]
    Reff, Veff, Area, Lambda, NR, NI = vars
    Cext, Csca, Pizer, Cosbar = vars2
    K_e = 3. / (2.*rho_ice*Reff*(1e-4)) # MASS EXTINCTION COEFFICIENT
    Kext = rho_snow * K_e
    tau = Kext*thk_snow
    # print("CEXT:", Kext)
    # print("thicknessL",thickness)
    return tau

class getIOPs:
    """
    class for inherent optical properties of scattering particles
    
    """
    
    def calcTau_imp(self, f_imp, rho_snow, rho_imp, wl, nr_imp, ni_imp, thickness):
        """
        Calculates the optical depth of impurity for a slab with a given thickness 
        
        optical depth τ is given by
          τ = N × sigma_ext × 10^−8 × z
          , where N = f/v x rho_snow/rho_imp x 10^12, 
        sigma_ext = extinction cross-section of the impurity particles in micron^2
        z = thickness of the slab in cm 
        f = impurity mass concentration in ppmw
        v = average volume of impurity particles
        rho_snow = snow water equivalent in g/cm^3
        rho_imp = impurity particle mass density in g/cm^3

        Parameters
        ----------
        f_imp : float
            mass fraction of snow impurity in ppmw (parts per million weight).
        rho_snow : float
            snow water equivalent of the snow pack in g/cm^3.
        rho_imp : float
            impurity particle mass density in g/cm^3.
        wl : float
            wavelength of the calculation.
        nr_imp : float
            real part of the impurity particle refractive index.
        ni_imp : float
            imaginary part of the impurity particle refractive index.
        thickness : float
            thickness of the slab in meters.

        Returns
        -------
        float
            optical depth of snow impurities within the slab.

        """
        
        ## unit conversions
        f_imp     = f_imp * 1.0e-6 # convert from ppmw to ppw
        thickness = thickness * 100. # convert to cm
    
        ## obtain the single scattering properties through mie calculations
        self = self.mie_calc(self, wl, nr_imp, ni_imp)
        
        ## calculate N 
        N = f_imp/self.volume * rho_snow/rho_imp *1.e12
        # print(N)
        
        ## calculate the optical depth
        tau = N * self.cext * 1.0e-8 * thickness
        # print(tau)
        
        return tau 
    
    def mie_calc(self, wl, nr, ni):
        """
        simple driver of the Mie code to calculate the ensemble average of 
        extinction crosssection, scattering crossection, single scattering albedo and 
        asymmetry parameter

        Parameters
        ----------
        wl : float
            wavelength for the calculation in micron.
        nr : float
            real part of the particle refractive index.
        ni : float
            imaginary part of the particle refractive index.

        Returns
        -------
        class
            the particle class with scattering properties.

        """

        import miepython
        
        ## Obtain the size distribution function and weights
        r_arr = self.x
        pdf_arr = self.y 
        w_arr = self.w 
        
        qext_arr = []
        qsca_arr = []
        qback_arr = []
        g_arr = []
        
        for r in r_arr:
            
           x = 2.0*np.pi*r/wl # size parameter
           m = complex(nr, -ni)
           qext, qsca, qback, g = miepython.mie(m, x)
           # print(qext, qsca, qback, g)
           qext_arr.append(qext)
           qsca_arr.append(qsca)
           qback_arr.append(qback)
           g_arr.append(g)
        
        qext_arr = np.array(qext_arr)
        qsca_arr = np.array(qsca_arr)  
        qback_arr = np.array(qback_arr)  
        g_arr = np.array(g_arr)  

        qext_mean = sum(qext_arr*pdf_arr*w_arr)
        qsca_mean = sum(qsca_arr*pdf_arr*w_arr)
        qback_mean = sum(qback_arr*pdf_arr*w_arr)
        g_mean = sum(g_arr*qsca_arr*r_arr*r_arr*pdf_arr*w_arr)*np.pi
        
        cext_mean = sum(qext_arr*r_arr*r_arr*pdf_arr*w_arr)*np.pi
        csca_mean = sum(qsca_arr*r_arr*r_arr*pdf_arr*w_arr)*np.pi
        g_mean = g_mean/csca_mean
        ssa_mean = csca_mean/cext_mean
        # print(cext_mean, ssa_mean, g_mean )        
        self.cext = cext_mean
        self.csca = csca_mean
        self.ssa  = ssa_mean
        self.g    = g_mean       
        
        return self
    
    @classmethod
    def sizedis_spher(cls, aa, bb, dist='lognormal', r1=1.0e-5, r2=1.0e6, \
                      resetr=True, cutoff=1.0e-8, n=1000):
        """
        Returns the size distribution for a collection of spherical particles

        Parameters
        ----------
        cls : class
            the particle class.
        aa : float
            distribution parameter 1.
        bb : float
            distribution parameter 2.
        dist : str, optional
            distribution type. The default is 'lognormal'.
        r1 : float, optional
            lower limit of integration. The default is 1.0e-5.
        r2 : float, optional
            upper limit of integration. The default is 1.0e6.
        resetr : bool, optional
            switch to reset the upper and lower limit. The default is True.
        cutoff : float, optional
            cut off threshold for upper and lower limit reset. The default is 1.0e-8.
        n : int, optional
            number of quadratures. The default is 1000.

        Returns
        -------
        class
            the particle class.

        """
        
        if dist=='lognormal': # log-normal distribution
        
            mu    = aa
            sigma = bb
            
            ## automatically reset the upper and lower limit of integration if resetr = True
            if resetr:
                sigma_sq = sigma * sigma
                ymax = 1./(np.sqrt(2.0*np.pi*sigma_sq)*mu*np.exp(-1.0*sigma_sq*0.5))
                ycut = ymax * cutoff
                r1_reset = np.exp( (np.log(mu)- sigma_sq) - \
                                  np.sqrt((np.log(mu)- sigma_sq)*(np.log(mu)- sigma_sq) - \
                                          np.log(mu)*np.log(mu) - sigma_sq*np.log(2.0*np.pi*sigma_sq) - \
                                          2.0*sigma_sq*np.log(ycut) ) )
                r2_reset = np.exp( (np.log(mu)- sigma_sq) + \
                                  np.sqrt((np.log(mu)- sigma_sq)*(np.log(mu)- sigma_sq) - \
                                          np.log(mu)*np.log(mu) - sigma_sq*np.log(2.0*np.pi*sigma_sq) - \
                                          2.0*sigma_sq*np.log(ycut) ) )
                    
                x, w = getIOPs.gauleg(r1_reset, r2_reset, n)    
                y = getIOPs.lognormal_pdf(x, mu, sigma)
                
            else:
                x, w = getIOPs.gauleg(r1, r2, n)
                y = getIOPs.lognormal_pdf(x, mu, sigma)            
                        
        else:
            
            print('Size distribution ', dist, ' not implemented, exit...')
            sys.exit()
        
        ## renormalize the distribution function
        sum_y = sum(y*w)
        y = y/sum_y
        
        ## now integrate over size distribution to get reff, veff, avg volume
        # second and third mode
        mode1 = sum(x*y*w)
        mode2 = sum(x*x*y*w)
        mode3 = sum(x*x*x*y*w)
        # effective radius
        reff = mode3/mode2
        # effective variance 
        xi = x-reff
        veff = sum(xi*xi*x*x*y*w)/(mode2*reff*reff)
        # mean area
        area = mode2*np.pi
        # mean volume
        volume = mode3*4.0*np.pi/3.0   
        # mean radius
        rmean  = mode1         

        cls.y = y
        cls.x = x
        cls.w = w        
        cls.reff = reff
        cls.veff = veff
        cls.area = area
        cls.volume = volume
        cls.rmean = rmean

        return cls
    
    @staticmethod
    def lognormal_pdf(x, mu, sigma):
        """
        Returns a log-normally distributed probability density function (pdf)
        
        Parameters
        ----------
        x : float array
            input array for x.
        mu : float
            mode radius in linear space
        sigma : float
            geometric standard deviation
    
        Returns
        -------
        y : float array
           pdf as a function of x
    
        """
        
        y1 = np.log(x) - np.log(mu)
        y = np.exp(-y1*y1*0.5/(sigma*sigma))/x
        
        return y
    
    @staticmethod
    def gauleg(x1, x2, n):
        """
        Returns the gaussian quadrature points and weights 
        for integration between x1 and x2
    
        Parameters
        ----------
        x1 : float
            Lower limit of x
        x2 : float
            Upper limit of x
        n : int
            Number of quadrature points 
    
        Returns
        -------
        x_out : float array
            Double Gaussian quadrature points 
        w_out : float array
            Weights
    
        """
        from scipy.special.orthogonal import p_roots
    
        # find the roots for Legendre Polynomials between [-1, 1]
        [x, w] = p_roots(n)
        x = np.flip(x)
    
        xm=0.5*(x2+x1)
        xl=0.5*(x2-x1)
        
        x_out = xm - xl*x # scale the roots to the desired range
        w_out = xl*w      # scale the weights
        
        return x_out, w_out    
    
if __name__ == '__main__':
    dens_t = 0.315
    dep_t = 0.00001
    AS_t = 0.457
    reff_t = 56.0
    tau = calcTau("/home/accurt/Downloads/ice_reff_555/" + f"AS" + ("%06.3f" % AS_t) + "_REFF" + (
                        "%08.2f" % reff_t) + "_WL0555_D0.05_K0.75.dat", dens_t, dep_t)
    print(tau)
    quit()

    ang, p34 = twoColRead("/home/accurt/Downloads/ice_reff_670/AS22.865_REFF01280.00_WL0670_D0.35_K0.75.p34")
    print(p34)
    print(p34[0] == -0.0)
    quit()
    allAR = []
    allReff = []
    allD = []
    createScatMatr(sourceFile="/home/accurt/Downloads/raytrace1589/AS00.020_R005.00_WL1589_D0.00",
                Lambda = 15.89E-2, NR = 1.3038E0, NI = 0.0E0, Reff = 4.7E2, Area = 4.48902E6, AR= 1.0, Delta = 0.0,
                Cext = 5.20041E5, Csca = 5.19059E5, ssa=0.998112, g = 0.713641, Veff = 1.0, Pizero = 0.998112,
                cosb = 0.713641)

    fList = glob.glob(r"/home/accurt/Downloads/raytrace1589/*.1dr")
    fList = glob.glob(r"/home/accurt/Downloads/raytrace1589/*.1dr")
    print(fList[:5])
    print(re.split("AS|_R|_WL|_D", "/home/accurt/Downloads/raytrace1589/AS00.020_R005.00_WL1589_D0.00"))
    AR  = float(re.split("AS|_R|_WL|_D", "/home/accurt/Downloads/raytrace1589/AS00.020_R005.00_WL1589_D0.00")[1])
    Reff= float(re.split("AS|_R|_WL|_D", "/home/accurt/Downloads/raytrace1589/AS00.020_R005.00_WL1589_D0.00")[2])
    WL  = float(re.split("AS|_R|_WL|_D", "/home/accurt/Downloads/raytrace1589/AS00.020_R005.00_WL1589_D0.00")[3])
    D   = float(re.split("AS|_R|_WL|_D", "/home/accurt/Downloads/raytrace1589/AS00.020_R005.00_WL1589_D0.00")[4])
    WL = np.format_float_scientific(WL, precision = 4, exp_digits = 2)
    print(AR, Reff, WL, D)
    for f in fList:
        # f.removesuffix(".dr")
        AR = float(re.split("AS|_R|_WL|_D|\.1dr", f)[1])
        Reff = float(re.split("AS|_R|_WL|_D|\.1dr", f)[2])
        WL = float(re.split("AS|_R|_WL|_D|\.1dr", f)[3])
        D = float(re.split("AS|_R|_WL|_D|\.1dr", f)[4])
        # WL = np.format_float_scientific(WL, precision=4, exp_digits=2) AS22.865_R056.00_WL1589_D0.85.dat
        allAR.append(AR), allReff.append(Reff), allD.append(D)
        print(AR, Reff, WL, D)
        print(type(WL))
        print(f[:-4])
        createScatMatr(sourceFile=f[:-4],
                       Lambda=WL, NR=1.3038E0, NI=0.0E0, Reff=Reff, Area=4.48902E6, AR=AR, Delta=D,
                       Cext=5.20041E5, Csca=5.19059E5, ssa=0.998112, g=0.713641, Veff=1.0, Pizero=0.998112,
                       cosb=0.713641)

    start_time = time.time()
    allAR  = list(set(allAR))
    allReff = list(set(allReff))
    allD = list(set(allD))
    allAR.sort(), allReff.sort(), allD.sort()
    print("_____AR____", allAR)
    print("______Reff____", allReff)
    print("_______D_______", allD)
    fdatList = glob.glob(r"/home/accurt/Downloads/raytrace1589/*.dat")
    print("_______EQUAL________", len(fList) == len(fdatList))
    # x1 = 0.0
    # x2 = 1.0e3
    # n = 1000
    # x, w = getIOPs.gauleg(x1, x2, n)    
    # # G = sum(w * x)  # G is the integral of y = x
    # # G = sum(w) # G is the integral of y = 1

    # mu = 0.046
    # sigma = 1.5
    # y = getIOPs.lognormal_pdf(x, mu, sigma)
    
    # fig, ax = plt.subplots()
    # ax.plot(x, y)
    
    ## define the impurity particle size distribution
    mu = 0.046
    sigma = 1.5
    imp = getIOPs.sizedis_spher(mu, sigma)
    
    # reff =    imp.reff
    # veff =    imp.veff 
    # area =    imp.area 
    # volume =  imp.volume   
    # rmean = imp.rmean
    
    ## calculate the optical depth of the impurity for a slab of snow 
    f_imp = 1.0  # ppmw
    rho_snow = 1.0 # g/cm^3
    rho_imp  = 1.8 # g/cm^3
    wl       = 0.86351
    nr_imp   = 1.85 # 
    ni_imp   = 0.71 #
    thickness = 0.98 # meters
    
    tau = getIOPs.calcTau_imp(imp, f_imp, rho_snow, rho_imp, wl, nr_imp, ni_imp, thickness)
    print(tau)
    
    ## perform mie calculation separately 
    # imp = getIOPs.mie_calc(imp, wl, nr_imp, ni_imp)
    # print(imp.g)
    
    print('calculation finished in %s seconds' % (time.time() - start_time) )
    
    
## below is the size distribution routine to be converted to Python     
# subroutine SIZDIS_MIE(NDISTR,AA,BB,GAM,R1,R2,YY,WY,NNK,NP,NK,N,NNPK, &
#                       REFF,VEFF,RVW,RMEAN,VOLUME,AREA)
# implicit none
# integer, intent(in) :: NDISTR,NNK,NP,NK,N,NNPK
# real(kind=8),intent(in) :: AA,BB,GAM,R1,R2
# real(kind=8),dimension(NNK),intent(inout) :: YY,WY
# real(kind=8),dimension(NNK) :: YTMP
# real(kind=8) :: SUM,XD,YD,DA,PI,TEMP,R1TMP,R2TMP,B2,DAB,A2,DB
# real(kind=8) :: G,REFF,VEFF,RVW,RMEAN,XI,VOLUME,AREA
# integer      :: I,J
# real(kind=8),dimension(NK)  :: X,W
# real(kind=8),dimension(NNK) :: XX
# real(kind=8) :: Z1,Z2,ZJ
# integer :: J1,J2,IJ,I1

# PI = DACOS(-1.0D0)
# CALL GAUSS (NK,0,0,X,W)
#       IF (NDISTR.NE.5) GO TO 30
#       Z1=R1/DFLOAT(NP)                                                      
#       Z2=Z1*0.5D0                                                               
#       DO 28 I=1,NK                                                              
#          XX(I)=Z2*X(I)+Z2                     
#    28 CONTINUE                                                                  
#       DO 29 J=1,NP                                                              
#          J1=J-1                                                                 
#          ZJ=Z1*DFLOAT(J1)                                                    
#          IJ=J1*NK                                                               
#          DO 29 I=1,NK                                                           
#             I1=I+IJ                                                             
#             YY(I1)=XX(I)+ZJ
#             WY(I1)=W(I)*Z2                                                      
#    29 CONTINUE                                                                  
#    30 Z1=(R2-R1)/DFLOAT(N)                                                      
#       Z2=Z1*0.5D0                                                               
#       DO 32 I=1,NK                                                              
#          XX(I)=Z2*X(I)+Z2                        
#    32 CONTINUE                                                                  
#       DO 34 J=1,N                                                               
#          J1=J-1                                                                 
#          ZJ=Z1*DFLOAT(J1)+R1                                                    
#          IJ=J1*NK+NNPK                                                         
#          DO 34 I=1,NK                                                           
#             I1=I+IJ                                                             
#             YY(I1)=XX(I)+ZJ  
#             WY(I1)=W(I)*Z2                                                      
#    34 CONTINUE 

# SELECT CASE (NDISTR)
# CASE(1)
# !     PRINT 1000,AA,BB,GAM 
# 1000 FORMAT('MODIFIED GAMMA DISTRIBUTION, ALPHA=',F6.4,'  r_c=',   &            
#        F6.4,'  GAMMA=',F6.4)                                                   
#  A2=AA/GAM                                                                 
#  DB=1D0/BB
#  DO I=1,NNK                                                             
#   XD=YY(I)                                                             
#   YD=XD**AA                                                                
#   XD=XD*DB
#   YD=YD*DEXP(-A2*(XD**GAM))
#   WY(I)=WY(I)*YD
#  ENDDO
 
# CASE(2)
# ! PRINT 2000,AA,BB,R1,R2                                                               
#  2000 FORMAT(' LOG NORMAL DISTRIBUTION',/, &            
#              ' r_g=',F9.4,', [ln(sigma_g)]**2=',F9.4,', R1=',F9.4,', R2=',F9.4) 
# ! CALL GAULEG(R1,R2,YY,W,NNK)
#  !R1TMP = DLOG(R1)
#  !R2TMP = DLOG(R2)
#  !CALL GAULEG(R1TMP,R2TMP,YTMP,W,NNK)
#  !DO I = 1, NNK
#  !YY(I) = EXP(YTMP(I))
#  !ENDDO
#  DA=1.0D0/AA
#  DO I = 1, NNK
#  XD=YY(I)
#  YD=DLOG(XD*DA)
#  YD=DEXP(-YD*YD*0.5D0/BB)/XD
#  WY(I)=WY(I)*YD 
#  !PRINT*, 'X = ',X,', Y = ', Y, ', WY(I)=',WY(I) 
#  ENDDO
 
# CASE(3)
#  CALL POWER (AA,BB,R1,R2) 
# ! CALL GAULEG(R1,R2,YY,W,NNK)
# ! PRINT 3000,AA,BB,R1,R2                                                                
#  3000 FORMAT(' POWER LAW DISTRIBUTION OF HANSEN & TRAVIS (1974)',/,&
#              ' REFF=',F9.4,', VEFF=',F9.4,', R1=',F9.4,', R2=',F9.4)              
#  DO I = 1, NNK                                                            
#   XD=YY(I)                                                                
#   WY(I)=WY(I)/(XD*XD*XD)                                                 
#  ENDDO      
          
# CASE(4)
# ! PRINT 4000,AA,BB     
#  4000 FORMAT ('GAMMA DISTRIBUTION,  a=',F8.4,'  b=',F6.4)
#  B2=(1D0-3D0*BB)/BB                                                        
#  DAB=1D0/(AA*BB)                                                          
#  DO I=1,NNK                                                            
#   XD=YY(I)                                                                
#   XD=(XD**B2)*DEXP(-XD*DAB)                                                 
#   WY(I)=WY(I)*XD
#   !YI(I) = XD 
#  ENDDO  
                                                  
# CASE(5)
# ! PRINT 5000,BB                                                             
#  5000 FORMAT ('MODIFIED POWER LAW DISTRIBUTION,  ALPHA=',D10.4)
#  DO I=1,NNK                                                            
#   XD=YY(I)                                                                
#   IF (XD.LE.R1) WY(I)=WY(I)
#   !IF (XD.LE.R1) YI(I)=1.0D0 
#   IF (XD.GT.R1) WY(I)=WY(I)*(XD/R1)**BB
#   !IF (XD.GT.R1) YI(I)=(XD/R1)**BB 
#  ENDDO 
 
# END SELECT

#  SUM=0D0
#  DO I = 1, NNK
#   SUM=SUM+WY(I)
#  ENDDO
#  SUM=1D0/SUM
#  DO I = 1, NNK
#   WY(I)=WY(I)*SUM
#  ENDDO   
 
#  G=0D0
#  DO I=1,NNK
#     XD=YY(I)
#     G=G+XD*XD*WY(I)
#  ENDDO
 
#  REFF=0D0
#  DO I=1,NNK
#     XD=YY(I)
#     REFF=REFF+XD*XD*XD*WY(I)
#  ENDDO
 
#  REFF=REFF/G
#  VEFF=0D0
#  VOLUME=0D0
#  RVW=0D0
#  RMEAN=0D0
 
#  DO I=1,NNK
#     XD=YY(I)
#     XI=XD-REFF
#     VEFF=VEFF+XI*XI*XD*XD*WY(I)
#     VOLUME=VOLUME+XD*XD*XD*WY(I)
#     RVW=RVW+XD*XD*XD*XD*WY(I)
#     RMEAN=RMEAN+XD*WY(I)
#  ENDDO
 
#  VEFF=VEFF/(G*REFF*REFF)
#  AREA=G*PI  
#  RVW=RVW/VOLUME
#  VOLUME=VOLUME*4D0*PI/3D0

# return  
# end subroutine SIZDIS_MIE        
