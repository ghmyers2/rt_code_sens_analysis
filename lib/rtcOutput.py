'''
rtcOutput class defined for use with plotting
Instances of the class are initialized with the name of the file to be unpacked
fileHandling function is used to unpack variables, and class variables are initialized based on the output
of that function
'''

from scipy.io import FortranFile
import numpy as np


def readRSP(fileName):
    """Define Variables of Interest from a given fortran file.
        Keyword arguments:
        fileName -- Name of fortran file for unpacking (string)
        [MAXVIEW,MAXLAYER,MAXKERN,NVIEW,NLAYER,NKERN,IINT,ISRF_PERT,PHI0,XMU0,THETAV,C22RA,S22RA,RV11,RV21,RV31] -- All
        desired variables packed into a single list (list)
    """
    f = FortranFile(fileName, 'r')

    #       CALLING THE read_ints AND read_reals FUNCTIONS IS A WAY OF READING THE FORTRAN FILE LINE BY LINE
    #       WE THEN UNPACK EACH LINE BASED ON THE VARIABLES CONTAINED IN THAT LINE (AS OUTLINED IN README)

    firstLine = f.read_ints(np.int32)  # CONTAINS MAXVIEW,MAXLAYER,MAXKERN,NVIEW,NLAYER,NKERN,IINT,ISRF_PERT
    secondLine = f.read_reals(np.float64)  # CONTAINS PHI0,XMU0
    thirdLine = f.read_reals(np.float64)  # CONTAINS THETAV,C22RA,S22RA
    fourthLine = f.read_reals(np.float64)  # CONTAINS RV11,RV21,RV31
    fifthLine = f.read_reals(np.float64)
    sixthLine = f.read_reals(np.float64)
    seventhLine = f.read_reals(np.float64)
    eigthLine = f.read_reals(np.float64)

    # Note that we must use extra array slicing when defining THETAV, C22RA, S22RA, RV11, RV21, and RV31, as
    # these variables have trailing zeros at the end of them due to a difference in the true number of angles
    # used and the space allotted to the array based on MAXVIEW.

    MAXVIEW, MAXLAYER, MAXKERN, NVIEW, NLAYER, NKERN, IINT, ISRF_PERT = firstLine
    PHI0, XMU0 = secondLine
    THETAV = np.trim_zeros(thirdLine[:MAXVIEW], "b")
    C22RA = thirdLine[MAXVIEW:2 * MAXVIEW]
    C22RA = C22RA[:len(THETAV)]
    S22RA = thirdLine[2 * MAXVIEW:3 * MAXVIEW]
    S22RA = S22RA[:len(THETAV)]
    RV11 = fourthLine[:MAXVIEW]
    RV11 = RV11[:len(THETAV)]
    RV21 = fourthLine[MAXVIEW:2 * MAXVIEW]
    RV21 = RV21[:len(THETAV)]
    RV31 = fourthLine[2 * MAXVIEW:3 * MAXVIEW]
    RV31 = RV31[:len(THETAV)]
    fif = fifthLine
    firstSec = fif[:MAXVIEW*MAXLAYER]
    secondSec = fif[MAXVIEW*MAXLAYER:2*MAXVIEW*MAXLAYER]
    thirdSec = fif[2*MAXVIEW*MAXLAYER:3*MAXVIEW*MAXLAYER]
    UV11 = [firstSec[i:i + MAXVIEW][:len(THETAV)] for i in range(0, len(firstSec), MAXVIEW)]
    UV21 = [secondSec[i:i + MAXVIEW][:len(THETAV)] for i in range(0, len(secondSec), MAXVIEW)]
    UV31 = [thirdSec[i:i + MAXVIEW][:len(THETAV)] for i in range(0, len(thirdSec), MAXVIEW)]
    RZ11 = eigthLine[:MAXVIEW]
    RZ11 = RZ11[:len(THETAV)]
    RZ21 = eigthLine[MAXVIEW:2 * MAXVIEW]
    RZ21 = RZ21[:len(THETAV)]
    RZ31 = eigthLine[2 * MAXVIEW:3 * MAXVIEW]
    RZ31 = RZ31[:len(THETAV)]
    return MAXVIEW, MAXLAYER, MAXKERN, NVIEW, NLAYER, NKERN, IINT, ISRF_PERT, PHI0, XMU0, THETAV, C22RA, S22RA, RV11, RV21, RV31, UV11, UV21, UV31, RZ11, RZ21, RZ31

class rtcOutput:
    def __init__(self, fileName):
        self.fileName = fileName
        self.MAXVIEW, self.MAXLAYER, self.MAXKERN, self.NVIEW, self.NLAYER, \
        self.NKERN, self.IINT, self.ISRF_PERT, self.PHI0, self.XMU0, \
        self.THETAV, self.C22RA, self.S22RA, self.RV11, self.RV21,\
        self.RV31, self.UV11, self.UV21, self.UV31, self.RZ11, self.RZ21, self.RZ31 = readRSP(fileName) 

        

