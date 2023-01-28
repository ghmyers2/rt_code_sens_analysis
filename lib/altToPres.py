import numpy as np
import matplotlib.pyplot as plt

def alt_to_geopotential_height(alt):
    """
            Returns the geopotential height at a given geometric altitude (distance from the earth's surface).
            This function is necessary for the use in the barometric formula, as the reference variables are provided
            in terms of geopotential height. This function uses equation 18 in the below reference:
                    https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf
            :param alt: Altitude in meters (float)
            :return: Geopotential height in meters' (float)
    """

    EARTH_RADIUS = 6.356766e6
    H = (EARTH_RADIUS * alt) / (EARTH_RADIUS + alt)
    return H


def getAtmRefsFromAlt(alt):
    """
        Returns the appropriate reference height, temperature lapse rate, pressure, and temperature
        at a given altitude. Values of reference variables in different atmospheric regimes are given by table 4
        and table 1 of the appendix in the below reference:
                https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf
        :param alt: Altitude in meters (float)
        :return: Reference Altitude in meters (float), Reference Temperature Lapse Rate in Kelvin / meters (float),
         Reference Pressure in milibar (float), Reference Temperature in Kelvin (float)
    """
    refAlts = [0, 11000, 20000, 32000, 47000, 51000, 71000, 84852.04584490575] # m
    refLapseRates = [-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002, 0.0] # K / m
    refPressures = [1013.25, 226.32, 54.748, 8.6801, 1.1090, 0.66938, 0.039564, 0.0037338] # mb
    refTemps =     [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.87] # K

    # Index of the reference values for the altitude argument
    refAltIndex = [i for i in range(len(refAlts)) if refAlts[i] <= alt][-1]

    return refAlts[refAltIndex], refLapseRates[refAltIndex], refPressures[refAltIndex], refTemps[refAltIndex]


def barometric(alt):
    """
            Returns the pressure in milibars at a given altitude using equations 33a and 33b in the
            below reference:
                    https://ntrs.nasa.gov/api/citations/19770009539/downloads/19770009539.pdf
            :param alt: Altitude in meters (float)
            :return: pressure: Pressure in mb (float)
    """
    g = 9.80665 # m/s^2
    M = 0.0289644 # kg/mol
    R = 8.3144598 # J/(mol·K)

    if alt <= 86000:
        hb, Lb, Pb, Tb = getAtmRefsFromAlt(alt)
        alt = alt_to_geopotential_height(alt)
        if Lb != 0:
            exp = ((-g * M) / (R * Lb))
            inner = ((Tb + Lb * (alt - hb)) / (Tb))
            pressure = Pb * (inner ** exp)
        if Lb == 0:
            inner = (-g*(M*(alt - hb)) / (R*Tb))
            pressure = Pb * np.exp(inner)

    elif 86000 < alt <= 91292.53270347098:
        pressure = 0.0030269737500209838

    elif 91292.53270347098 < alt <= 96441.28669132637:
        pressure = 1.610E-03

    elif 96441.28669132637 < alt <= 101598.26977707013:
        pressure = 6.060E-04

    elif 101598.26977707013 < alt <= 106763.50170495825:
        pressure = 2.480E-04

    elif 106763.50170495825 < alt <= 111937.00228246104:
        pressure = 1.130E-04

    elif 111937.00228246104 < alt <= 117118.79138051634:
        pressure = 6.000E-05

    elif 117118.79138051634 < alt <= 122308.88893378395:
        pressure = 3.540E-05

    else:
        pressure = 2.260E-05

    return pressure


def altToPressure(layer_m):
    """
    Conversion function of list of altitudes to list of pressures
    :param layer_m: List of altitudes in meters (list of floats)
    :return: List of pressures in mb (list of floats)
    """
    retPres = []
    for i in range(len(layer_m)):
        if i < len(layer_m) - 1:
            diff = barometric(layer_m[i]) - barometric(layer_m[i + 1])

        else:
            diff = barometric(layer_m[i]) - 0.
        retPres.append(diff)

    return retPres



if __name__ == "__main__":
    alt_list = [0.0, 0.98, 1.0, 2000., 4000., 6000., 10000., 100000.]
    dep_t = 0.0005
    alt_list = [0.0, (1.0 - dep_t), 1, 2000., 4000., 6000., 10000., 100000.]
    print(altToPressure(alt_list))
    pressure_list = [barometric(alt) for alt in alt_list]
    # print(sum(pressure_list))
    # print(pressure_list)
    print("Total pressure: " ,sum(altToPressure(alt_list)))
    plt.plot(alt_list, pressure_list)
    plt.ylabel("Pressure (mb)")
    plt.xlabel("Altitude (m)")
    # plt.show()
