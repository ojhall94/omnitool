#!/bin/env python
# -*- coding: utf-8 -*-

"""
A simple Python script that holds constants and other coefficients along with
their sources

.. versioncreated:: 0.1
.. codeauthor:: Oliver James Hall <ojh251@student.bham.ac.uk>
"""

import numpy as np
import pandas as pd
from astropy import constants as const

'''Extinction coefficient transformation values from [Green et al. 2018] and others'''
Av_coeffs = pd.DataFrame({'g':3.384, #Green 2018
            'r':[2.483],     #Green 2018
            'i':[1.838],     #Green 2018
            'z':[1.414],     #Green 2018
            'y':[1.126],     #Green 2018
            'J':[0.650],     #Green 2018
            'H':[0.327],     #Green 2018
            'K':[0.161],    #Green 2018
            'G':[2.85],     #Jordi et al. (2010) (as refd in Hawkins et al. 2017)
            'W1':[0.18],    #Yuan, Liu & Xiang(2013) (as refd in Hawkins et al. 2017)
            'W2':[0.16],    #Yuan, Liu & Xiang(2013) (as refd in Hawkins et al. 2017)
            'W3':[0.16],    #Xue et al. (2016) (as refd in Hawkins et al. 2017)
            'W4':[0.11],    #Xue et al. (2016) (as refd in Hawkins et al. 2017)
            })

'''Other defined constants'''
stefboltz = const.sigma_sb.value
Lzp = 3.0128e28 #Watts, zero point luminosity [E. E. Mamajek et al. 2015]f

'''Solar parameters [!Still need literature sources!!]'''
'''Basic solar parameters'''
Rsol = const.R_sun.value #meters
Tsol = 5777 #K [Williams 2013]
Msol = const.M_sun.value #Kg
gsol = 100*(const.G.value * const.M_sun.value)/(const.R_sun.value)**2 #cms^2
Zsol = 0.01756 #[T. Rodrigues et al. 2017]
Lsol = const.L_sun.value
Asol = 4 * np.pi * (Rsol*100)**2

'''Asteroseismic parameters'''
Numaxsol = 3090 #[Huber et al 2011]
eNumaxsol = 30
Dnusol = 135.1  #[Huber et al 2011]
eDnusol = 0.1

'''Photometric parameters'''
Mbolsol = 4.75 #[Casagrande & Vandenberg 2014, 2018a, b]
err_bc = 0.02 #Assumed uncertainty in bolometric corrections [Huber et al. 2017]
err_av = 0.02 #Assumed uncertainty in extinction coefficients [Huber et al. 2017]

'''Apparent magnitudes for the sun from various sources. See:
http://mips.as.arizona.edu/~cnaw/sun.html'''

Mbandsol = pd.DataFrame({'J':[3.67],    #Cohen, Wheaton & Megeath 2003
                        'H':[3.32],     #Cohen, Wheaton & Megeath 2003
                        'Ks':[3.27],    #Cohen, Wheaton & Megeath 2003
                        'W1':[3.26],    #Jarrett et al. 2011
                        'W2':[2.38],    #Jarrett et al. 2011
                        'W3':[3.26],    #Jarrett et al. 2011
                        'W4':[3.27],    #Jarrett et al. 2011
                        'B' :[5.44],    #Mann & von Braun, 2015
                        'V' :[4.81],    #Mann & von Braun, 2015
                        'R' :[4.43],    #Mann & von Braun, 2015
                        'I' :[4.10],    #Mann & von Braun, 2015
                        'Gaia':[4.68],   #Andrae et al. 2018
#Calculated using apparent magnitudes from http://mips.as.arizona.edu/~cnaw/sun.html & https://keplergo.arc.nasa.gov/CalibrationZeropoint.shtml
                        'Kepler':[4.82]
                        })

'''Bolometric correction information from Gaia papers'''
BCgsol = 0.060

'''Red Clump locations as given in [Hawkins et al. 2017]'''
hawkvals = dict({'K':-1.61,'J':-0.93,'H':-1.46, 'GAIA':0.44})
hawkerr = 0.01
