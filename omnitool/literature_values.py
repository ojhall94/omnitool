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
            'J':[0.72],     #Yuan, Liu & Xiang(2013) (as refd in Hawkins et al. 2017)
            'H':[0.46],     #Yuan, Liu & Xiang(2013) (as refd in Hawkins et al. 2017)
            'Ks':[0.30],    #Yuan, Liu & Xiang(2013) (as refd in Hawkins et al. 2017)
            'G':[2.85],     #Jordi et al. (2010) (as refd in Hawkins et al. 2017)
            'W1':[0.18],    #Yuan, Liu & Xiang(2013) (as refd in Hawkins et al. 2017)
            'W2':[0.16],    #Yuan, Liu & Xiang(2013) (as refd in Hawkins et al. 2017)
            'W3':[0.16],    #Xue et al. (2016) (as refd in Hawkins et al. 2017)
            'W4':[0.11],    #Xue et al. (2016) (as refd in Hawkins et al. 2017)
            })

'''Other defined constants'''
stefboltz = const.sigma_sb.value
Lzp = 3.0128e28 #Watts, zero point luminosity [E. E. Mamajek et al. 2015]

'''Solar parameters [!Still need literature sources!!]'''
'''Basic solar parameters'''
Rsol = const.R_sun.value #meters
Tsol = 5778 #K
Msol = 1.989e30 #Kg
gsol = 27400. #cms^2
Zsol = 0.01756 #[T. Rodrigues et al. 2017]
Lsol = 4 * np.pi * stefboltz * Rsol**2 * Tsol**4
Asol = 4 * np.pi * (Rsol*100)**2

'''Asteroseismic parameters'''
Numaxsol = 3090 #[Huber et al 2011]
Dnusol = 135.1  #[Huber et al 2011]

'''Photometric parameters'''
Mbolsol = 4.74 #[E. E. Mamajek et al. 2015]
err_bc = 0.02 #Assumed uncertainty in bolometric corrections [Huber et al. 2017]
err_av = 0.02 #Assumed uncertainty in extinction coefficients [Huber et al. 2017]

'''Apparent magnitudes for the sun from [Bohlin & Gilliland 2004]'''
Mbandsol = pd.DataFrame({'J':[3.64],
        'H':[3.32],
        'Ks':[3.28]})

'''Red Clump locations as given in [Hawkins et al. 2017]'''
hawkvals = dict({'Ks':-1.61,'J':-0.93,'H':-1.46})
hawkerr = 0.01
