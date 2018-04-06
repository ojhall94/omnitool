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

'''Extinction coefficient transformation values from Green et al. 2018'''
Av_coeffs = pd.DataFrame({'g':3.384,
            'r':[2.483],
            'i':[1.838],
            'z':[1.414],
            'y':[1.126],
            'J':[0.650],
            'H':[0.327],
            'Ks':[0.161]})

'''Other defined constants'''
stefboltz = 5.670367e-8 #Wm-2K-4
Lzp = 3.0128e28 #Watts, zero point luminosity

'''Solar parameters [!Still need literature sources!!]'''
'''Basic solar parameters'''
Rsol = 695700e3 #meters
Tsol = 5778 #K
Msol = 1.989e30 #Kg
gsol = 27400. #cms^2
Zsol = 0.01756
Lsol = 4 * np.pi * stefboltz * Rsol**2 * Tsol**4
Asol = 4 * np.pi * (Rsol*100)**2

'''Asteroseismic parameters'''
Numaxsol = 3090 #Huber et al 2011ish
Dnusol = 135.1

'''Photometric parameters'''
Mbolsol = 4.74  #Torres 2010?
Mabssol_err = 0.09

'''Apparent magnitudes for the sun from Bohlin & Gilliland 2004'''
Mbandsol = pd.DataFrame({'J':3.64,
        'H',3.32,
        'Ks',3.28})
