#!/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

from astropy.coordinates import SkyCoord
import astropy.units as units
from dustmaps.bayestar import BayestarWebQuery

import sys

from .literature_values import *

class spyglass:
    '''
    An object class that stores and calls all astrophysical components
    necessary for finding a targets absolute magnitude in a given band.

    Input:
        _ID: List of star IDs for identification

    Returns:
        M (pandas dataframe): Absolute magnitude of a target
        M_err (pandas dataframe): Error on above

    .. codeauthor:: Oliver James Hall
    '''
    def __init__(self):
        self.oo = 0.
        self.r = None
        self.m = 0.
        self.ra = 9999.
        self.dec = 9999.
        self.band = None
        self.frame = None
        self.bailerjones = False

    def check_contents(self):
        kill = False
        try:
            if self.ra == 9999.:
                print('Please pass a RA (or appropriate longitude) in degrees.')
                kill = True
            if self.dec == 9999.:
                print('Please pass a Dec (or appropriate latitude) in degrees.')
                kill = True
        except ValueError:
            pass
        if self.band == None:
            print('You did not pass a band string when passing magnitude.\nThe band has been set to Ks.')
        if self.frame == None:
            print('You did not pass a frame string when passing position.\nThe frame has been set to ICRS.')
        return kill

    def pass_parallax(self, par, err=None):
        '''Feed in a value of parallax (in mas) and an optional error on parallax (in mas).
        '''
        #Give a sanity check that we're dealing with mas
        if any(par < 0.01) or any(par > 5.):
            print(r'Are you 100% certain that (all) parallax(es) is/are in units of milliarcsec?')
        self.oo = par   #in mas
        self.oo_err = err #in mas

        if type(self.r) != type(None):
            print('Youve already passed in a value of distance.')
            print('Be mindful that these parallaxes will overwrite this as 1000./parallax, which is technically incorrect.')

        self.r = 1000./self.oo
        self.r_err = None
        if type(self.oo_err) != type(None):
            self.r_err = np.sqrt((-1000./self.oo)**2*self.oo_err**2)

    def pass_distance(self, dist, err=None):
        '''Feed in a properly treated value of distance (in pc) and an optional error on distance (in pc).
        '''
        self.r = dist #in pc
        self.r_err = err #in pc
        self.bailerjones = True #Using proper distances

    def pass_position(self, ra, dec, frame='icrs'):
        '''Feed in the position and the reference frame.

        Input:
            ra: the longitutde value (can be GLON for frame='galactic', for example)
            dec: the latitude value (can be GLAT for frame='galactic', for example)
        '''
        self.ra = ra
        self.dec = dec
        self.frame = frame

    def pass_magnitude(self, mag, err=None, band='Ks'):
        '''Pass in the apparent magnitues of a target in a given band.
        '''
        self.m = mag
        self.m_err = err
        self.band = band

    def get_mu0(self):
        '''Calculate the distance modulus from the parallax.
        '''
        return 5*np.log10(self.r) - 5

    def get_Ebv(self):
        '''Send a request to the online Bayestar catalogue for the extinction
        coefficients of all the targets.
        '''
        #Call the Bayestar catalogue
        bayestar = BayestarWebQuery(version='bayestar2017')
        coords = SkyCoord(self.ra.values*units.deg, self.dec.values*units.deg,
                distance=(self.r.values)*units.pc,frame=self.frame)
        #Find the extinction coefficient
        try:
            Ebv = bayestar(coords, mode='median')
        except:
            print('The Av values cant be downloaded for some reason.')
            print('No Av values for this star. Set Av to 0. for now.')
            Ebv = np.ones(self.ra.shape)
        return Ebv

    def get_Aband(self):
        '''Calculates the extinction coefficient in the chosen band by using the
        values from Green et al. 2018.
        '''
        #Check we have the conversion for this band of magnitudes
        if not any(self.band in head for head in list(Av_coeffs)):
            print('The class cant handle this specific band yet.')
            print('Please double check youve input the string correctly. The list of possible bands is:')
            print(list(coeffs))
            print('And you passed in: '+self.band)
            sys.exit()

        return self.get_Ebv() * Av_coeffs[self.band].values

    def get_error(self):
        '''Calculate the error on the absolute magnitude, given possible errors in
        apparent magnitude and parallax.
        '''
        #Case 0: No errors
        try:
            if (self.m_err is None) & (self.r_err is None):
                print('No errors given, error on M set to 0.')
                M_err = 0.
        except ValueError:
            pass

        #Case 1: Distance error only
        if self.m_err is None:
                M_err = np.sqrt((5/(self.r*np.log(10)))**2 * self.r_err**2)

        #Case 2: Magnitude error only
        elif self.oo_err is None:
            M_err = self.m_err

        #Case 3: Errors on both values
        else:
            M_err = np.sqrt(self.m_err**2 + (5/(self.r*np.log(10)))**2 * self.r_err**2)

        #Add the assumed error on Av
        return np.sqrt(M_err**2 + err_av**2)

    def get_M(self):
        '''Calculate the absolute magnitude given apparent magnitude, parallax
        and the Bayestar extiction coefficients.
        '''
        if self.check_contents():
            print('Im going to kill the run so you can pass the values correctly.')
            sys.exit()
            return None
        m = self.m       #Call in the magnitude
        mu0 = self.get_mu0()  #Call in the distance modulus
        Aband = self.get_Aband() #Call in the appropriate extinction coeff

        M_err = self.get_error()   #Propagate the error through

        return m - mu0 - Aband, M_err
