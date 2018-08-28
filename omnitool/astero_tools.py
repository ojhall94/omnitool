#!/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys

from .literature_values import *

class scalings:
    '''
    An object class that stores asteroseismic and Temperature values (with error)
    and returns results from asteroseismic scaling relations upon request, as well
    as asteroseismic absolute magnitudes.

    .. codeauthor:: Oliver James Hall
    '''
    def __init__(self, _numax, _dnu, _Teff, _numax_err = None, _dnu_err = None, _Teff_err = None):
        self.numax = _numax
        self.dnu = _dnu
        self.numax_err = _numax_err
        self.dnu_err = _dnu_err
        self.Teff = _Teff
        self.Teff_err = _Teff_err
        self.Rbool = False
        self.Mbool = False
        self.fdnu = np.ones(len(self.numax))

    def give_corrections(self, fdnu = None, Mcorr=None, Mcorr_err=None, Rcorr=None, Rcorr_err=None):
        if type(fdnu) != type(None):
            self.fdnu = fdnu
            print('You have passed corrections to the Delta Nu scaling relation')

        if type(Mcorr) != type(None):
            self.Mcorr = Mcorr
            self.MCorr_err = Mcorr_err
            self.Mbool = True
            print('You have passed your own selection of masses.')

        if type(Rcorr) != type(None):
            self.Rcorr = Rcorr
            self.Rcorr_err = Rcorr_err
            self.Rbool = True
            print('You have a passed your own selection of radii.')

    def get_radius(self):
        '''Calculates radius using the uncorrected asteroseismic scaling relations.
        '''
        if self.Rbool:
            print('Please note that this does not return the corrected R')

        R = Rsol * (self.numax / Numaxsol) * (self.dnu / (self.fdnu * Dnusol))**(-2) * (self.Teff / Tsol)**(0.5)
        return R

    def get_radius_err(self):
        try:
            term = (Rsol/Numaxsol)*(self.dnu/(self.fdnu * Dnusol))**(-2)*(self.Teff/Tsol)**(0.5)
            drdnumax = term**2 * self.numax_err**2
        except TypeError: drdnumax = 0.

        try:
            term = (Rsol/((self.fdnu * Dnusol)**(-2)))*(self.numax/Numaxsol)*(self.Teff/Tsol)**(0.5) * (-2*self.dnu**(-3))
            drdnu = term**2 * self.dnu_err**2
        except TypeError: drdnu = 0.

        try:
            term = (Rsol/Tsol**(0.5))*(self.numax/Numaxsol)*(self.dnu / (self.fdnu * Dnusol))**(-2) * 0.5*self.Teff**(-0.5)
            drdt = term**2 * self.Teff_err**2
        except TypeError: drdt = 0.

        term = Rsol * (self.dnu/ (self.fdnu * Dnusol))**(-2) * (self.Teff/Tsol)**(0.5) * (-1*self.numax / (Numaxsol**2))
        drdnumaxsol = term**2 * eNumaxsol**2

        term = Rsol * (self.numax/Numaxsol) * (self.Teff/Tsol)**(0.5) * (2*self.fdnu*Dnusol)/(self.Dnu**2)
        drddnusol = term**2 * eDnusol**2

        sigR = np.sqrt(drdnumax + drdnu + drdt + drdnumaxsol + drddnusol)
        return sigR

    def get_mass(self):
        '''Calculates mass using the uncorrected asteroseismic scaling relations.
        '''
        if self.Mbool:
            print('Please note that this does not return the corrected R')

        M = Msol * (self.numax / Numaxsol)**3 * (self.dnu / (self.fdnu * Dnusol))**(-4) * (self.Teff / Tsol)**(1.5)
        return M

    def get_mass_err(self):
        try:
            term = (Msol/Numaxsol**3)*(self.dnu/(self.fdnu * Dnusol))**(-4)*(self.Teff/Tsol)**(1.5) * 3*self.numax**2
            dmdnumax = term**2 * self.numax_err**2
        except TypeError: dmdnumax = 0.

        try:
            term = (Msol/((self.fdnu * Dnusol)**(-4)))*(self.numax/Numaxsol)**3*(self.Teff/Tsol)**(1.5) * (-4*self.dnu**(-5))
            dmdnu = term**2 * self.dnu_err**2
        except TypeError: dmdnu = 0.

        try:
            term = (Msol/Tsol**(1.5))*(self.numax/Numaxsol)**3*(self.dnu /(self.fdnu * Dnusol))**(-4) * 1.5*self.Teff**(0.5)
            dmdt = term**2 * self.Teff_err**2
        except TypeError: dmt = 0.

        term = Msol * (self.dnu/(self.fdnu * Dnusol))**(-4)*(self.Teff/Tsol)**(1.5) * (-3*self.numax**3)/(Numaxsol**4)
        dmdnumaxsol = term**2 * eNumaxsol**2

        term = Msol * (self.numax/Numaxsol)**3*(self.Teff/Tsol)**(1.5) * (4.*self.fdnu*Dnusol**3)/(self.dnu**4)
        dmddnusol = term**2 * eDnusol**2

        sigM = np.sqrt(dmdnumax + dmdnu + dmdt + dmdnumaxsol + dmddnusol)
        return sigM

    def get_logg(self):
        '''Calculates log(g) using the uncorrected asteroseismic scaling relations.
        '''
        g = gsol * (self.numax/Numaxsol) * (self.Teff/Tsol)**0.5
        return np.log10(g)

    def get_logg_err(self):
        #First get error in g space
        try:
            dgdnumax = ((gsol/Numaxsol)*(self.Teff/Tsol)**0.5)**2 * self.numax_err**2
        except TypeError: term1 = 0.
        try:
            dgdteff = ((gsol/Tsol**(0.5)) * (self.numax/Numaxsol) * 0.5*self.Teff**(-0.5))**2 * self.Teff_err**2
        except TypeError: term2 = 0.

        dgdnumaxsol = (gsol * self.numax * (self.Teff/Tsol)**0.5 * (-1./Numaxsol**2))**2 * eNumaxsol**2

        sigg = np.sqrt(dgdnumax + dgdteff + dgdnumaxsol)

        #Now convert to logg space
        g = gsol * (self.numax/Numaxsol) * (self.Teff/Tsol)**0.5
        siglogg = sigg / (g * np.log(10))

        return siglogg

    def get_luminosity(self):
        '''Calculates luminosity using the asteroseismically determined radius
        and given effective Temperature.
        '''
        if self.Rbool:
            print('Calculating L using a given radius')
            L = 4 * np.pi * stefboltz * self.Rcorr**2 * self.Teff**4

        else:
            print('Calculating luminosity using basic asteroseismic radius')
            L = 4 * np.pi * stefboltz * self.get_radius()**2 * self.Teff**4

        return L

    def get_luminosity_err(self, give_radius_err = ''):
        if self.Rbool:
            R = self.Rcorr
            Rerr = self.Rcorr_err
        else:
            R = self.get_radius()
            Rerr = self.get_radius_err()

        term1 = (8*np.pi*stefboltz*R*self.Teff**4)**2 * Rerr**2
        try:
            term2 = (16*np.pi*stefboltz*R**2*self.Teff**3)**2 * self.Teff_err**2
        except TypeError: term2 = 0.

        sigL = np.sqrt(term1 + term2)
        return sigL

    def get_bolmag(self):
        '''Calculates the bolometric magnitude of the target given its
        luminosity.
        '''
        # Mbol = -2.5*np.log10(self.get_luminosity()/Lzp) #IAU accepted measure
        Mbol = -2.5 * np.log10(self.get_luminosity()/Lsol) + Mbolsol #Casagrande & Vandenburg proposed measure
        return Mbol

    def get_bolmag_err(self):
        nLum = self.get_luminosity()/Lsol
        nLume = self.get_luminosity_err()/Lsol
        sigMbol = np.sqrt( (-2.5/(nLum*np.log(10.)))**2*nLume**2)
        return sigMbol
