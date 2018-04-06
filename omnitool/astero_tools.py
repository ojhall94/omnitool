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
    def __init__(self, _core_df, _numax, _dnu, _Teff, _numax_err = None, _dnu_err = None, _Teff_err = None):
        self.core_df = _core_df
        self.numax = _numax
        self.dnu = _dnu
        self.numax_err = _numax_err
        self.dnu_err = _dnu_err
        self.Teff = _Teff
        self.Teff_err = _Teff_err

    def get_radius(self):
        '''Calculates radius using the uncorrected asteroseismic scaling relations.
        '''
        R = Rsol * (self.numax / Numaxsol) * (self.dnu / Dnusol)**(-2) * (self.Teff / Tsol)**(0.5)
        return R

    def get_radius_err(self):
        try:
            term = (Rsol/Numaxsol)*(self.dnu/Dnusol)**(-2)*(self.Teff/Tsol)**(0.5)
            drdnumax = term**2 * self.numax_err**2
        except TypeError: drdnumax = 0.

        try:
            term = (Rsol/Dnusol**(-2))*(self.numax/Numaxsol)*(self.Teff/Tsol)**(0.5) * (-2*self.dnu**(-3))
            drdnu = term**2 * self.dnu_err**2
        except TypeError: drdnu = 0.

        try:
            term = (Rsol/Tsol**(0.5))*(self.numax/Numaxsol)*(self.dnu / Dnusol)**(-2) * 0.5*self.Teff**(-0.5)
            drdt = term**2 * self.Teff_err**2
        except TypeError: drdt = 0.

        sigR = np.sqrt(drdnumax + drdnu + drdt)
        return sigR

    def get_mass(self):
        '''Calculates mass using the uncorrected asteroseismic scaling relations.
        '''
        M = Msol * (self.numax / Numaxsol)**3 * (self.dnu / Dnusol)**(-4) * (self.Teff / Tsol)**(1.5)
        return M

    def get_mass_err(self):
        try:
            term = (Msol/Numaxsol**3)*(self.dnu/Dnusol)**(-4)*(self.Teff/Tsol)**(1.5) * 3*self.numax**2
            drdnumax = term**2 * self.numax_err**2
        except TypeError: drdnumax = 0.

        try:
            term = (Msol/Dnusol**(-4))*(self.numax/Numaxsol)**3*(self.Teff/Tsol)**(1.5) * (-4*self.dnu**(-5))
            drdnu = term**2 * self.dnu_err**2
        except TypeError: drdnu = 0.

        try:
            term = (Msol/Tsol**(1.5))*(self.numax/Numaxsol)**3*(self.dnu / Dnusol)**(-4) * 1.5*self.Teff**(0.5)
            drdt = term**2 * self.Teff_err**2
        except TypeError: drt = 0.

        sigM = np.sqrt(drdnumax + drdnu + drdt)
        return sigM

    def get_logg(self):
        '''Calculates log(g) using the uncorrected asteroseismic scaling relations.
        '''
        g = gsol * (self.numax/Numaxsol) * (self.Teff/Tsol)**0.5
        return np.log10(g)

    def get_logg_err(self):
        #First get error in g space
        try:
            term1 = ((gsol/Numaxsol)*(self.Teff/Tsol)**0.5) **2 * self.numax_err**2
        except TypeError: term1 = 0.
        try:
            term2 = ((gsol/Tsol**(0.5)) * (self.numax/Numaxsol) * 0.5*self.Teff**(-0.5))**2 * self.Teff_err**2
        except TypeError: term2 = 0.
        sigg = np.sqrt(term1 + term2)

        #Now convert to logg space
        g = gsol * (self.numax/Numaxsol) * (self.Teff/Tsol)**0.5
        siglogg = sigg / (g * np.log(10))

        return siglogg

    def get_luminosity(self):
        '''Calculates luminosity using the asteroseismically determined radius
        and given effective Temperature.
        '''
        L = 4 * np.pi * stefboltz * self.get_radius()**2 * self.Teff**4
        return L

    def get_luminosity_err(self):
        term1 = (8*np.pi*stefboltz*self.get_radius()*self.Teff**4)**2 * self.get_radius_err()**2
        try:
            term2 = (16*np.pi*stefboltz*self.get_radius()**2*self.Teff**3)**2 * self.Teff_err**2
        except TypeError: term2 = 0.

        sigL = np.sqrt(term1 + term2)
        return sigL

    def get_bolmag(self):
        '''Calculates the bolometric magnitude of the target given its
        luminosity.
        '''
        Mbol = -2.5*np.log10(self.get_luminosity()/Lzp)
        return Mbol

    def get_bolmag_err(self):
        nLum = self.get_luminosity()/Lsol
        nLume = self.get_luminosity_err()/Lsol
        sigMbol = np.sqrt( (-2.5/(nLum*np.log(10.)))**2*nLume**2)
        return sigMbol

    def get_M(self, BC, band='Ks'):
        '''Calculates the magnitude in a given band using an inverse bolometric
        correction.
        '''
        Mabs = self.get_bolmag() - BC
        return Mabs

    def get_M_err(self, BC_err band='Ks'):
        M_err = np.sqrt(self.get_bolmag_err()**2 + BC_err**2)
        return M_err
