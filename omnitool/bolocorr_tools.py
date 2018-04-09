#!/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sys

from .literature_values import *

from astropy.table import Table
from scipy import integrate
from pystellibs import Kurucz
from synphot import SpectralElement

from tqdm import tqdm

class bolometric_correction:
    '''
    A class that calculates a bolometric correction for a star in a band, given
    its Temperature, logg, iron content, and mass.

    .. codeauthor:: Oliver James Hall
    '''
    def __init__(self, _Teff, _logg, _L, _Z):
        self.Teff = _Teff
        self.logg = _logg
        self.L = _L
        self.Z = _Z

    def get_spectra(self):
        '''Find synthetic spectra for each target given its Temperature, logg,
        Luminosity and metallicity. This makes use of the pystellibs package, which
        integrates from a choice of synthetic spectra (in this case, Kurucz)
        '''
        self.spect = Kurucz() #Call in the Kurucz spectrum

        # Define the interpolation data
        logteff = np.log10(self.Teff)
        logg = self.logg
        logL = np.log10(self.L)
        Z = self.Z

        #Find the spectra
        ap = (logteff, logg, logL, Z)
        ap = Table(ap, names=('logT','logg','logL','Z'))

        _, spectra = self.spect.generate_individual_spectra(ap)
        self.spectra = np.array(spectra)

    def set_band(self, band):
        '''You can use this line to set the band of the bolometric correction
        dynamically for a hardcoded number of bands.
        '''
        if band == 'Ks':
            self.B = SpectralElement.from_filter('bessel_k')
        if band == 'J':
            self.B = SpectralElement.from_filter('bessel_j')
        if band == 'H':
            self.B = SpectralElement.from_filter('bessel_h')
        if band == 'Rc':
            self.B = SpectralElement.from_filter('cousins_r')
        if band == 'Ic':
            self.B = SpectralElement.from_filter('cousins_i')
        self.band = band

    def integrate_spectra(self):
        '''This integrates the area under the spectra for the given passband.
        '''
        #Follow the equation prescribed in Torres 2010
        fl = np.zeros(len(self.Teff))
        slfl = np.zeros(len(self.Teff))

        for idx in tqdm(range(len(self.Teff))):
            fl[idx] = integrate.simps(self.spectra[idx], x=self.spect._wavelength)
            slfl[idx] = integrate.simps(self.spectra[idx] * self.B(self.spect._wavelength), x=self.spect._wavelength)

        return fl, slfl

    def solve_for_C(self):
        '''This calculates the integration constant for the bolometric correction
        integration function, using known Solar values.
        '''
        #Build solar spectrum
        ap = (np.log10(Tsol), np.log10(gsol), np.log10(Lsol), Zsol)
        sb = self.spect.generate_stellar_spectrum(*ap)

        #Perform the integration
        fl = integrate.simps(sb, x=self.spect._wavelength)
        slfl = integrate.simps(sb * self.B(self.spect._wavelength), x=self.spect._wavelength)

        #Calculate the integration constant given known values
        C = Mbolsol - Mbandsol[self.band] - 2.5*np.log10(slfl/fl)

        return C.values[0]

    def calculate_bolocorr(self):
        '''This calculates the bolometric correction for a given passband, given
        proper calculation of the integration constant for the solar values.
        '''
        fl, slfl = self.integrate_spectra()
        C = self.solve_for_C()
        BC = 2.5*np.log10(slfl/fl) + C
        return BC

    def __call__(self, band='Ks'):
        self.set_band(band)
        self.get_spectra()
        BC = self.calculate_bolocorr()
        return BC
