#!/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

import sys

from .literature_values import *

class bailerjones:
    '''
    A class that takes parallax, parallax error and a chosen value L for the
    length scale on the exponentially decaying space density prior, and returns
    the mode of the posterior probability distribution to obtain the dsitance given
    the parallax and parallax uncertainty as a distance estimate.

    For more information, see Bailer-Jones15, Bailer-Jones+18, and/or the
    codeauthors blog post on the topic: ojhall94.github.io/blogs/distances_blog.html

    .. codeauthor:: Oliver James Hall <ojh251@student.bham.ac.uk>
    '''
    def __init__(self, _oo, _oo_err, _L):
        '''
        Parameters:
            _oo (ndarray): A single float or array containing parallax values in
                milliarcsec (mas)(!)
            _oo_err (ndarray): A single float or array containing uncertainty on
                parallax values in milliarcsec (mas)(!)

            _L (ndarray): A single float or array with the values for the length
                scale of the exponentially decaying space density prior.

        Returns:
            ndarray: Estimated mode distance in parsec for each target.
        '''
        self.oo = _oo/1000.
        self.oo_err = _oo_err/1000.
        self.L = _L*np.ones_like(_oo)

    def get_roots(self,L, oo, oo_err):
        '''Calculate the mode distance using the roots for a given L, oo and oo_err.'''
        p = np.array([1./L, -2, oo/(oo_err**2.), -1./(oo_err**2.)])
        fullroots = np.roots(p) #Calculating the roots of the derivative of the posterior distribution
        roots = fullroots[np.isreal(fullroots)] #Make sure we take only the true values

        #This follows the prescriptions set out in Bailer-Jones 2015
        if len(roots) == 1:
            return float(roots[0])
        if len(roots) == 3:
            if oo >= 0.:
                r = float(np.min(roots))
            if oo < 0.:
                r = float(roots[roots > 0][0])
            return r
        else:
            print('You shouldnt be here, printing roots below for diagnostic:')
            print(roots)
            print('L should be > 0, as should oo_err. Printing both below:')
            print('L: '+str(np.round(L,2))+' | oo_err: '+str(np.round(oo_err,2)))

    def get_modes(self):
        '''Call the get_roots function one by one for the given data.'''
        try:
            if len(self.oo) > 1:
                return np.array([self.get_roots(L, o, err) for o, err, L in zip(self.oo, self.oo_err, self.L)])
        except TypeError:
            return np.array(self.get_roots(self.L, self.oo, self.oo_err))

    def __call__(self):
        return self.get_modes()
