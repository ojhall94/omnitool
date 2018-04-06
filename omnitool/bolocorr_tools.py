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

class bolometric:
    '''
    A class that calculates a bolometric correction for a star in a band, given
    its Temperature, logg, iron content, and mass.

    .. codeauthor:: Oliver James Hall
    '''
