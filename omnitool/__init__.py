#!/bin/env python
# -*- coding: utf-8 -*-

"""
A Python package that compiles any classes I've written to help with
asteroseismic or photmetric analysis of data

.. versioncreated:: 0.1

.. codeauthor:: Oliver James Hall <ojh251@student.bham.ac.uk>
"""

from .spyglass import spyglass
from .astero_tools import scalings
from .bolocorr_tools import bolometric_correction
from .literature_values import *
