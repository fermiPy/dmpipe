#!/usr/bin/env python
#

""" 
Classes and utilities that manage spectral models
"""

__facility__ = "dmpipe.dmp_spectral"
__abstract__ = __doc__
__author__ = "E. Charles"
__date__ = "$Date: 2016/07/20 23:11:18 $"
__version__ = "$Revision: 1.2 $"
__release__ = "$Name:  $"

import sys
import os

import yaml


class SpectralLibrary(object):
    """
    """

    def __init__(self, spectral_dict={}):
        """
        """
        self.spectral_dict = spectral_dict

    def update(self, spectral_dict):
        """ """
        self.spectral_dict.update(spectral_dict)

    def __getitem__(self, key):
        """ """
        return self.spectral_dict.get(key, {})

    @staticmethod
    def create_from_yaml(yamlfile):
        """
        """
        spectral_dict = yaml.load(open(yamlfile))
        return SpectralLibrary(spectral_dict)
