#!/usr/bin/env python
#

"""
Utilities to fit dark matter spectra to castro data
"""

import os
import collections

import numpy as np
from fermipy import spectrum

def reverse_dict(din):
    """ Build and return a dict that maps value to key """
    dout = {}
    for key, val in din.items():
        for val2 in val:
            dout[val2] = key
    return dout


def find_gammac_file():
    """ Find the data file with the tables that DMFIT expects
    """
    inst_path = os.path.expandvars("${INST_DIR}/data/Likelihood/gammamc_dif.dat")
    if os.path.exists(inst_path):
        return inst_path
    base_path = os.path.expandvars("${BASE_DIR}/data/Likelihood/gammamc_dif.dat")
    if os.path.exists(base_path):
        return base_path
    raise RuntimeError("Missing DMFIT gammamc_dif.dat file")


class DMFitFunction(spectrum.SpectralFunction):
    """ Wrap the pylikelihood DMFitFunction interface.

    """
    channel_mapping = {
        1:  ["e+e-", "ee"],
        2:  ["mu+mu-", "mumu", "musrc"],
        3:  ["tau+tau-", "tautau", "tausrc"],
        4:  ["bb-bar", "bb", "bbbar", "bbsrc"],
        5:  ["tt-bar", "tt"],
        6:  ["gluons", "gg"],
        7:  ["W+W-", "w+w-", "ww", "wwsrc"],
        8:  ["ZZ", "zz"],
        9:  ["cc-bar", "cc"],
        10:  ["uu-bar", "uu"],
        11:  ["dd-bar", "dd"],
        12:  ["ss-bar", "ss"]}

    channel_rev_map = reverse_dict(channel_mapping)

    channel_tex = {
        1:  r'$e^{+}e^{-}$',
        2:  r'$\mu^{+}\mu^{-}$',
        3:  r'$\tau^{+}\tau^{-}$',
        4:  r'$b \bar b$',
        5:  r'$t \bar t$',
        6:  r'$gg$',
        7:  r'$W^{+}W^{-}$',
        8:  r'$ZZ$',
        9:  r'$c \bar c$',
        10:  r'$u \bar u$',
        11:  r'$d \bar d$',
        12:  r'$s \bar s$'}

    default_fit_par_names = ['sigmav', 'mass']

    default_params = dict(sigmav=1e-26, mass=100.,
                          norm=1e20, bratio=1.0,
                          channel0=4, channel1=1)

    default_gammac_file = find_gammac_file()
    base_gammac_file = '$(INST_DIR)/data/Likelihood/gammamc_dif.dat'

    static_parnames = None
    static_param_dict = None
    static_gammac_file = None
    static_dmf = None

    def __init__(self, params,
                 parnames=default_fit_par_names,
                 gammac_file=default_gammac_file,
                 default_pars=default_params):
        """ C'tor
        """
        super(DMFitFunction, self).__init__(params)
        if len(parnames) != len(params):
            raise ValueError("Number of parameters (%i) != number of parameter names (%i)" %
                             (len(parnames), len(params)))
        DMFitFunction.static_parnames = []
        DMFitFunction.static_parnames += parnames
        DMFitFunction.static_param_dict = default_pars.copy()
        for pname in parnames:
            DMFitFunction.static_param_dict.pop(pname)

        DMFitFunction.build_pylike(gammac_file)
        DMFitFunction._update(params)  # update all parameters in DMFitFunction

    @staticmethod
    def build_pylike(gammac_file):
        """ Build the pylikelihood DMFitFunction() object that evaluates the spectrum
        """
        import pyLikelihood

        if DMFitFunction.static_gammac_file == gammac_file:
            return DMFitFunction.static_dmf

        DMFitFunction.static_dmf = pyLikelihood.DMFitFunction()

        # unbound all parameters in gtlike
        for pname in DMFitFunction.static_param_dict.keys():
            DMFitFunction.static_dmf.getParam(pname).setBounds(-float('inf'), float('inf'))

            DMFitFunction.static_dmf.readFunction(os.path.expandvars(gammac_file))
            DMFitFunction.static_gamma_file = gammac_file

    @staticmethod
    def set_default_par_value(parname, parvalue):
        """ Set the default value for a parameter by name
        """
        DMFitFunction.static_param_dict[parname] = parvalue

    @staticmethod
    def _update(params):
        """ Update the DMFitFunction internally.
            This function should be called
            automatically when necessary.
        """
        if DMFitFunction.static_dmf is None:
            raise ValueError(
                "DMFitFunction instance must created before calling DMFitFunction._update()")

        for i, param_name in enumerate(DMFitFunction.static_parnames):
            if isinstance(params[i], collections.Iterable):
                val = float(np.squeeze(params[i]))
                DMFitFunction.static_dmf.setParam(param_name, val)
            else:
                try:
                    DMFitFunction.static_dmf.setParam(param_name, params[i])
                except:
                    raise ValueError("Failed to set parameter %s %.2f" % (param_name, params[i]))

        # Set the parameters which are not fixed explicitly
        for key, val in DMFitFunction.static_param_dict.items():
            DMFitFunction.static_dmf.setParam(key, val)

    @staticmethod
    def channel2int(cname):
        """ Convert a channel name to an integer """
        if DMFitFunction.channel_rev_map.has_key(cname):
            return DMFitFunction.channel_rev_map[cname]
        raise ValueError("Can't find channel named %s" % cname)

    @staticmethod
    def channel2tex(chan):
        """ Convert a channel name or integer to a laTex string """
        if DMFitFunction.channel_tex.has_key(chan):
            return DMFitFunction.channel_tex[chan]
        elif DMFitFunction.channel_rev_map.has_key(chan):
            return DMFitFunction.channel_tex[DMFitFunction.channel_rev_map[chan]]
        else:
            raise ValueError("Can't find channel %s" % chan)

    @staticmethod
    def int2channel(i):
        """ Convert an integer to a channel name """
        return DMFitFunction.channel_mapping[i][0]

    @staticmethod
    def channels():
        """ Return all available DMFit channel strings """
        return DMFitFunction.channel_rev_map.keys()

    @staticmethod
    def channels_unique():
        """ Return all available DMFit channel strings """
        return DMFitFunction.channels_unique

    @staticmethod
    def call_pylike_spectrum(spec, evals):
        """ Method to call a pylikelihood spectrum given
            either a python numer or a numpy array. """
        from pyLikelihood import dArg
        ex = evals
        if isinstance(ex, collections.Iterable):
            ret_val = np.zeros(ex.shape)
            for i, xval in np.ndenumerate(ex):
                ret_val[i] = spec(dArg(xval))
            return ret_val
            # return np.asarray([spec(dArg(i)) for i,x in np.ndenumerate(ex)])
        else:
            return spec(dArg(ex))

    @staticmethod
    def _eval_dnde(x, params, scale=1.0, extra_params=None):
        """ Return energy in MeV. This could be vectorized. """
        DMFitFunction._update(params)
        if isinstance(x, collections.Iterable):
            return np.asarray([DMFitFunction.call_pylike_spectrum(DMFitFunction.static_dmf, xx/scale)
                               for xx in x])
        else:
            return DMFitFunction.call_pylike_spectrum(DMFitFunction.static_dmf, x/scale)


if __name__ == "__main__":

    import argparse

    # Argument defintion
    USAGE = "usage: %(prog)s"
    DESCRIPTION = "Test program for dm_fit_spec.py"

    PARSER = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)

    # Argument parsing
    ARGS = PARSER.parse_args()

    INIT_PARS = [1e-26, 100.]
    DMF = DMFitFunction(INIT_PARS)
    E_BINS = np.logspace(2, 5, 13)
    E_VALS = np.sqrt(E_BINS[1:] * E_BINS[0:-1])

    DFDE = DMF.dfde(E_VALS)
    FLUX = DMF.flux(E_BINS[0:-1], E_BINS[1:])
    EFLUX = DMF.eflux(E_BINS[0:-1], E_BINS[1:])

    print ("dF/dE [ph MeV-1 cm^-2 s-1] : ", DFDE)
    print ("Flux  [ph cm^-2 s-1]       : ", FLUX)
    print ("EFlux [MeV cm^-2 s-1]      : ", EFLUX)
