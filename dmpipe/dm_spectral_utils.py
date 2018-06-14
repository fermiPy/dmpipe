#!/usr/bin/env python
#

"""
Interface to Dark Matter spectra
"""
from __future__ import absolute_import, division, print_function

import sys

import yaml
import numpy as np

from astropy.table import Table, Column

from dmsky.priors import create_prior_functor

from fermipy import castro
from fermipy import fits_utils
from fermipy import skymap
from fermipy.utils import load_yaml
from fermipy.jobs.utils import is_null, is_not_null

from fermipy.spectrum import DMFitFunction

from dmpipe.lnl_norm_prior import LnLFn_norm_prior

REF_J = 1.0e19
REF_SIGV = 1.0e-26
REF_D = 1.0e17
REF_TAU = 1.0e26

class DMCastroData(castro.CastroData_Base):
    """ This class wraps the data needed to make a "Castro" plot,
    namely the log-likelihood as a function of <sigma v> (or tau) * J (or D) for a
    series of DM masses
    """

    def __init__(self, norm_vals, nll_vals,
                 channel, masses, astro_value,
                 **kwargs):
        """ C'tor

        Parameters
        ----------
        norm_vals : `~numpy.ndarray`
           The normalization values ( n_mass X N array, where N is the
           number of sampled values for each bin )
           Note that these should be the true values, with the 
           reference J-value (or D-value) included, and _NOT_ the values w.r.t. to the 
           reference J-value (or D-value) spectrum.

        nll_vals : `~numpy.ndarray`
           The log-likelihood values ( n_mass X N array, where N is
           the number of sampled values for each bin )

        channel : int or str
           The DM decay or annihilation channel

        masses : '~numpy.ndarray'
           The masses (in GeV)

        astro_value : float
           The J-factor (or D-factor) of the target


        Keyword arguments
        -----------------

        astro_prior : `dmsky.priors.PriorFunctor`
           The prior on the J-factor (or D-factor)

        prior_applied : bool
           If true, then the prior has already been applied

        ref_astro : float
           Reference value for J-factor (or D-factor)

        ref_inter : float
           Reference value for sigmav (or tau)

        norm_type : str
           Normalization type: 'sigmav' or 'tau' or 'norm'
        
        decay : bool
           Trye for decay instead of annihilation

        """
        self._masses = masses
        self._astro_value = astro_value
        self._astro_prior = kwargs.get('astro_prior', None)
        self._prior_applied = kwargs.get('prior_applied', True)
        self._decay = kwargs.get('decay', False)
        if self._decay:
            self._ref_astro = kwargs.get('ref_astro', REF_D)
            self._ref_inter = kwargs.get('ref_inter', REF_TAU)
            norm_type = kwargs.get('norm_type', 'tau')
        else:
            self._ref_astro = kwargs.get('ref_astro', REF_J)
            self._ref_inter = kwargs.get('ref_inter', REF_SIGV)
            norm_type = kwargs.get('norm_type', 'sigmav')

        if isinstance(channel, str):
            self._channel = DMFitFunction.channel_rev_map[channel]
        else:
            self._channel = channel
        super(DMCastroData, self).__init__(norm_vals,
                                           nll_vals,
                                           norm_type=norm_type)

    @property
    def n_masses(self):
        """ The number of masses tested """
        return self._nx

    @property
    def masses(self):
        """ The masses tested (in GeV) """
        return self._masses

    @property
    def astro_value(self):
        """ The astrophysical J or D factor """
        return self._astro_value

    @property
    def astro_prior(self):
        """ The prior on the astrophysical J or D factor """
        return self._astro_prior

    @property
    def prior_mean(self):
        """ The mean on the prior on the astrophysical J or D factor.

        Note that this is a multiplicative prior, so we except a mean of 1
        """
        if self._astro_prior is None:
            return 1.0
        return self.astro_prior.mean()

    @property
    def prior_sigma(self):
        """ The width of the prior on the astrophysical J or D factor.

        This is actually the width parameter given to the prior function.
        The exact meaning depends on the function type being used.
        """
        if self._astro_prior is None:
            return 0.0
        return self.astro_prior.sigma()

    @property
    def prior_type(self):
        """ The function type of the astrophysical J or D factor.

        See '~fermipy.stats_utils.create_prior_functor' for recognized types.
        """
        if self._astro_prior is None:
            return "none"
        return self.astro_prior.funcname

    @property
    def prior_applied(self):
        """ Has the prior already been applied """
        return self._prior_applied

    @property
    def channel(self):
        """ The string specifying the decay channel """
        return self._channel

    @property
    def ref_astro(self):
        """ Reference value for J (or D) """
        return self._ref_astro

    @property
    def ref_inter(self):
        """ Reference value for <sigmav> (or tau)"""
        return self._ref_inter

    @property
    def decay(self):
        """ Return True for decay instead of annihilation"""
        return self._decay

    @classmethod
    def create_from_yamlfile(cls, yamlfile, channel, prior=None, decay=False):
        """ Create a DMCastroData object from a yaml file 
 
        Parameters
        ----------
        yamlfile : str
           Input file

        channel : str
           The DM decay or annihilation channel

        prior : str or None
           Type of prior on the J-factor (or D-factor)

        decay : bool
           True for decay instead of annihilation

        Returns
        -------

        output : `DMCastroData`
            Object with the DM-space likelihoods

        """
        data = load_yaml(yamlfile)
        masses = np.array([float(v) for v in data['param']])

        if decay:
            astro_str = 'dvalue'
            rvalues_str = 'rdvalues'
            sigma_str = 'dsigma'
            ref_str = 'd_ref'
            norm_type = 'tau'
        else:
            astro_str = 'jvalue'
            rvalues_str = 'rjvalues'
            sigma_str = 'jsigma'
            ref_str = 'j_ref'
            norm_type = 'sigmav'

        if prior is not None:
            try:
                astro_value = data[rvalues_str][prior]
                sigma = data[sigma_str]
            except KeyError:
                astro_value = 1.
                sigma = 1.
            lnlstr = 'p1lnl'
            prior_dict = dict(functype=jprior,
                              mu=astro_value,
                              sigma=sigma)
            prior_dict[ref_str] = astro_value
            prior = create_prior_functor(prior_dict)
        else:
            try:
                astro_value = data[astro_str]
            except KeyError:
                astro_value = 1.            
            jsigma = None
            lnlstr = 'lnl'
            prior = None
            
        prior_applied = True
        
        norm_list = []
        nll_list = []

        masses =  np.array(sorted ([ float(v) for v in data['param'] ]))
        masses_st =[ "%0.1f" % v for v in masses ]
        for mass in masses_st:
            norm_list.append(data['lnldata'][mass]['norm'])
            ll_vals = data['lnldata'][mass][lnlstr]
            nll_vals = ll_vals.max() - ll_vals
            nll_list.append(nll_vals)
            
        norm_vals = np.vstack(norm_list)
        nll_vals = np.vstack(nll_list)

        return cls(norm_vals, nll_vals, channel, masses, astro_value,
                   astro_prior=prior, 
                   prior_applied=prior_applied,
                   ref_astro=astro_value,
                   norm_type=norm_type,
                   decay=decay)


    @classmethod
    def create_from_stack(cls, components, 
                          **kwargs):
        """ Create a DMCastroData object by stacking a series of DMCastroData objects

        Parameters
        ----------
        components : list
           List of `DMCastroData` objects we are stacking

        Keyword Arguments
        -----------------

        nystep : int
            Number of steps in <sigmav> to take in sampling.

        ylims : tuple
            Limits of range of <sigmav> to scan

        weights : list or None
            List of weights to apply to components

        ref_astro : float
            Reference J-factor value

        ref_inter : float
            Refernece <sigmav> value

        decay : bool
            True for decay instead of annihilation
            
        Returns
        -------

        output : `DMCastroData`
            Object with the DM-space likelihoods

        """
        if not components:
            return None

        nystep = kwargs.get('nystep', 200)
        weights = kwargs.get('weights', None)
        decay = kwargs.get('decay', False)
        if decay:
            ref_astro = kwargs.get('ref_astro', REF_D)
            ref_inter = kwargs.get('ref_inter', REF_TAU)
            ylims = kwargs.get('nystep', (1e+20, 1e+30))
        else:
            ref_astro = kwargs.get('ref_astro', REF_J)
            ref_inter = kwargs.get('ref_inter', REF_SIGV)
            ylims = kwargs.get('nystep', (1e-30, 1e-20))

        shape = (components[0].nx, nystep)
        norm_vals, nll_vals = castro.CastroData_Base.stack_nll(shape, components,
                                                               ylims, weights)
        return cls(norm_vals, nll_vals, components[0].channel, components[0].masses,
                   astro_value=None, ref_astro=ref_astro, ref_inter=ref_inter, 
                   norm_type='norm', decay=decay)

    @classmethod
    def create_from_tables(cls, tab_s, tab_m, norm_type, decay=False):
        """ Create a DMCastroData object from likelihood scan and mass tables

        Parameters
        ----------

        tab_s : `astropy.table.Table`
            Table with likelihood scan data

        tab_s : `astropy.table.Table`
            Table with masses of corresponding spectra

        norm_type : str
            Type of normalization to use.  Valid options are:

            * norm : Self-normalized
            * sigmav : Reference values of sigmav and J
            * tau : Reference values of tau and D
            
        decay : bool
            True for decay instead of annihilation

        Returns
        -------

        output : `DMCastroData`
            Object with the DM-space likelihoods

        """
        if decay:
            astro_str = 'ref_D'
            inter_str = 'ref_tau'
        else:
            astro_str = 'ref_J'
            inter_str = 'ref_sigmav'

        ref_astro = np.squeeze(np.array(tab_s[astro_str]))
        ref_inter = np.squeeze(np.array(tab_s[inter_str]))

        if norm_type in ['sigmav', 'tau']:
            norm_vals = np.squeeze(np.array(tab_s['norm_scan']*ref_inter))
        elif norm_type in ['norm']:
            norm_vals = np.squeeze(np.array(tab_s['norm_scan']))
        else:
            raise ValueError('Unrecognized normalization type: %s' % norm_type)

        nll_vals = -np.squeeze(np.array(tab_s['dloglike_scan']))

        masses = np.squeeze(np.array(tab_m['masses']))
        channel = np.squeeze(np.array(tab_m['channel']))

        astro_value = np.squeeze(np.array(tab_s['astro_value']))
        try:
            astro_priortype = tab_s['prior_type']
        except KeyError:
            astro_priortype = None

        if astro_priortype is not None:
            prior_mean = np.squeeze(np.array(tab_s['prior_mean']))
            prior_sigma = np.squeeze(np.array(tab_s['prior_sigma']))
            prior_applied = np.squeeze(np.array(tab_s['prior_applied']))
            prior_dict = dict(functype=astro_priortype,
                              mu=prior_mean,
                              sigma=prior_sigma)
            prior_dict[astro_str] = astro_value
            prior = create_prior_functor(prior_dict)
        else:
            prior = None
            prior_applied = True

        return cls(norm_vals, nll_vals, channel,
                   masses, astro_value,
                   astro_prior=prior, 
                   prior_applied=prior_applied,
                   ref_astro=ref_astro, ref_inter=ref_inter,
                   norm_type=norm_type, decay=decay)


    @classmethod
    def create_from_fitsfile(cls, filepath, channel, norm_type=None):
        """ Create a DMCastroData object likelihood scan and mass tables in FITS file

        Parameters
        ----------

        filepath : str
            Path to likelihood scan data.

        channel : str
            DM interaction channel
 
        norm_type : str
            Type of normalization to use.  Valid options are:

            * norm : Self-normalized
            * sigmav : Reference values of sigmav and J
            * tau : Reference values of tau and D

        Returns
        -------

        output : `DMCastroData`
            Object with the DM-space likelihoods

        """
        decay = channel.find('_decay') >= 0
        if norm_type is None:
            if decay:
                norm_type = 'tau'
            else:
                norm_type = 'sigmav'
        tab_s = Table.read(filepath, hdu=channel)
        tab_m = Table.read(filepath, hdu='masses')
        return cls.create_from_tables(tab_s, tab_m,
                                      norm_type=norm_type, 
                                      decay=decay)


    def build_lnl_fn(self, normv, nllv):
        """ Build a function to return the likelihood value arrays of
        normalization and likelihood values.

        Parameters
        ----------

        normv : `numpy.array`
            Set of test normalization values

        nllv : `numpy.array`
            Corresponding set of negative log-likelihood values

        Returns
        -------

        output : FIXME
            Object that can compute the negative log-likelihood with
            the prior included.

        """
        lnlfn = castro.LnLFn(normv, nllv, self._norm_type)
        if self._astro_prior is None:
            return lnlfn
        if self._prior_applied:
            return lnlfn
        return LnLFn_norm_prior(lnlfn, self._astro_prior)

    def build_scandata_table(self, norm_type=None):
        """Build a FITS table with likelihood scan data

        Parameters
        ----------
        
        norm_type : str or None
            Type of normalization to use.  Valid options are:

            * norm : Self-normalized
            * sigmav : Reference values of sigmav and J
            * tau : Reference values of tau and D

        Returns
        -------

        table : `astropy.table.Table`
            The table has these columns

        astro_value : float
            The astrophysical J-factor (or D-factor) for this target
        ref_astro : float
            The reference J-factor (or D-factor) used to build `DMSpecTable`
        ref_inter : float
            The reference <sigmav> (or tau) used to build `DMSpecTable`

        norm_scan : array
            The test values of <sigmav> (or tau)
        dloglike_scan : array
            The corresponding values of the negative log-likelihood

        """
        if self.decay:
            astro_str = 'ref_D'
            inter_str = 'ref_tau'
            if norm_type is None:
                norm_type = 'tau'
        else:
            astro_str = 'ref_J'
            inter_str = 'ref_sigmav'
            if norm_type is None:
                norm_type = 'sigmav'

        shape = self._norm_vals.shape
        #dtype = 'f%i'%self._norm_vals.size
        
        col_normv = Column(name="norm_scan", dtype=float,
                           shape=shape)
        col_dll = Column(name="dloglike_scan", dtype=float,
                         shape=shape)

        col_astro_val = Column(name="astro_value", dtype=float)
        col_ref_astro = Column(name=astro_str, dtype=float)
        col_ref_inter = Column(name=inter_str, dtype=float)

        collist = [col_normv, col_dll, col_astro_val, col_ref_astro,
                   col_ref_inter]

        if norm_type in ['sigmav', 'tau']:
            norm_vals = self._norm_vals / self.ref_inter
        elif norm_type in ['norm']:
            norm_vals = self._norm_vals
        else:
            raise ValueError('Unrecognized normalization type: %s' % norm_type)

        valdict = {"norm_scan": norm_vals,
                   "dloglike_scan": -1 * self._nll_vals,
                   "astro_value": self.astro_value,
                   astro_str: self.ref_astro,
                   inter_str: self.ref_inter}

        if self._astro_prior is not None:
            col_prior_type = Column(name="prior_type", dtype="S16")
            col_prior_mean = Column(name="prior_mean", dtype=float)
            col_prior_sigma = Column(name="prior_sigma", dtype=float)
            col_prior_applied = Column(name="prior_applied", dtype=bool)
            collist += [col_prior_type, col_prior_mean,
                        col_prior_sigma, col_prior_applied]
            valdict["prior_type"] = self.prior_type
            valdict["prior_mean"] = self.prior_mean
            valdict["prior_sigma"] = self.prior_sigma
            valdict["prior_applied"] = self.prior_applied

        tab = Table(data=collist)
        tab.add_row(valdict)
        return tab

    def build_mass_table(self):
        """Build a FITS table with mass values

        Returns
        -------

        table : `astropy.table.Table`
            The table has these columns

        masses : array
            The masses of the spectra, in GeV
        channel : int
            The index of the channel of the DM interaction

        """
        col_masses = Column(name="masses", dtype=float,
                            shape=self._masses.shape)

        col_channel = Column(name="channel", dtype=int)
        tab = Table(data=[col_masses, col_channel])
        tab.add_row({"masses": self._masses,
                     "channel": self._channel})
        return tab

    def build_limits_table(self, limit_dict):
        """Build a FITS table with limits data

        Paramters
        ---------

        limit_dict : dict
            Dictionary from limit names to values


        Returns
        -------

        table : `astropy.table.Table`
            The table has these columns

        astro_value : float
            The astrophysical J-factor for this target
        ref_astro : float
            The reference J-factor used to build `DMSpecTable`
        ref_inter : float
            The reference <sigmav> used to build `DMSpecTable`
        <LIMIT> : array
            The upper limits

        If a prior was applied these additional colums will be present

        prior_type : str
            Key specifying what kind of prior was applied
        prior_mean : float
            Central value for the prior
        prior_sigma : float
            Width of the prior
        prior_applied : bool
            Flag to indicate that the prior was applied

        """
        if self.decay:
            astro_str = "ref_D"
            inter_str = "ref_tau"
        else:
            astro_str = "ref_J"
            inter_str = "ref_sigmav"

        col_astro_val = Column(name="astro_value", dtype=float)
        col_ref_astro = Column(name=astro_str, dtype=float)
        col_ref_inter = Column(name=inter_str, dtype=float)
        collist = [col_astro_val, col_ref_astro, col_ref_inter]
        valdict = {"astro_value": self.astro_value,
                   astro_str: self.ref_astro,
                   inter_str: self.ref_inter}

        for k, v in limit_dict.items():
            collist.append(Column(name=k, dtype=float, shape=v.shape))
            valdict[k] = v

        if self._astro_prior is not None:
            col_prior_type = Column(name="prior_type", dtype="S16")
            col_prior_mean = Column(name="prior_mean", dtype=float)
            col_prior_sigma = Column(name="prior_sigma", dtype=float)
            col_prior_applied = Column(name="prior_applied", dtype=bool)
            collist += [col_prior_type, col_prior_mean,
                        col_prior_sigma, col_prior_applied]
            valdict["prior_type"] = self.prior_type
            valdict["prior_mean"] = self.prior_mean
            valdict["prior_sigma"] = self.prior_sigma
            valdict["prior_applied"] = self.prior_applied

        tab = Table(data=collist)
        tab.add_row(valdict)
        return tab

    
    def x_edges(self):
        """ Make a reasonable set of bin edges for plotting 

        To do this we expand relative to the mass points by half the bid width in either direction
        """
        log_masses = np.log10(self.masses)
        log_half_widths = (log_masses[1:] - log_masses[0:-1]) /2
        last_mass = log_masses[-1] + log_half_widths[-1]
        log_masses[0:-1] -= log_half_widths
        log_masses[-1] -= log_half_widths[-1]
        log_masses = np.append(log_masses, last_mass)
        return np.power(10, log_masses)



class DMSpecTable(object):
    """ Version of the DM spectral tables in tabular form
    """

    def __init__(self, e_table, s_table, ref_vals):
        """ C'tor to build this object from energy binning and spectral values tables.
        """
        self._e_table = e_table
        self._s_table = s_table
        self._ref_vals = ref_vals
        self._channel_map, self._channel_names = DMSpecTable.make_channel_map(s_table)

    @staticmethod
    def make_channel_map(s_table):
        """ Construct a map from channel number to list of table rows """
        data = s_table['ref_chan']
        o = {}
        for i, d in enumerate(data):
            if d in o:
                o[d].append(i)
            else:
                o[d] = [i]
        l = []
        for k in sorted(o.keys()):
            is_decay = k >= 100
            if is_decay:
                nm = DMFitFunction.channel_shortname_mapping[k-100]
                nm += '_decay'
            else:
                nm = DMFitFunction.channel_shortname_mapping[k]
            l.append(nm)
        return o, l

    @staticmethod
    def make_ebounds_table(emin, emax, eref):
        """ Construct the energy bounds table

        Returns
        -------

        table : `astropy.table.Table`
            The table has these columns and one row per energy bin

        E_MIN : float
            Energy bin lower edge

        E_MAX : float
            Energy bin upper edge

        E_REF : float
            Reference energy for bin, typically geometric mean
            of bin edges

        """
        from astropy import table

        col_emin = table.Column(
            name="E_MIN", dtype=float, unit="MeV", data=emin)
        col_emax = table.Column(
            name="E_MAX", dtype=float, unit="MeV", data=emax)
        col_eref = table.Column(
            name="E_REF", dtype=float, unit="MeV", data=eref)

        tab = table.Table(data=[col_emin, col_emax, col_eref])
        return tab

    @staticmethod
    def make_spectra_table(nebins, spec_dict):
        """Construct the spectral values table

        Returns
        -------

        table : `astropy.table.Table`
            The table has these columns

        ref_mass : float
            The DM mass for this row in the table

        ref_chan : int
            The index of the DM interaction channel for this row
 
        ref_dnde : array
            The reference differential photon flux fpr each energy [ph / (MeV cm2 s)]

        ref_flux : array
            The reference integral photon flux for each energy [ph / (cm2 s)]

        ref_eflux : array
            The reference integral energy flux for each energy [MeV / (cm2 s)]

        """
        from astropy import table

        col_masses = table.Column(name="ref_mass", dtype=float, unit="GeV",
                                  data=spec_dict['mass'])
        col_chans = table.Column(name="ref_chan", dtype=int,
                                 data=spec_dict['chan'])
        col_dnde = table.Column(name="ref_dnde", dtype=float, shape=nebins, unit="ph / (MeV cm2 s)",
                                data=spec_dict['dnde'])
        col_flux = table.Column(name="ref_flux", dtype=float, shape=nebins, unit="ph / (cm2 s)",
                                data=spec_dict['flux'])
        col_eflux = table.Column(name="ref_eflux", dtype=float, shape=nebins, unit="MeV / (cm2 s)",
                                 data=spec_dict['eflux'])

        table = table.Table(
            data=[col_masses, col_chans, col_dnde, col_flux, col_eflux])
        return table

    @property
    def ebounds_table(self):
        """Return the energy binning table """
        return self._e_table

    @property
    def spectra_table(self):
        """Return the spectral values table """
        return self._s_table

    @property
    def ref_vals(self):
        """Return the spectral reference values """
        return self._ref_vals

    @property
    def channel_map(self):
        """Return the channel to index mapping """
        return self._channel_map

    @property
    def channel_names(self):
        """Return the channel to index mapping """
        return self._channel_names


    def spectrum(self, chan, mass, spec_type):
        """Return the spectrum for a particular channel, mass and spectral type
        """
        mask = (self._s_table["ref_chan"] == chan) & (
            np.abs(self._s_table["ref_mass"] - mass) < 1e-9)
        spec_vals = self._s_table[mask]["ref_%s" % spec_type].data
        return spec_vals

    def masses(self, chan):
        """Return the array of masses for a given channel
        """
        mask = (self._s_table["ref_chan"] == chan)
        return self._s_table[mask]["ref_mass"]

    def ebin_edges(self):
        """Return an array with the energy bin edges
        """
        return np.hstack([self._e_table["E_MIN"].data,
                          self._e_table["E_MAX"].data])

    def ebin_refs(self):
        """Return an array with the energy bin reference energies
        """
        return self._e_table["E_REF"].data

    def write_fits(self, filepath, clobber=False):
        """ Write this object to a FITS file

        Paramters
        ---------

        filepath : str
            Path to output file

        clobber : bool
            Flag to allow overwriting existing files

        """
        fits_utils.write_tables_to_fits(filepath, [self._e_table, self._s_table],
                                        clobber=clobber,
                                        namelist=["EBOUNDS", "SPECDATA"],
                                        cardslist=[{}, self._ref_vals])

    @staticmethod
    def get_ref_vals(filepath):
        """ Extract the reference values from a FITS header

        Paramters
        ---------

        filepath : str
            Path to input file


        Returns
        -------

        ref_sigv : float
            Reference value of <sigmav>

        ref_j : float
            Refernce value of J-factor

        """
        import astropy.io.fits as pf
        fin = pf.open(filepath)
        hdu = fin["SPECDATA"]
        hin = hdu.header
        dref = {"REF_SIGV": hin["REF_SIGV"],
                "REF_J": hin["REF_J"],
                "REF_TAU": hin["REF_TAU"],
                "REF_D": hin["REF_D"]}
        fin.close()
        return dref

    @staticmethod
    def create_from_fits(filepath):
        """ Build a DMSpecTable object from a FITS file

        Paramters
        ---------

        filepath : str
            Path to the input file

        Returns
        -------

        output : `DMSpecTable`
            The newly created object

        """
        from astropy import table
        e_table = table.Table.read(filepath, "EBOUNDS")
        s_table = table.Table.read(filepath, "SPECDATA")
        dref = DMSpecTable.get_ref_vals(filepath)
        return DMSpecTable(e_table, s_table, dref)

    @staticmethod
    def create_from_data(emins, emaxs, refs, data_dict, ref_vals):
        """ Build a DMSpecTable object from energy binning, reference values and spectral values
        """
        e_table = DMSpecTable.make_ebounds_table(emins, emaxs, refs)
        s_table = DMSpecTable.make_spectra_table(len(emins), data_dict)
        return DMSpecTable(e_table, s_table, ref_vals)

    @classmethod
    def create(cls, emin, emax, channels, masses):
        """Create a DM spectrum table.

        Parameters
        ----------

        emin : `~numpy.ndarray`
            Low bin edges.

        emax : `~numpy.ndarray`
            High bin edges.

        channels : list
            List of channel names.

        masses : `~numpy.ndarray`
            Mass points at which to evaluate the spectrum.

        Returns
        -------

        output : `DMSpecTable`
            The newly created object

        """
        ebin_edges = np.concatenate((emin, emax[-1:]))
        evals = np.sqrt(ebin_edges[:-1] * ebin_edges[1:])

        init_params_ann = np.array([REF_SIGV, 100.])
        init_params_dec = np.array([REF_TAU, 100.])
        dmf = DMFitFunction(init_params_ann, jfactor=REF_J, dfactor=REF_D)
        ref_vals = {"REF_J": REF_J,
                    "REF_SIGV": REF_SIGV,
                    "REF_D": REF_D,
                    "REF_TAU": REF_TAU}

        nebins = len(ebin_edges) - 1
        nrow = len(channels) * len(masses)
        dnde = np.ndarray((nrow, nebins))
        flux = np.ndarray((nrow, nebins))
        eflux = np.ndarray((nrow, nebins))
        masses_out = np.ndarray((nrow))
        chan_ids = np.ndarray((nrow), int)

        for i, chan in enumerate(channels):
            ichan = DMFitFunction.channel_rev_map[chan]
            dmf.set_channel(ichan)
            if dmf.decay:
                params = init_params_dec
            else:
                params = init_params_ann            
            s = slice(i * len(masses), (i + 1) * len(masses))
            dnde[s] = dmf.dnde(evals, (params[0], masses)).T
            flux[s] = dmf.flux(emin, emax, (params[0], masses)).T
            eflux[s] = dmf.eflux(emin, emax, (params[0], masses)).T
            chan_ids[s] = ichan
            masses_out[s] = masses

        spec_dict = {"dnde": dnde,
                     "flux": flux,
                     "eflux": eflux,
                     "mass": masses_out,
                     "chan": chan_ids}

        return cls.create_from_data(ebin_edges[0:-1], ebin_edges[1:],
                                    evals, spec_dict, ref_vals)

    @classmethod
    def create_from_config(cls, configfile, channels, masses):
        """ Build a DMSpecTable object from a yaml config file

        Parameters
        ----------

        configfile : str
            Fermipy yaml configuration file

        channels : list
            List of str with the names of the channels to consider

        masses : numpy.array
            DM masses to test, in GeV


        Returns
        -------

        output : `DMSpecTable`
            The newly created object

        """
        config = yaml.safe_load(open(configfile))

        # look for components
        components = config.get('components', [config])
        
        emins = np.array([])
        emaxs = np.array([])

        binsperdec = config['binning'].get('binsperdec', None)

        for comp in components:
            try:
                emin = comp['selection']['emin']
                emax = comp['selection']['emax']
                logemin = np.log10(emin)
                logemax = np.log10(emax)
            except AttributeError:
                logemin = comp['selection']['logemin']
                logemax = comp['selection']['logemax']
                emin = np.power(10., logemin)
                emax = np.power(10., logemax)

            try:
                nebins = comp['binning'].get('enumbins', None)
            except KeyError:
                nebins = np.round(binsperdec * np.log10(emax / emin))

            if nebins is None:
                nebins = np.round(comp['binning']['binsperdec'] * np.log10(emax / emin))

            ebin_edges = np.logspace(logemin, logemax, nebins + 1)
            emins = np.append(emins, ebin_edges[:-1])
            emaxs = np.append(emaxs, ebin_edges[1:])

        return cls.create(emins, emaxs, channels, masses)


    def check_energy_bins(self, ref_spec, tol=1e-3):
        """ Make sure that the energy binning matches the reference spectrum

        Parameters
        ----------

        ref_spec :
            The reference spectrum

        tol : float
            The maximum allowed relative difference in bin edges


        Returns
        -------

        status : bool
            True if the energy binnings match

        """
        emin_local = self._e_table['E_MIN']
        emax_local = self._e_table['E_MAX']
        try:
            if str(emin_local.unit) == 'keV':
                emin_local = emin_local / 1000.
        except KeyError:
            pass
        try:
            if str(emax_local.unit) == 'keV':
                emax_local = emax_local / 1000.
        except KeyError:
            pass

        if len(emin_local) != len(ref_spec.emin):
            return False
        if len(emax_local) != len(ref_spec.emax):
            return False
        if (np.abs(emin_local - ref_spec.emin) > tol * emin_local).any():
            return False
        if (np.abs(emax_local - ref_spec.emax) > tol * emax_local).any():
            return False
        return True

    def convert_castro_data(self, castro_data, channel,
                            norm_type, astro_factor=None):
        """ Convert CastroData object, i.e., Likelihood as a function of
        flux and energy flux, to a DMCastroData object, i.e., Likelihood as
        a function of DM mass and sigmav

        Parameters
        ----------

        castro_data : `CastroData`
            Input data

        channel : int
            Index for the DM interaction channel

        norm_type : str
            Type of spectral normalization to use for the conversion

        astro_factor : dict or float or None
            Nuisance factor used to make the conversion.

            If astro_factor is None, it will use the reference value
            If astro_factor is a float, it will use that
            If astro_factor is a dict, it will use that to create a prior

        Returns
        -------

        output : `DMCastroData`
            The DM-space likelihood curves


        """
        if not self.check_energy_bins(castro_data.refSpec):
            raise ValueError("CastroData energy binning does not match")
        mask = self._s_table['ref_chan'] == channel
        masses = self._s_table[mask]['ref_mass'].data
        spec_vals = self._s_table[mask]['ref_%s' % norm_type].data
        nmass = len(masses)

        is_decay = channel >= 100

        # Get the reference values
        if is_decay:
            ref_inter = self._ref_vals["REF_TAU"]
            ref_norm = self._ref_vals["REF_D"]
            astro_str = 'd_value'
        else:
            ref_inter = self._ref_vals["REF_SIGV"]
            ref_norm = self._ref_vals["REF_J"]
            astro_str = 'j_value'

        astro_prior = None
        if is_null(astro_factor):
            # Just use the reference values
            astro_value = ref_norm
            norm_factor = 1.
        elif isinstance(astro_factor, float):
            # Rescale the normalization values
            astro_value = astro_factor
            norm_factor = ref_norm / astro_factor
        elif isinstance(astro_factor, dict):
            astro_value = astro_factor.get(astro_str)
            norm_factor = ref_norm / astro_value
            astro_factor['scale'] = astro_value
            astro_functype = astro_factor.get('functype', None)
            if is_null(astro_functype):
                astro_prior = None
            else:
                astro_prior = create_prior_functor(astro_factor)
        else:
            sys.stderr.write(
                "Did not recoginize Astro factor %s %s\n" %
                (astro_factor, type(astro_factor)))

        norm_limits = castro_data.getLimits(1e-5)
        # This puts the spectrum in units of the reference spectrum
        # This means that the scan values will be expressed
        # In units of the reference spectra as well        
        spec_vals /= norm_factor
        n_scan_pt = 200

        norm_vals = np.ndarray((nmass, n_scan_pt))
        dll_vals = np.ndarray((nmass, n_scan_pt))
        mle_vals = np.ndarray((nmass))

        mass_mask = np.ones((nmass), bool)
        # for i, mass in enumerate(masses):
        for i in range(nmass):
            max_ratio = 1. / ((spec_vals[i] / norm_limits).max())
            log_max_ratio = np.log10(max_ratio)
            norm_vals[i][0] = 10**(log_max_ratio - 5)
            norm_vals[i][1:] = np.logspace(log_max_ratio - 4, log_max_ratio + 4,
                                           n_scan_pt - 1)
            test_vals = (np.expand_dims(
                spec_vals[i], 1) * (np.expand_dims(norm_vals[i], 1).T))
            dll_vals[i, 0:] = castro_data(test_vals)
            mle_vals[i] = norm_vals[i][dll_vals[i].argmin()]
            mle_ll = dll_vals[i].min()
            dll_vals[i] -= mle_ll

            msk = np.isfinite(dll_vals[i])
            if not msk.any():
                print (
                    "Skipping mass %0.2e for channel %s" %
                    (masses[i], channel))
                mass_mask[i] = False
                continue

            if is_not_null(astro_prior):
                try:
                    lnlfn = castro.LnLFn(norm_vals[i], dll_vals[i], 'dummy')
                    lnlfn_prior = LnLFn_norm_prior(lnlfn, astro_prior)
                    dll_vals[i, 0:] = lnlfn_prior(norm_vals[i])
                except ValueError:
                    print (
                        "Skipping mass %0.2e for channel %s" %
                        (masses[i], channel))
                    mass_mask[i] = False
                    dll_vals[i, 0:] = np.nan * np.ones((n_scan_pt))

        # Here we convert the normalization values to standard units                    
        norm_vals *= (ref_inter)
        dm_castro = DMCastroData(norm_vals[mass_mask], dll_vals[mass_mask],
                                 channel, masses[mass_mask], astro_value,
                                 astro_prior=astro_prior, ref_astro=ref_norm,
                                 ref_inter=ref_inter, decay=is_decay)
        return dm_castro

