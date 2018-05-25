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


class DMCastroData(castro.CastroData_Base):
    """ This class wraps the data needed to make a "Castro" plot,
    namely the log-likelihood as a function of <sigma v> * J (or D) for a
    series of DM masses
    """

    def __init__(self, norm_vals, nll_vals, channel, masses, astro_value,
                 astro_prior=None, prior_applied=True, ref_j=REF_J, ref_sigmav=REF_SIGV, 
                 norm_type='sigmav'):
        """ C'tor

        Parameters
        ----------
        norm_vals : `~numpy.ndarray`
           The normalization values ( n_mass X N array, where N is the
           number of sampled values for each bin )
           Note that these should be the true values, with the 
           reference J-value included, and _NOT_ the values w.r.t. to the 
           reference J-value spectrum.

        nll_vals : `~numpy.ndarray`
           The log-likelihood values ( n_mass X N array, where N is
           the number of sampled values for each bin )

        channel : int or str
           The DM decay or annihilation channel

        masses : '~numpy.ndarray'
           The masses (in GeV)

        astro_value : float
           The J-factor (or D-factor) of the target

        astro_prior : '~fermipy.stats_utils.nuiscance_functor'
           The prior on the J-factor (or D-factor)

        prior_applied : bool
           If true, then the prior has already been applied

        ref_j : float
           Reference value for J-factor

        ref_sigmav : float
           Reference value for sigma
        """
        self._masses = masses
        self._astro_value = astro_value
        self._astro_prior = astro_prior
        self._prior_applied = prior_applied
        self._ref_j = ref_j
        self._ref_sigmav = ref_sigmav

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
    def ref_j(self):
        """ Reference value for J """
        return self._ref_j

    @property
    def ref_sigmav(self):
        """ Reference value for <sigmav> """
        return self._ref_sigmav

    @classmethod
    def create_from_yamlfile(cls, yamlfile, channel, jprior=None):
        """ Create a DMCastroData object from a yaml file 
 
        Parameters
        ----------
        yamlfile : str
           Input file

        channel : str
           The DM decay or annihilation channel

        jprior : str or None
           Type of prior on the J-factor

        Returns
        -------

        output : `DMCastroData`
            Object with the DM-space likelihoods

        """
        data = load_yaml(yamlfile)
        masses = np.array([float(v) for v in data['param']])
        if jprior is not None:
            try:
                astro_value = data['rjvalues'][jprior]
                jsigma = data['jsigma']
            except KeyError:
                astro_value = 1.
                jsigma = 1.
            lnlstr = 'p1lnl'
            prior_dict = dict(functype=jprior,
                              mu=astro_value,
                              sigma=jsigma,
                              j_ref=astro_value)
            prior = create_prior_functor(prior_dict)
        else:
            try:
                astro_value = data['jvalue']
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
                   prior, prior_applied, ref_j=astro_value)


    @classmethod
    def create_from_stack(cls, components, nystep=200, ylims=(1e-30, 1e-20),
                          weights=None, ref_j=REF_J, ref_sigmav=REF_SIGV):
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

        ref_j : float
            Refernce J-factor value

        ref_sigmav : float
            Refernec <sigmav> value

        Returns
        -------

        output : `DMCastroData`
            Object with the DM-space likelihoods

        """
        if not components:
            return None
        shape = (components[0].nx, nystep)
        norm_vals, nll_vals = castro.CastroData_Base.stack_nll(shape, components,
                                                               ylims, weights)
        return cls(norm_vals, nll_vals, components[0].channel, components[0].masses,
                   astro_value=None, ref_j=ref_j, ref_sigmav=ref_sigmav, 
                   norm_type='norm')

    @classmethod
    def create_from_tables(cls, tab_s, tab_m, norm_type):
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

        Returns
        -------

        output : `DMCastroData`
            Object with the DM-space likelihoods

        """
        ref_j = np.squeeze(np.array(tab_s['ref_J']))
        ref_sigmav = np.squeeze(np.array(tab_s['ref_sigmav']))

        if norm_type in ['sigmav']:
            norm_vals = np.squeeze(np.array(tab_s['norm_scan']*ref_sigmav))
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
                              sigma=prior_sigma,
                              j_ref=astro_value)
            prior = create_prior_functor(prior_dict)
        else:
            prior = None
            prior_applied = True

        return cls(norm_vals, nll_vals, channel,
                   masses, astro_value, prior, prior_applied,
                   ref_j=ref_j, ref_sigmav=ref_sigmav)

    @classmethod
    def create_from_fitsfile(cls, filepath, channel, norm_type='sigmav'):
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

        Returns
        -------

        output : `DMCastroData`
            Object with the DM-space likelihoods

        """
        tab_s = Table.read(filepath, hdu=channel)
        tab_m = Table.read(filepath, hdu='masses')
        return cls.create_from_tables(tab_s, tab_m, norm_type=norm_type)

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

    def build_scandata_table(self, norm_type='sigmav'):
        """Build a FITS table with likelihood scan data

        Parameters
        ----------
        
        norm_type : str
            Type of normalization to use.  Valid options are:

            * norm : Self-normalized
            * sigmav : Reference values of sigmav and J

        Returns
        -------

        table : `astropy.table.Table`
            The table has these columns

        astro_value : float
            The astrophysical J-factor for this target
        ref_J : float
            The reference J-factor used to build `DMSpecTable`
        ref_sigmav : float
            The reference <sigmav> used to build `DMSpecTable`

        norm_scan : array
            The test values of <sigmav>
        dloglike_scan : array
            The corresponding values of the negative log-likelihood

        """
        shape = self._norm_vals.shape
        #dtype = 'f%i'%self._norm_vals.size

        col_normv = Column(name="norm_scan", dtype=float,
                           shape=shape)
        col_dll = Column(name="dloglike_scan", dtype=float,
                         shape=shape)

        col_astro_val = Column(name="astro_value", dtype=float)
        col_ref_j = Column(name="ref_J", dtype=float)
        col_ref_sigmav = Column(name="ref_sigmav", dtype=float)

        collist = [col_normv, col_dll, col_astro_val, col_ref_j,
                   col_ref_sigmav]

        if norm_type in ['sigmav']:
            norm_vals = self._norm_vals / self.ref_sigmav
        elif norm_type in ['norm']:
            norm_vals = self._norm_vals
        else:
            raise ValueError('Unrecognized normalization type: %s' % norm_type)

        valdict = {"norm_scan": norm_vals,
                   "dloglike_scan": -1 * self._nll_vals,
                   "astro_value": self.astro_value,
                   "ref_J": self.ref_j,
                   "ref_sigmav": self.ref_sigmav}

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
        ref_J : float
            The reference J-factor used to build `DMSpecTable`
        ref_sigmav : float
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
        col_astro_val = Column(name="astro_value", dtype=float)
        col_ref_j = Column(name="ref_J", dtype=float)
        col_ref_sigmav = Column(name="ref_sigmav", dtype=float)
        collist = [col_astro_val, col_ref_j, col_ref_sigmav]
        valdict = {"astro_value": self.astro_value,
                   "ref_J": self.ref_j,
                   "ref_sigmav": self.ref_sigmav}

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
        self._channel_map, self._channel_names = DMSpecTable.make_channel_map(
            s_table)

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
        dref = {"ref_sigv": hin["ref_sigv"],
                "ref_J": hin["ref_J"]}
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

        init_params = np.array([REF_SIGV, 100.])
        dmf = DMFitFunction(init_params, jfactor=REF_J)

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
            s = slice(i * len(masses), (i + 1) * len(masses))
            dnde[s] = dmf.dnde(evals, (init_params[0], masses)).T
            flux[s] = dmf.flux(emin, emax, (init_params[0], masses)).T
            eflux[s] = dmf.eflux(emin, emax, (init_params[0], masses)).T
            chan_ids[s] = ichan
            masses_out[s] = masses

        spec_dict = {"dnde": dnde,
                     "flux": flux,
                     "eflux": eflux,
                     "mass": masses_out,
                     "chan": chan_ids}

        ref_vals = {"ref_J": REF_J,
                    "ref_sigv": REF_SIGV}

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

        emin = config['selection']['emin']
        emax = config['selection']['emax']
        log_emin = np.log10(emin)
        log_emax = np.log10(emax)
        ndec = log_emax - log_emin
        binsperdec = config['binning']['binsperdec']
        nebins = int(np.ceil(binsperdec * ndec))
        ebin_edges = np.logspace(log_emin, log_emax, nebins + 1)
        return cls.create(ebin_edges[:-1], ebin_edges[1:], channels, masses)

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
                            norm_type, jfactor=None):
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

        jfactor : dict or float or None
            J-factor used to make the conversion.

            If jfactor is None, it will use the reference value
            If jfactor is a float, it will use that
            If jfactor is a dict, it will use that to create a prior

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

        # Get the reference values
        ref_sigmav = self._ref_vals["ref_sigv"]
        ref_norm = self._ref_vals["ref_J"]

        j_prior = None
        if is_null(jfactor):
            # Just use the reference values
            j_value = ref_norm
            norm_factor = 1.
        elif isinstance(jfactor, float):
            # Rescale the normalization values
            j_value = jfactor
            norm_factor = ref_norm / jfactor
        elif isinstance(jfactor, dict):
            j_value = jfactor.get('j_value')
            norm_factor = ref_norm / j_value
            jfactor['scale'] = j_value
            j_functype = jfactor.get('functype', None)
            if is_null(j_functype):
                j_prior = None
            else:
                j_prior = create_prior_functor(jfactor)
        else:
            sys.stderr.write(
                "Did not recoginize J factor %s %s\n" %
                (jfactor, type(jfactor)))

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

            if is_not_null(j_prior):
                try:
                    lnlfn = castro.LnLFn(norm_vals[i], dll_vals[i], 'dummy')
                    lnlfn_prior = LnLFn_norm_prior(lnlfn, j_prior)
                    dll_vals[i, 0:] = lnlfn_prior(norm_vals[i])
                except ValueError:
                    print (
                        "Skipping mass %0.2e for channel %s" %
                        (masses[i], channel))
                    mass_mask[i] = False
                    dll_vals[i, 0:] = np.nan * np.ones((n_scan_pt))

        # Here we convert the normalization values to standard units                    
        norm_vals *= (ref_sigmav)
        dm_castro = DMCastroData(norm_vals[mass_mask], dll_vals[mass_mask],
                                 channel, masses[mass_mask], j_value, j_prior,
                                 ref_j=ref_norm, ref_sigmav=ref_sigmav)
        return dm_castro

    def convert_tscube(self, tscube, channel, norm_type):
        """ Convert TSCube object, i.e., Likelihood as a function of
        flux and energy flux for a set of pixels to a set of DMCastroData objects,
        i.e., Likelihood as a function of DM mass and sigmav for those pixels.

        Parameters
        ----------

        tscube :
            Input data

        channel : int
            Index for the DM interaction channel

        norm_type : str
            Type of spectral normalization to use for the conversion

        Returns
        -------

        dm_table : `astropy.table.Table`
            Table with the log-likelihood data.
            One row per pixel.

        mass_table : `astropy.table.Table`
            Table with the DM masses used

        dm_ts_cube : `skymap.Map`
            Maps of the TS per pixel in each energy bin

        dm_ul_cube : `skymap.Map`
            Maps of the upper limit on <sigma v> per pixel in each energy bin.

        dm_mle_cube : `skymap.Map`
            Maps of the maximum likelihood estimate per pixel in each energy bin.

        """
        if not self.check_energy_bins(tscube.refSpec):
            raise ValueError("TSCube energy binning does not match")

        wcs = tscube.tscube.wcs
        mask = self._s_table['ref_chan'] == channel
        masses = self._s_table[mask]['ref_mass'].data
        spec_vals = self._s_table[mask]['ref_%s' % norm_type].data
        nmasses = len(masses)

        # Get the reference values
        ref_sigmav = self._ref_vals["ref_sigv"]
        ref_norm = self._ref_vals["ref_J"]
        j_ref = ref_norm

        nvals = tscube.nvals
        cube_shape = (
            len(masses),
            tscube.tscube.counts.shape[1],
            tscube.tscube.counts.shape[2])

        dm_ts_cube = np.zeros(cube_shape)
        dm_ul_cube = np.zeros(cube_shape)
        dm_mle_cube = np.zeros(cube_shape)

        dm_table = None
        mass_table = None

        prog = int(nvals / 20)
        for ipix in xrange(nvals):
            xpix, ypix = tscube.tsmap.ipix_to_xypix(ipix, colwise=True)
            if ipix % prog == 0:
                sys.stdout.write('.')
                sys.stdout.flush()

            castro_data = tscube.castroData_from_ipix(ipix)
            norm_limits = castro_data.getLimits(1e-5)

            n_scan_pt = 100

            norm_vals = np.ndarray((nmasses, n_scan_pt))
            dll_vals = np.ndarray((nmasses, n_scan_pt))
            mle_vals = np.ndarray((nmasses))

            # for i, m in enumerate(masses):
            for i in range(nmasses):
                max_ratio = 1. / ((spec_vals[i] / norm_limits).max())
                log_max_ratio = np.log10(max_ratio)
                norm_vals[i][0] = 0.
                norm_vals[i][1:] = np.logspace(log_max_ratio - 2,
                                               log_max_ratio + 2, n_scan_pt - 1)
                test_vals = (np.expand_dims(spec_vals[i], 1) *
                             (np.expand_dims(norm_vals[i], 1).T))
                dll_vals[i, 0:] = castro_data(test_vals)
                mle_vals[i] = norm_vals[i][dll_vals[i].argmin()]
                mle_ll = dll_vals[i].min()
                dll_vals[i] -= mle_ll

            norm_vals *= (ref_sigmav)
            dm_castro = DMCastroData(norm_vals, dll_vals,
                                     channel, masses, j_ref)
            tab = dm_castro.build_scandata_table()
            if dm_table is None:
                dm_table = tab
            else:
                dm_table.add_row(tab[0])
            if mass_table is None:
                mass_table = dm_castro.build_mass_table()

            dm_ts_cube[0:, xpix, ypix] = dm_castro.ts_vals()
            dm_ul_cube[0:, xpix, ypix] = dm_castro.getLimits(0.05)
            dm_mle_cube[0:, xpix, ypix] = dm_castro.mles()

        dm_ts_cube = skymap.Map(dm_ts_cube, wcs)
        dm_ul_cube = skymap.Map(dm_ul_cube, wcs)
        dm_mle_cube = skymap.Map(dm_mle_cube, wcs)

        return (dm_table, mass_table, dm_ts_cube, dm_ul_cube, dm_mle_cube)
