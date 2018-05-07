#!/usr/bin/env python
#

"""
Interface to Dark Matter spectra
"""
from __future__ import absolute_import, division, print_function

import sys
import os

import yaml
import numpy as np

from astropy.table import Table, Column

from fermipy import castro
from fermipy import fits_utils
from fermipy import stats_utils
from fermipy import skymap
from fermipy.castro import CastroData

from fermipy.utils import load_yaml
from fermipy.jobs.utils import is_null, is_not_null
from fermipy.jobs.link import Link
from fermipy.jobs.scatter_gather import ScatterGather
from fermipy.jobs.slac_impl import make_nfs_path

from fermipy.spectrum import DMFitFunction

from dmpipe.name_policy import NameFactory
from dmpipe import defaults

NAME_FACTORY = NameFactory(basedir='.')

REF_J = 1.0e19
REF_SIGV = 1.0e-26


class DMCastroData(castro.CastroData_Base):
    """ This class wraps the data needed to make a "Castro" plot,
    namely the log-likelihood as a function of <sigma v> * J (or D) for a
    series of DM masses
    """

    def __init__(self, norm_vals, nll_vals, norm_value, channel, masses, astro_value,
                 astro_prior=None, prior_applied=True, ref_j=REF_J, ref_sigmav=REF_SIGV):
        """ C'tor

        Parameters
        ----------
        norm_vals : `~numpy.ndarray`
           The normalization values ( n_mass X N array, where N is the
           number of sampled values for each bin )

        nll_vals : `~numpy.ndarray`
           The log-likelihood values ( n_mass X N array, where N is
           the number of sampled values for each bin )

        norm_value : float
           The normalization value

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
        self._norm_value = norm_value
        if self._norm_value is None:
            self._norm_value = 1.0
        self._astro_value = astro_value
        self._astro_prior = astro_prior
        self._prior_applied = prior_applied
        self._ref_j = ref_j
        self._ref_sigmav = ref_sigmav

        if self._astro_value is not None:
            test_norm = self._ref_j / self._astro_value
        else:
            test_norm = self._norm_value
        if np.fabs(test_norm - self._norm_value) > 0.01:
            sys.stderr.write(
                "WARNING, normalization value does not match expected value: %.2E %.2E\n" %
                (self._norm_value, test_norm))

        if isinstance(channel, str):
            self._channel = DMFitFunction.channel_rev_map[channel]
        else:
            self._channel = channel
        super(DMCastroData, self).__init__(norm_vals,
                                           nll_vals,
                                           norm_type="sigmav")

    @property
    def n_masses(self):
        """ The number of masses tested """
        return self._nx

    @property
    def masses(self):
        """ The masses tested (in GeV) """
        return self._masses

    @property
    def norm_value(self):
        """ The global normalization value

        The array off normalization values should be multiplied by this factor """
        return self._norm_value

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

    @staticmethod
    def create_from_stack(components, nystep=200, ylims=(1e-30, 1e-20),
                          weights=None, ref_j=REF_J, ref_sigmav=REF_SIGV):
        """ Create a DMCastroData object by stacking a series of DMCastroData objects
        """
        if not components:
            return None
        shape = (components[0].nx, nystep)
        norm_vals, nll_vals = castro.CastroData_Base.stack_nll(shape, components,
                                                               ylims, weights)
        return DMCastroData(norm_vals, nll_vals, 1.0,
                            components[0].channel, components[0].masses,
                            astro_value=None, ref_j=ref_j, ref_sigmav=ref_sigmav)

    @staticmethod
    def create_from_tables(tab_s, tab_m):
        """ Create a DMCastroData object from likelihood scan and mass tables
        """
        ref_j = np.squeeze(np.array(tab_s['REF_J']))
        ref_sigmav = np.squeeze(np.array(tab_s['REF_SIGMAV']))

        norm = np.squeeze(np.array(tab_s['NORM']))
        norm_vals = np.squeeze(np.array(tab_s['NORM_SCAN'])) * norm
        nll_vals = -np.squeeze(np.array(tab_s['DLOGLIKE_SCAN']))

        masses = np.squeeze(np.array(tab_m['MASSES']))
        channel = np.squeeze(np.array(tab_m['CHANNEL']))

        astro_value = np.squeeze(np.array(tab_s['ASTRO_VALUE']))
        try:
            astro_priortype = tab_s['PRIOR_TYPE']
        except KeyError:
            astro_priortype = None

        if astro_priortype is not None:
            prior_mean = np.squeeze(np.array(tab_s['PRIOR_MEAN']))
            prior_sigma = np.squeeze(np.array(tab_s['PRIOR_SIGMA']))
            prior_applied = np.squeeze(np.array(tab_s['PRIOR_APPLIED']))
            prior_dict = dict(functype=astro_priortype,
                              mu=prior_mean,
                              sigma=prior_sigma,
                              j_ref=astro_value)
            prior = stats_utils.create_prior_functor(prior_dict)
        else:
            prior = None
            prior_applied = True

        return DMCastroData(norm_vals, nll_vals, norm, channel,
                            masses, astro_value, prior, prior_applied,
                            ref_j=ref_j, ref_sigmav=ref_sigmav)

    @staticmethod
    def create_from_fitsfile(filepath, channel):
        """ Create a DMCastroData object likelihood scan and mass tables in FITS file """
        tab_s = Table.read(filepath, hdu=channel)
        tab_m = Table.read(filepath, hdu='MASSES')
        return DMCastroData.create_from_tables(tab_s, tab_m)

    def build_lnl_fn(self, normv, nllv):
        """ Build a function to return the likelihood value arrays of
        normalization and likelihood values
        """
        lnlfn = castro.LnLFn(normv, nllv, self._norm_type)
        if self._astro_prior is None:
            return lnlfn
        if self._prior_applied:
            return lnlfn
        return stats_utils.LnLFn_norm_prior(lnlfn, self._astro_prior)

    def build_scandata_table(self):
        """Build a FITS table with likelihood scan data
        """
        shape = self._norm_vals.shape
        #dtype = 'f%i'%self._norm_vals.size

        col_norm = Column(name="NORM", dtype=float)
        col_normv = Column(name="NORM_SCAN", dtype=float,
                           shape=shape)
        col_dll = Column(name="DLOGLIKE_SCAN", dtype=float,
                         shape=shape)

        col_astro_val = Column(name="ASTRO_VALUE", dtype=float)
        col_ref_j = Column(name="REF_J", dtype=float)
        col_ref_sigmav = Column(name="REF_SIGMAV", dtype=float)

        collist = [
            col_norm,
            col_normv,
            col_dll,
            col_astro_val,
            col_ref_j,
            col_ref_sigmav]
        valdict = {"NORM": self._norm_value,
                   "NORM_SCAN": self._norm_vals * self._norm_value,
                   "DLOGLIKE_SCAN": -1 * self._nll_vals,
                   "ASTRO_VALUE": self.astro_value,
                   "REF_J": self.ref_j,
                   "REF_SIGMAV": self.ref_sigmav}

        if self._astro_prior is not None:
            col_prior_type = Column(name="PRIOR_TYPE", dtype="S16")
            col_prior_mean = Column(name="PRIOR_MEAN", dtype=float)
            col_prior_sigma = Column(name="PRIOR_SIGMA", dtype=float)
            col_prior_applied = Column(name="PRIOR_APPLIED", dtype=bool)
            collist += [col_prior_type, col_prior_mean,
                        col_prior_sigma, col_prior_applied]
            valdict["PRIOR_TYPE"] = self.prior_type
            valdict["PRIOR_MEAN"] = self.prior_mean
            valdict["PRIOR_SIGMA"] = self.prior_sigma
            valdict["PRIOR_APPLIED"] = self.prior_applied

        tab = Table(data=collist)
        tab.add_row(valdict)
        return tab

    def build_mass_table(self):
        """Build a FITS table with mass values
        """
        col_masses = Column(name="MASSES", dtype=float,
                            shape=self._masses.shape)

        col_channel = Column(name="CHANNEL", dtype=int)
        tab = Table(data=[col_masses, col_channel])
        tab.add_row({"MASSES": self._masses,
                     "CHANNEL": self._channel})
        return tab

    def build_limits_table(self, limit_dict):
        """Build a FITS table with limits data
        """
        col_norm = Column(name="NORM", dtype=float)
        col_astro_val = Column(name="ASTRO_VALUE", dtype=float)
        col_ref_j = Column(name="REF_J", dtype=float)
        col_ref_sigmav = Column(name="REF_SIGMAV", dtype=float)
        collist = [col_norm, col_astro_val, col_ref_j, col_ref_sigmav]
        valdict = {"NORM": self._norm_value,
                   "ASTRO_VALUE": self.astro_value,
                   "REF_J": self.ref_j,
                   "REF_SIGMAV": self.ref_sigmav}

        for k, v in limit_dict.items():
            collist.append(Column(name=k, dtype=float, shape=v.shape))
            valdict[k] = v

        if self._astro_prior is not None:
            col_prior_type = Column(name="PRIOR_TYPE", dtype="S16")
            col_prior_mean = Column(name="PRIOR_MEAN", dtype=float)
            col_prior_sigma = Column(name="PRIOR_SIGMA", dtype=float)
            col_prior_applied = Column(name="PRIOR_APPLIED", dtype=bool)
            collist += [col_prior_type, col_prior_mean,
                        col_prior_sigma, col_prior_applied]
            valdict["PRIOR_TYPE"] = self.prior_type
            valdict["PRIOR_MEAN"] = self.prior_mean
            valdict["PRIOR_SIGMA"] = self.prior_sigma
            valdict["PRIOR_APPLIED"] = self.prior_applied

        tab = Table(data=collist)
        tab.add_row(valdict)
        return tab



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
        """ Construct the spectral values table
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
        """ Return the energy binning table """
        return self._e_table

    @property
    def spectra_table(self):
        """ Return the spectral values table """
        return self._s_table

    @property
    def ref_vals(self):
        """ Return the spectral reference values """
        return self._ref_vals

    @property
    def channel_map(self):
        """ Return the channel to index mapping """
        return self._channel_map

    @property
    def channel_names(self):
        """ Return the channel to index mapping """
        return self._channel_names

    def spectrum(self, chan, mass, spec_type):
        """ Return the spectrum for a particular channel an mass
        """
        mask = (self._s_table["ref_chan"] == chan) & (
            np.abs(self._s_table["ref_mass"] - mass) < 1e-9)
        spec_vals = self._s_table[mask]["ref_%s" % spec_type].data
        return spec_vals

    def masses(self, chan):
        """ Return the array of masses for a given channel
        """
        mask = (self._s_table["ref_chan"] == chan)
        return self._s_table[mask]["ref_mass"]

    def ebin_edges(self):
        """ Return an array with the energy bin edges
        """
        return np.hstack([self._e_table["E_MIN"].data,
                          self._e_table["E_MAX"].data])

    def ebin_refs(self):
        """ Return an array with the energy bin reference energies
        """
        return self._e_table["E_REF"].data

    def write_fits(self, filepath, clobber=False):
        """ Write this object to a FITS file
        """
        fits_utils.write_tables_to_fits(filepath, [self._e_table, self._s_table],
                                        clobber=clobber,
                                        namelist=["EBOUNDS", "SPECDATA"],
                                        cardslist=[{}, self._ref_vals])

    @staticmethod
    def get_ref_vals(filepath):
        """ Extract the reference values from a FITS header """
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
        """
        ebin_edges = np.concatenate((emin, emax[-1:]))
        evals = np.sqrt(ebin_edges[:-1] * ebin_edges[1:])
        ichans = DMFitFunction.channel_name_mapping.keys()

        init_params = np.array([REF_SIGV, 100.])
        dmf = DMFitFunction(init_params, jfactor=REF_J)

        nebins = len(ebin_edges) - 1
        nrow = len(ichans) * len(masses)
        dnde = np.ndarray((nrow, nebins))
        flux = np.ndarray((nrow, nebins))
        eflux = np.ndarray((nrow, nebins))
        masses_out = np.ndarray((nrow))
        channels = np.ndarray((nrow), int)

        for i, ichan in enumerate(ichans):
            dmf.set_channel(ichan)
            s = slice(i * len(masses), (i + 1) * len(masses))
            dnde[s] = dmf.dnde(evals, (init_params[0], masses)).T
            flux[s] = dmf.flux(emin, emax, (init_params[0], masses)).T
            eflux[s] = dmf.eflux(emin, emax, (init_params[0], masses)).T
            channels[s] = ichan
            masses_out[s] = masses

        spec_dict = {"dnde": dnde,
                     "flux": flux,
                     "eflux": eflux,
                     "mass": masses_out,
                     "chan": channels}

        ref_vals = {"ref_J": REF_J,
                    "ref_sigv": REF_SIGV}

        return cls.create_from_data(ebin_edges[0:-1], ebin_edges[1:],
                                    evals, spec_dict, ref_vals)

    @classmethod
    def create_from_config(cls, configfile, channels, masses):
        """ Build a DMSpecTable object from a yaml config file
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
            jfactor['j_ref'] = ref_norm
            j_value = jfactor.get('j_value')
            norm_factor = ref_norm / j_value
            j_functype = jfactor.get('functype', None)
            if is_null(j_functype):
                j_prior = None
            else:
                j_prior = stats_utils.create_prior_functor(jfactor)
        else:
            sys.stderr.write(
                "Did not recoginize J factor %s %s\n" %
                (jfactor, type(jfactor)))

        norm_limits = castro_data.getLimits(1e-5)
        spec_vals *= norm_factor
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
                    lnlfn_prior = stats_utils.LnLFn_norm_prior(lnlfn, j_prior)
                    dll_vals[i, 0:] = lnlfn_prior(norm_vals[i])
                except ValueError:
                    print (
                        "Skipping mass %0.2e for channel %s" %
                        (masses[i], channel))
                    mass_mask[i] = False
                    dll_vals[i, 0:] = np.nan * np.ones((n_scan_pt))

        norm_vals *= (ref_sigmav)
        dm_castro = DMCastroData(norm_vals[mass_mask], dll_vals[mass_mask], norm_factor,
                                 channel, masses[mass_mask], j_value, j_prior,
                                 ref_j=ref_norm, ref_sigmav=ref_sigmav)
        return dm_castro

    def convert_tscube(self, tscube, channel, norm_type):
        """ Convert TSCube object, i.e., Likelihood as a function of
        flux and energy flux for a set of pixels to a set of DMCastroData objects,
        i.e., Likelihood as a function of DM mass and sigmav for those pixels.
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
        norm_factor = 1.

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
            dm_castro = DMCastroData(
                norm_vals,
                dll_vals,
                norm_factor,
                channel,
                masses,
                j_ref)
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

        return dm_table, mass_table, dm_ts_cube, dm_ul_cube, dm_mle_cube


class ConvertCastro(Link):
    """Small class to convert CastroData to a DMCastroData"""

    appname = 'dmpipe-convert-castro'
    linkname_default = 'convert-castro'
    usage = '%s [options]' % (appname)
    description = "Convert SED to DMCastroData"

    default_options = dict(specfile=defaults.common['specfile'],
                           sed_file=defaults.common['sed_file'],
                           j_value_file=defaults.common['j_value_file'],
                           jprior=defaults.common['jprior'],
                           outfile=defaults.generic['outfile'],
                           limitfile=defaults.generic['limitfile'],
                           # Note that this defaults to -1
                           nsims=defaults.common['nsims'],
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'],
                           clobber=defaults.common['clobber'])

    @staticmethod
    def convert_sed_to_dm(spec_table, sed, channels, norm_type, j_val):
        """ Convert an SED file to a DMCastroData object """
        c_list = []
        t_list = []
        n_list = []

        mass_table = None

        for chan in channels:
            chan_idx = DMFitFunction.channel_rev_map[chan]
            # try:
            dm_castro = spec_table.convert_castro_data(
                sed, chan_idx, norm_type, j_val)
            tab_castro = dm_castro.build_scandata_table()

            if mass_table is None:
                mass_table = dm_castro.build_mass_table()
            # except IndexError, msg:
            #    raise IndexError("Skipping channel %s" % msg)
            #    continue
            c_list.append(dm_castro)
            t_list.append(tab_castro)
            n_list.append(chan)

        t_list.append(mass_table)
        n_list.append("MASSES")
        return c_list, t_list, n_list

    @staticmethod
    def extract_dm_limits(dm_castro_list, channels, alphas, mass_table):
        """ Convert an SED file to a DMCastroData object """
        l_list = []
        t_list = []
        n_list = []

        for castro_data, chan in zip(dm_castro_list, channels):
            norm = castro_data.norm_value
            mles = norm * castro_data.mles()
            limit_dict = dict(MLES=mles)
            for alpha in alphas:
                limits = norm * castro_data.getLimits(alpha)
                limit_dict['UL_%.02f' % alpha] = limits

            tab_limits = castro_data.build_limits_table(limit_dict)
            l_list.append(limit_dict)
            t_list.append(tab_limits)
            n_list.append(chan)

        t_list.append(mass_table)
        n_list.append("MASSES")
        return l_list, t_list, n_list

    @staticmethod
    def convert_sed(spec_table, sed_file, norm_type, channels,
                    j_factor, outfile, limitfile, clobber):
        """Convert a single SED to DM"""
        sed = CastroData.create_from_sedfile(sed_file, norm_type)
        c_list, t_list, n_list = ConvertCastro.convert_sed_to_dm(
            spec_table, sed, channels, norm_type, j_factor)

        if is_not_null(outfile):
            fits_utils.write_tables_to_fits(
                outfile, t_list, clobber=clobber, namelist=n_list)

        if is_not_null(limitfile):
            mass_table = t_list[-1]
            c_list_lim, t_list_lim, n_list_lim = ConvertCastro.extract_dm_limits(
                c_list, channels, [0.68, 0.95], mass_table)
            fits_utils.write_tables_to_fits(limitfile, t_list_lim,
                                            clobber=clobber, namelist=n_list_lim)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        norm_type = 'eflux'
        channels = None

        spec_table = DMSpecTable.create_from_fits(args.specfile)
        profile = load_yaml(args.j_value_file)

        if channels is None:
            channels = spec_table.channel_names

        j_value = profile.get('j_integ')
        j_sigma = profile.get('j_sigma', None)
        if is_null(args.jprior) or is_null(j_sigma) or j_sigma == 0.0:
            j_factor = j_value
        else:
            j_factor = dict(functype=args.jprior,
                            j_value=j_value,
                            mu=j_value, sigma=j_sigma)

        if args.nsims < 0:
            seedlist = [None]
        else:
            seedlist = range(args.seed, args.seed + args.nsims)

        for seed in seedlist:
            sedfile = args.sed_file
            outfile = args.outfile
            limitfile = args.limitfile
            if seed is not None:
                sedfile = sedfile.replace('_SEED.fits', '_%06i.fits' % seed)
                if is_not_null(outfile):
                    outfile = outfile.replace(
                        '_SEED.fits',
                        '_%06i.fits' %
                        seed)
                if is_not_null(limitfile):
                    limitfile = limitfile.replace(
                        '_SEED.fits', '_%06i.fits' % seed)

            self.convert_sed(spec_table, sedfile, norm_type,
                             channels, j_factor, outfile,
                             limitfile, args.clobber)


class SpecTable(Link):
    """ Version of the DM spectral tables in tabular form
    """
    appname = 'dmpipe-spec-table'
    linkname_default = 'spec-table'
    usage = '%s [options]' % (appname)
    description = "Build a table with the spectra for DM signals"

    default_options = dict(ttype=defaults.common['ttype'],
                           config=defaults.common['config'],
                           specconfig=defaults.common['specconfig'],
                           specfile=defaults.common['specfile'],
                           dry_run=defaults.common['dry_run'],
                           clobber=defaults.common['clobber'])


    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if args.ttype is not None:
            name_keys = dict(target_type=args.ttype,
                             fullpath=True)
            config_file = NAME_FACTORY.ttypeconfig(**name_keys)
            spec_config = NAME_FACTORY.specconfig(**name_keys)
            spec_file = NAME_FACTORY.specfile(**name_keys)
        else:
            config_file = None
            spec_config = None
            spec_file = None

        if is_not_null(args.config):
            config_file = args.config_file
        if is_not_null(args.specconfig):
            spec_config = args.specconfig
        if is_not_null(args.specfile):
            spec_file = args.specfile

        if config_file is None:
            sys.stderr.write('No input configuration file is specified')
            return -1

        if spec_config is None:
            sys.stderr.write('No input spectra configurate file is specified')
            return -1

        if spec_file is None:
            sys.stderr.write('No output spectra file is specified')
            return -1

        spec_config = load_yaml(spec_config)
        channels = spec_config['channels']
        masses = np.logspace(np.log10(spec_config['masses']['mass_min']),
                             np.log10(spec_config['masses']['mass_max']),
                             spec_config['masses']['mass_nstep'])

        dm_spec_table = DMSpecTable.create_from_config(
            config_file, channels, masses)
        dm_spec_table.write_fits(spec_file, args.clobber)
        return 0


class StackLikelihood(Link):
    """Small class to convert stack DMCastroData """
    appname = 'dmpipe-stack-likelihood'
    linkname_default = 'stack-likelihood'
    usage = '%s [options]' % (appname)
    description = "Stack the likelihood from a set of targets"

    default_options = dict(ttype=defaults.common['ttype'],
                           specconfig=defaults.common['specconfig'],
                           rosterlist=defaults.common['rosterlist'],
                           jprior=defaults.common['jprior'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'],
                           clobber=defaults.common['clobber'])

    @staticmethod
    def stack_roster(rost, ttype,
                     channels, jprior_key, sim, seed):
        """ Stack all of the DMCastroData in a roster
        """
        component_dict = {}
        out_dict = {}
        for chan in channels:
            component_dict[chan] = []

        for target_key in rost:
            tokens = target_key.split(':')
            name_keys = dict(target_type=ttype,
                             target_name=tokens[0],
                             profile=tokens[1],
                             fullpath=True,
                             sim_name=sim,
                             seed="%06i" % seed,
                             jprior=jprior_key)

            if is_not_null(sim):
                dmlike_path = NAME_FACTORY.sim_dmlikefile(**name_keys)
            else:
                dmlike_path = NAME_FACTORY.dmlikefile(**name_keys)

            tab_m = Table.read(dmlike_path, hdu="MASSES")

            for chan in channels:
                try:
                    tab_s = Table.read(dmlike_path, hdu=chan)
                except KeyError:
                    continue
                dm_castro = DMCastroData.create_from_tables(tab_s, tab_m)
                component_dict[chan].append(dm_castro)

        for chan, comps in component_dict.items():
            if not comps:
                continue
            stacked = DMCastroData.create_from_stack(comps, ref_j=REF_J,
                                                     ref_sigmav=REF_SIGV)
            out_dict[chan] = stacked

        return out_dict

    @staticmethod
    def write_stacked(ttype, roster_name, stacked_dict,
                      jprior_key, sim, seed, clobber):
        """ Write the stacked DMCastroData object to a FITS file
        """
        name_keys = dict(target_type=ttype,
                         target_name="stacked",
                         fullpath=True,
                         roster_name=roster_name,
                         sim_name=sim,
                         seed="%06i" % seed,
                         jprior=jprior_key)

        if is_not_null(sim):
            outdir = NAME_FACTORY.sim_targetdir(**name_keys)
            outpath = NAME_FACTORY.sim_resultsfile(**name_keys)
        else:
            outdir = NAME_FACTORY.targetdir(**name_keys)
            outpath = NAME_FACTORY.resultsfile(**name_keys)

        try:
            os.makedirs(outdir)
        except OSError:
            pass

        limitfile = outpath.replace('results', 'limits')
        print("Writing stacked results %s" % outpath)
        channels = stacked_dict.keys()
        t_list = []
        n_list = []
        lim_list = []
        lim_table_list = []
        mass_table = None
        alphas = [0.68, 0.95]
        for chan in channels:
            stacked = stacked_dict[chan]
            norm = stacked.norm_value
            mles = norm * stacked.mles()
            limit_dict = dict(MLES=mles)
            for alpha in alphas:
                limits = norm * stacked.getLimits(alpha)
                limit_dict['UL_%.02f' % alpha] = limits
            tab_limits = stacked.build_limits_table(limit_dict)
            if mass_table is None:
                mass_table = stacked.build_mass_table()
            t_list.append(stacked.build_scandata_table())
            n_list.append(chan)
            lim_list.append(limit_dict)
            lim_table_list.append(tab_limits)

        t_list.append(mass_table)
        lim_table_list.append(mass_table)
        n_list.append("MASSES")
        fits_utils.write_tables_to_fits(outpath, t_list,
                                        clobber=clobber, namelist=n_list)
        fits_utils.write_tables_to_fits(limitfile, lim_table_list,
                                        clobber=clobber, namelist=n_list)

    @staticmethod
    def stack_rosters(roster_dict, ttype, channels,
                      jprior_key, sim, seed, clobber):
        """ Stack all of the DMCastroData in a dictionary of rosters
        """
        for roster_name, rost in roster_dict.items():
            stacked_dict = StackLikelihood.stack_roster(rost, ttype,
                                                        channels, jprior_key, sim, seed)
            StackLikelihood.write_stacked(ttype, roster_name, stacked_dict,
                                          jprior_key, sim, seed, clobber)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if args.ttype is None:
            raise RuntimeError('Target type must be specified')

        name_keys = dict(target_type=args.ttype,
                         rosterlist='roster_list.yaml',
                         sim_name=args.sim,
                         fullpath=True)

        spec_config = NAME_FACTORY.specconfig(**name_keys)
        if is_not_null(args.specconfig):
            spec_config = args.specconfig

        spec_config = load_yaml(spec_config)
        channels = spec_config['channels']

        if is_not_null(args.sim):
            roster_file = NAME_FACTORY.sim_rosterfile(**name_keys)
            sim_name = args.sim
            is_sim = True
        else:
            roster_file = NAME_FACTORY.rosterfile(**name_keys)
            is_sim = False
            sim_name = None

        if is_not_null(args.rosterlist):
            roster_file = args.rosterlist

        roster_dict = load_yaml(roster_file)

        if is_sim:
            seedlist = range(args.seed, args.seed + args.nsims)
        else:
            seedlist = [0]

        jprior = args.jprior
        if is_null(jprior):
            jprior = 'none'

        for seed in seedlist:
            StackLikelihood.stack_rosters(
                roster_dict,
                args.ttype,
                channels,
                jprior,
                sim_name,
                seed,
                args.clobber)


class ConvertCastro_SG(ScatterGather):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-convert-castro-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = ConvertCastro

    job_time = 600

    default_options = dict(ttype=defaults.common['ttype'],
                           specfile=defaults.common['specfile'],
                           targetlist=defaults.common['targetlist'],
                           jpriors=defaults.common['jpriors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           clobber=defaults.common['clobber'])

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args)
        if targets_yaml is None:
            return job_configs

        specfile = NAME_FACTORY.resolve_specfile(args)

        targets = load_yaml(targets_yaml)

        jpriors = args['jpriors']
        dry_run = args['dry_run']
        clobber = args['clobber']

        if is_not_null(sim):
            is_sim = True
            nsims = args['nsims']
            seed = args['seed']
        else:
            is_sim = False
            nsims = -1
            seed = -1

        for target_name, profile_list in targets.items():
            for profile in profile_list:
                for jprior in jpriors:
                    full_key = "%s:%s:%s" % (target_name, profile, jprior)
                    target_version = profile.split('_')[0]
                    name_keys = dict(target_type=ttype,
                                     target_name=target_name,
                                     target_version=target_version,
                                     profile=profile,
                                     jprior=jprior,
                                     fullpath=True)
                    if is_sim:
                        name_keys['sim_name'] = sim
                        sed_file = NAME_FACTORY.sim_sedfile(**name_keys)
                        j_value_yaml = NAME_FACTORY.sim_j_valuefile(
                            **name_keys)
                        outfile = NAME_FACTORY.sim_dmlikefile(**name_keys)
                        limitfile = NAME_FACTORY.sim_dmlimitsfile(**name_keys)
                        full_key += ":%s" % sim
                    else:
                        sed_file = NAME_FACTORY.sedfile(**name_keys)
                        j_value_yaml = NAME_FACTORY.j_valuefile(**name_keys)
                        outfile = NAME_FACTORY.dmlikefile(**name_keys)
                        limitfile = NAME_FACTORY.dmlimitsfile(**name_keys)

                    logfile = make_nfs_path(outfile.replace('.fits', '.log'))
                    job_config = dict(specfile=specfile,
                                      sed_file=sed_file,
                                      j_value_file=j_value_yaml,
                                      jprior=jprior,
                                      outfile=outfile,
                                      limitfile=limitfile,
                                      logfile=logfile,
                                      nsims=nsims,
                                      seed=seed,
                                      dry_run=dry_run,
                                      clobber=clobber)

                    job_configs[full_key] = job_config

        return job_configs


class StackLikelihood_SG(ScatterGather):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-stack-likelihood-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = StackLikelihood

    job_time = 120

    default_options = dict(ttype=defaults.common['ttype'],
                           specconfig=defaults.common['specfile'],
                           rosterlist=defaults.common['rosterlist'],
                           jpriors=defaults.common['jpriors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           clobber=defaults.common['clobber'])

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        jpriors = args['jpriors']
        dry_run = args['dry_run']
        clobber = args['clobber']
        sim = args['sim']

        if is_not_null(sim):
            is_sim = True
            nsims = args['nsims']
            seed = args['seed']
        else:
            is_sim = False
            nsims = -1
            seed = -1

        for jprior in jpriors:
            full_key = "%s:%s" % (jprior, sim)
            name_keys = dict(target_type=args['ttype'],
                             target_name='stacked',
                             jprior=jprior,
                             sim_name=sim,
                             fullpath=True)
            if is_sim:
                target_dir = NAME_FACTORY.sim_targetdir(**name_keys)
            else:
                target_dir = NAME_FACTORY.targetdir(**name_keys)

            logfile = os.path.join(
                target_dir, 'stack_%s_%s.log' %
                (jprior, sim))
            job_config = dict(ttype=args['ttype'],
                              specconfig=args['specconfig'],
                              rosterlist=args['rosterlist'],
                              jprior=jprior,
                              sim=sim,
                              logfile=logfile,
                              nsims=nsims,
                              seed=seed,
                              dry_run=dry_run,
                              clobber=clobber)

            job_configs[full_key] = job_config

        return job_configs


def register_classes():
    """Register these classes with the `LinkFactory` """
    ConvertCastro.register_class()
    ConvertCastro_SG.register_class()
    StackLikelihood.register_class()
    StackLikelihood_SG.register_class()
    SpecTable.register_class()
