#!/usr/bin/env python
#

"""
Interface to Dark Matter spectra
"""

import sys
import yaml
import numpy as np

from astropy.table import Table, Column

from fermipy import castro
from fermipy import fits_utils
from fermipy import stats_utils
from fermipy import skymap

from dmfit.dm_fit_spec import DMFitFunction



class DMCastroData(castro.CastroData_Base):
    """ This class wraps the data needed to make a "Castro" plot,
    namely the log-likelihood as a function of <sigma v> * J (or D) for a
    series of DM masses
    """

    def __init__(self, norm_vals, nll_vals, norm_value, channel, masses, astro_value,
                 astro_prior=None, prior_applied=True):
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
        """
        self._masses = masses
        self._norm_value = norm_value
        if self._norm_value is None:
            self._norm_value = 1.0
        self._astro_value = astro_value
        self._astro_prior = astro_prior
        self._prior_applied = prior_applied
        if isinstance(channel, str):
            self._channel = DMFitFunction.channel2int(channel)
        else:
            self._channel = channel
        super(DMCastroData, self).__init__(norm_vals, nll_vals, norm_type="sigmav")

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
            return "None"
        return self.astro_prior.funcname

    @property
    def prior_applied(self):
        """ Has the prior already been applied """
        return self._prior_applied

    @property
    def channel(self):
        """ The string specifying the decay channel """
        return self._channel

    @staticmethod
    def create_from_stack(components, nystep=100, ylims=(1e-28, 1e-22), weights=None):
        """ Create a DMCastroData object by stacking a series of DMCastroData objects
        """
        if len(components) == 0:
            return None
        shape = (components[0].nx, nystep)
        norm_vals, nll_vals = castro.CastroData_Base.stack_nll(shape, components, ylims, weights)
        return DMCastroData(norm_vals, nll_vals, 1.0,
                            components[0].channel, components[0].masses, astro_value=None)

    @staticmethod
    def create_from_tables(tab_s, tab_m):
        """ Create a DMCastroData object from likelihood scan and mass tables
        """
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
                              sigma=prior_sigma)
            prior = stats_utils.create_prior_functor(prior_dict)
        else:
            prior = None
            prior_applied = True

        return DMCastroData(norm_vals, nll_vals, norm, channel,
                            masses, astro_value, prior, prior_applied)

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

        collist = [col_norm, col_normv, col_dll, col_astro_val]
        valdict = {"NORM": self._norm_value,
                   "NORM_SCAN": self._norm_vals / self._norm_value,
                   "DLOGLIKE_SCAN": -1 * self._nll_vals,
                   "ASTRO_VALUE": self.astro_value}

        if self._astro_prior is not None:
            col_prior_type = Column(name="PRIOR_TYPE", dtype="S16")
            col_prior_mean = Column(name="PRIOR_MEAN", dtype=float)
            col_prior_sigma = Column(name="PRIOR_SIGMA", dtype=float)
            col_prior_applied = Column(name="PRIOR_APPLIED", dtype=bool)
            collist += [col_prior_type, col_prior_mean, col_prior_sigma, col_prior_applied]
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
        col_masses = Column(name="MASSES", dtype=float, shape=self._masses.shape)
        col_channel = Column(name="CHANNEL", dtype=int)
        tab = Table(data=[col_masses, col_channel])
        tab.add_row({"MASSES": self._masses,
                     "CHANNEL": self._channel})
        return tab


class DMSpecTable(object):
    """ Version of the DM spectral tables in tabular form
    """
    @staticmethod
    def make_ebounds_table(emin, emax, eref):
        """ Construct the energy bounds table
        """
        from astropy import table

        col_emin = table.Column(name="E_MIN", dtype=float, unit="MeV", data=emin)
        col_emax = table.Column(name="E_MAX", dtype=float, unit="MeV", data=emax)
        col_eref = table.Column(name="E_REF", dtype=float, unit="MeV", data=eref)

        tab = table.Table(data=[col_emin, col_emax, col_eref])
        return tab

    @staticmethod
    def make_spectra_table(nebins, spec_dict):
        """ Construct the spectral values table
        """
        from astropy import table

        col_masses = table.Column(name="REF_MASS", dtype=float, unit="GeV",
                                  data=spec_dict['MASS'])
        col_chans = table.Column(name="REF_CHAN", dtype=int,
                                 data=spec_dict['CHAN'])
        col_dfde = table.Column(name="REF_DFDE", dtype=float, shape=nebins, unit="ph / (MeV cm2 s)",
                                data=spec_dict['DFDE'])
        col_flux = table.Column(name="REF_FLUX", dtype=float, shape=nebins, unit="ph / (cm2 s)",
                                data=spec_dict['FLUX'])
        col_eflux = table.Column(name="REF_EFLUX", dtype=float, shape=nebins, unit="MeV / (cm2 s)",
                                 data=spec_dict['EFLUX'])

        table = table.Table(data=[col_masses, col_chans, col_dfde, col_flux, col_eflux])
        return table

    def __init__(self, e_table, s_table, ref_vals):
        """ C'tor to build this object from energy binning and spectral values tables.
        """
        self._e_table = e_table
        self._s_table = s_table
        self._ref_vals = ref_vals

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

    def spectrum(self, chan, mass, spec_type):
        """ Return the spectrum for a particular channel an mass
        """
        mask = (self._s_table["REF_CHAN"] == chan) & (
            np.abs(self._s_table["REF_MASS"] - mass) < 1e-9)
        spec_vals = self._s_table[mask]["REF_%s" % spec_type].data
        return spec_vals

    def masses(self, chan):
        """ Return the array of masses for a given channel
        """
        mask = (self._s_table["REF_CHAN"] == chan)
        return self._s_table[mask]["REF_MASS"]

    def ebin_edges(self):
        """ Return an array with the energy bin edges
        """
        return np.hstack([self._e_table["E_MIN"].data, self._e_table["E_MAX"].data])

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
        dref = {"REF_SIGV": hin["REF_SIGV"],
                "REF_J": hin["REF_J"]}
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

    @staticmethod
    def create_from_config(configfile, channels, masses):
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
        evals = np.sqrt(ebin_edges[0:-1] * ebin_edges[1:])

        ichans = [DMFitFunction.channel2int(c) for c in channels]

        init_params = np.array([1e-26, 100.])
        dmf = DMFitFunction(init_params)

        nrow = len(ichans) * len(masses)
        irow = 0

        dfde = np.ndarray((nrow, nebins))
        flux = np.ndarray((nrow, nebins))
        eflux = np.ndarray((nrow, nebins))
        masses_out = np.ndarray((nrow))
        channels = np.ndarray((nrow), int)

        #for i, chan in zip(ichans, channels):
        for i in ichans:
            DMFitFunction.set_default_par_value('channel0', i)
            for mass in masses:
                init_params[1] = mass
                masses_out[irow] = mass
                channels[irow] = i
                dfde[irow].flat = dmf.dfde(evals)
                flux[irow].flat = dmf.flux(ebin_edges[0:-1], ebin_edges[1:], init_params)
                eflux[irow].flat = dmf.eflux(ebin_edges[0:-1], ebin_edges[1:], init_params)
                irow += 1

        spec_dict = {"DFDE": dfde,
                     "FLUX": flux,
                     "EFLUX": eflux,
                     "MASS": masses_out,
                     "CHAN": channels}

        ref_vals = {"REF_J": DMFitFunction.default_params['norm'],
                    "REF_SIGV": DMFitFunction.default_params['sigmav']}

        return DMSpecTable.create_from_data(ebin_edges[0:-1], ebin_edges[1:],
                                            evals, spec_dict, ref_vals)

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

    def convert_castro_data(self, castro_data, channel, norm_type, jfactor=None):
        """ Convert CastroData object, i.e., Likelihood as a function of
        flux and energy flux, to a DMCastroData object, i.e., Likelihood as
        a function of DM mass and sigmav
        """
        if not self.check_energy_bins(castro_data.refSpec):
            raise ValueError("CastroData energy binning does not match")
        mask = self._s_table['REF_CHAN'] == channel
        masses = self._s_table[mask]['REF_MASS'].data
        spec_vals = self._s_table[mask]['REF_%s' % norm_type].data
        nmass = len(masses)

        # Get the reference values
        ref_sigmav = self._ref_vals["REF_SIGV"]
        ref_norm = self._ref_vals["REF_J"]

        j_prior = None
        if jfactor is None:
            # Just use the reference values
            j_ref = ref_norm
            norm_factor = 1.
        elif isinstance(jfactor, float):
            # Rescale the normalization values
            j_ref = jfactor
            norm_factor = jfactor / ref_norm
        elif isinstance(jfactor, dict):
            j_ref = jfactor.pop('jvalue')
            norm_factor = j_ref / ref_norm
            j_prior = stats_utils.create_prior_functor(jfactor)

        norm_limits = castro_data.getLimits(1e-5)
        spec_vals *= norm_factor

        n_scan_pt = 100

        norm_vals = np.ndarray((nmass, n_scan_pt))
        dll_vals = np.ndarray((nmass, n_scan_pt))
        mle_vals = np.ndarray((nmass))

        #for i, mass in enumerate(masses):
        for i in range(nmass):
            max_ratio = 1. / ((spec_vals[i] / norm_limits).max())
            log_max_ratio = np.log10(max_ratio)
            norm_vals[i][0] = 0.
            norm_vals[i][1:] = np.logspace(log_max_ratio - 2, log_max_ratio + 2, n_scan_pt - 1)
            test_vals = (np.expand_dims(spec_vals[i], 1) * (np.expand_dims(norm_vals[i], 1).T))
            dll_vals[i, 0:] = castro_data(test_vals)
            mle_vals[i] = norm_vals[i][dll_vals[i].argmin()]
            mle_ll = dll_vals[i].min()
            dll_vals[i] -= mle_ll

            if j_prior is not None:
                lnlfn = castro_data.LnLFn(norm_vals[i], dll_vals[i], 'dummy')
                lnlfn_prior = stats_utils.LnLFn_norm_prior(lnlfn, j_prior)
                dll_vals[i, 0:] = lnlfn_prior(norm_vals[i])

        norm_vals *= (ref_sigmav)
        dm_castro = DMCastroData(norm_vals, dll_vals, norm_factor,
                                 channel, masses, j_ref, j_prior)
        return dm_castro

    def convert_tscube(self, tscube, channel, norm_type):
        """ Convert TSCube object, i.e., Likelihood as a function of
        flux and energy flux for a set of pixels to a set of DMCastroData objects,
        i.e., Likelihood as a function of DM mass and sigmav for those pixels.
        """
        if not self.check_energy_bins(tscube.refSpec):
            raise ValueError("TSCube energy binning does not match")

        wcs = tscube.tscube.wcs
        mask = self._s_table['REF_CHAN'] == channel
        masses = self._s_table[mask]['REF_MASS'].data
        spec_vals = self._s_table[mask]['REF_%s' % norm_type].data
        nmasses = len(masses)

        # Get the reference values
        ref_sigmav = self._ref_vals["REF_SIGV"]
        ref_norm = self._ref_vals["REF_J"]
        j_ref = ref_norm
        norm_factor = 1.

        nvals = tscube.nvals
        cube_shape = (len(masses), tscube.tscube.counts.shape[1], tscube.tscube.counts.shape[2])

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

            #for i, m in enumerate(masses):
            for i in range(nmasses):
                max_ratio = 1. / ((spec_vals[i] / norm_limits).max())
                log_max_ratio = np.log10(max_ratio)
                norm_vals[i][0] = 0.
                norm_vals[i][1:] = np.logspace(log_max_ratio - 2, log_max_ratio + 2, n_scan_pt - 1)
                test_vals = (np.expand_dims(spec_vals[i], 1) * (np.expand_dims(norm_vals[i], 1).T))
                dll_vals[i, 0:] = castro_data(test_vals)
                mle_vals[i] = norm_vals[i][dll_vals[i].argmin()]
                mle_ll = dll_vals[i].min()
                dll_vals[i] -= mle_ll

            norm_vals *= (ref_sigmav)
            dm_castro = DMCastroData(norm_vals, dll_vals, norm_factor, channel, masses, j_ref)
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


if __name__ == "__main__":

    pass
