#!/usr/bin/env python
#

"""
Utilities to plot dark matter analyses
"""

import numpy as np

from fermipy import castro
from fermipy import stats_utils
from fermipy import sed_plotting
from fermipy.spectrum import DMFitFunction

from dmpipe.dm_spectral_utils import DMSpecTable


ENERGY_AXIS_LABEL = r'Energy [MeV]'
ENERGY_FLUX_AXIS_LABEL = r'Energy Flux [MeV s$^{-1}$ cm$^{-2}$]'
FLUX_AXIS_LABEL = r'Flux [ph s$^{-1}$ cm$^{-2}$]'
MASS_AXIS_LABEL = r"$m_{\chi}$ [GeV]"
SIGMAV_AXIS_LABEL = r'$\langle \sigma v \rangle$ [$cm^{3} s^{-1}$]'
TAU_AXIS_LABEL = r'$\tau_{\chi}$ [$s$]'
J_UNC_AXIS_LABEL = r'$\delta J / J$'
DELTA_LOGLIKE_AXIS_LABEL = r'$\Delta \log\mathcal{L}$'


def plot_dm_spectra_by_channel(dm_spec_table, mass=100.,
                               spec_type='eflux', ylims=(1e-12, 1e-8)):
    """ Make a plot of the DM spectra.

    Parameters
    ----------

    dm_spec_table : `DMSpecTable`
        The object with the spectra

    mass : float
        The DM particle mass in GeV

    spec_type : str
        Spectral type, one of 'flux' or 'eflux'

    ylims : tuple
        Y-axis limits for plot


    Returns
    -------

    fig : `matplotlib.Figure`
        The figure

    axis : `matplotlib.Axes`
        The plot axes

    leg : `matplotlib.Legend`
        The figure legend

    """
    import matplotlib.pyplot as plt

    chan_names = dm_spec_table.channel_names
    chan_ids = dm_spec_table.channel_map.keys()
    chan_idx_list = dm_spec_table.channel_map.values()
    energies = dm_spec_table.ebin_refs()

    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlim((energies[0], energies[-1]))
    axis.set_ylim(ylims)
    axis.set_xlabel(ENERGY_AXIS_LABEL)
    if spec_type == 'eflux':
        axis.set_ylabel(ENERGY_FLUX_AXIS_LABEL)
    elif spec_type == 'flux':
        axis.set_ylabel(FLUX_AXIS_LABEL)

    for chan, chan_id, idx_list in zip(chan_names, chan_ids, chan_idx_list):
        chan_masses = dm_spec_table.masses(chan_id).data
        mass_idx = np.abs(chan_masses - mass).argmin()
        table_idx = idx_list[mass_idx]
        spectrum = dm_spec_table._s_table[table_idx]["ref_%s" % spec_type]
        axis.plot(energies, spectrum, label=chan)

    leg = axis.legend(loc="best", ncol=2, fontsize=10)
    return fig, axis, leg


def plot_dm_spectra_by_mass(dm_spec_table, chan='bb',
                            spec_type='eflux', ylims=(1e-12, 1e-6)):
    """ Make a plot of the DM spectra.

    Parameters
    ----------

    dm_spec_table : `DMSpecTable`
        The object with the spectra

    chan : str
        The DM interaction channel

    spec_type : str
        Spectral type, one of 'flux' or 'eflux'

    ylims : tuple
        Y-axis limits for plot


    Returns
    -------

    fig : `matplotlib.Figure`
        The figure

    axis : `matplotlib.Axes`
        The plot axes

    leg : `matplotlib.Legend`
        The figure legend

    """
    import matplotlib.pyplot as plt

    chan_id = DMFitFunction.channel_rev_map[chan]
    chan_idx_list = dm_spec_table.channel_map[chan_id]
    energies = dm_spec_table.ebin_refs()

    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlim((energies[0], energies[-1]))
    axis.set_ylim(ylims)
    axis.set_xlabel(ENERGY_AXIS_LABEL)
    if spec_type == 'eflux':
        axis.set_ylabel(ENERGY_FLUX_AXIS_LABEL)
    else:
        axis.set_ylabel(FLUX_AXIS_LABEL)

    masses = dm_spec_table.masses(chan_id)
    for table_idx, mass in zip(chan_idx_list, masses):
        spectrum = dm_spec_table._s_table[table_idx]["ref_%s" % spec_type]
        axis.plot(energies, spectrum, label="%.1F GeV" % mass)

    leg = axis.legend(loc="best", ncol=2, fontsize=10)
    return fig, axis, leg


def plot_dm_castro(castro_dm, ylims=None, nstep=100, zlims=None, global_min=False):
    """ Make a color plot (1castro plot) of the delta log-likelihood as a function of
    DM particle mass and cross section.

    Parameters
    ----------

    castro_dm :  `DMCastroData`
        Object with the log-likelihood v. normalization for each mass 

    ylims      : tuple
        Y-axis limits for the plot

    nstep      : int
        Number of y-axis steps to plot for each energy bin

    zlims      : tuple
        z-axis limits

    global_min : bool
        If True plot likelihood w.r.t. the global minimimum.

    """
    if castro_dm.decay:
        ylabel = TAU_AXIS_LABEL
        if ylims is None:
            ylims = (1e+22, 1e+28)
    else:
        ylabel = SIGMAV_AXIS_LABEL
        if ylims is None:
            ylims = (1e-28, 1e-22)

    return sed_plotting.plotCastro_base(castro_dm,
                                        ylims=ylims,
                                        xlabel=MASS_AXIS_LABEL,
                                        ylabel=ylabel,
                                        nstep=nstep,
                                        zlims=zlims,
                                        global_min=global_min)


def plot_castro_nuiscance(xlims, ylims, zvals, zlims=None, decay=False):
    """ Make a castro plot including the effect of the nuisance parameter
    """
    import matplotlib.pyplot as plt
    from matplotlib import cm

    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.set_yscale('log')
    axis.set_xlim(xlims)
    axis.set_ylim(ylims)
    axis.set_xlabel(J_UNC_AXIS_LABEL)
    if decay:
        axis.set_ylabel(TAU_AXIS_LABEL)
    else:
        axis.set_ylabel(SIGMAV_AXIS_LABEL)

    if zlims is None:
        zmin = 0
        zmax = 10.
    else:
        zmin = zlims[0]
        zmax = zlims[1]

    image = axis.imshow(zvals, extent=[xlims[0], xlims[-1], ylims[0], ylims[-1]],
                        origin='lower', aspect='auto', interpolation='nearest',
                        vmin=zmin, vmax=zmax, cmap=cm.jet_r)
    return fig, axis, image


def plot_nll(nll_dict, xlims=None, nstep=50, ylims=None, decay=False):
    """ Plot the -log(L) as a function of sigmav for each object in a dict
    """
    import matplotlib.pyplot as plt

    if xlims is None:
        xmin = 1e-28
        xmax = 1e-24
    else:
        xmin = xlims[0]
        xmax = xlims[1]

    xvals = np.logspace(np.log10(xmin), np.log10(xmax), nstep)
    fig = plt.figure()
    axis = fig.add_subplot(111)

    axis.set_xlim((xmin, xmax))
    if ylims is not None:
        axis.set_ylim((ylims[0], ylims[1]))

    if decay:
        axis.set_xlabel(TAU_AXIS_LABEL)
    else:
        axis.set_xlabel(SIGMAV_AXIS_LABEL)

    axis.set_ylabel(DELTA_LOGLIKE_AXIS_LABEL)
    axis.set_xscale('log')

    for lab, nll in nll_dict.items():
        yvals = nll.interp(xvals)
        yvals -= yvals.min()
        axis.plot(xvals, yvals, label=lab)

    leg = axis.legend(loc="upper left")
    return fig, axis, leg


def plot_comparison(nll, nstep=25, xlims=None, decay=False):
    """ Plot the comparison between differnt version of the -log(L)
    """
    import matplotlib.pyplot as plt

    if xlims is None:
        xmin = nll._lnlfn.interp.xmin
        xmax = nll._lnlfn.interp.xmax
    else:
        xmin = xlims[0]
        xmax = xlims[1]

    xvals = np.linspace(xmin, xmax, nstep)
    yvals_0 = nll.straight_loglike(xvals)
    yvals_1 = nll.profile_loglike(xvals)
    yvals_2 = nll.marginal_loglike(xvals)

    ymin = min(yvals_0.min(), yvals_1.min(), yvals_2.min(), 0.)
    ymax = max(yvals_0.max(), yvals_1.max(), yvals_2.max(), 0.5)

    fig = plt.figure()
    axis = fig.add_subplot(111)

    axis.set_xlim((xmin, xmax))
    axis.set_ylim((ymin, ymax))

    if decay:
        axis.set_xlabel(SIGMAV_AXIS_LABEL)
    else:
        axis.set_xlabel(SIGMAV_AXIS_LABEL)
    axis.set_ylabel(DELTA_LOGLIKE_AXIS_LABEL)

    axis.plot(xvals, yvals_0, 'r', label=r'Simple $\log\mathcal{L}$')
    axis.plot(xvals, yvals_1, 'g', label=r'Profile $\log\mathcal{L}$')
    #axis.plot(xvals,yvals_2,'b', label=r'Marginal $\log\mathcal{L}$')

    leg = axis.legend(loc="upper left")

    return fig, axis, leg


def plot_stacked(sdict, xlims, ibin=0, decay=False):
    """ Stack a set of -log(L) curves and plot the stacked curve
    as well as the individual curves
    """
    import matplotlib.pyplot as plt
    ndict = {}

    for key, val in sdict.items():
        ndict[key] = val[ibin]

    #mles = np.array([n.mle() for n in ndict.values()])

    fig = plt.figure()
    axis = fig.add_subplot(111)

    xvals = np.linspace(xlims[0], xlims[-1], 100)

    axis.set_xlim((xvals[0], xvals[-1]))
    if decay:
        axis.set_xlabel(SIGMAV_AXIS_LABEL)
    else:
        axis.set_xlabel(SIGMAV_AXIS_LABEL)
    axis.set_ylabel(DELTA_LOGLIKE_AXIS_LABEL)

    for key, val in ndict.items():
        yvals = val.interp(xvals)
        if key.lower() == "stacked":
            axis.plot(xvals, yvals, lw=3, label=key)
        else:
            yvals -= yvals.min()
            axis.plot(xvals, yvals, label=key)
    leg = axis.legend(loc="upper left", fontsize=10, ncol=2)
    return fig, axis, leg


def plot_limits_from_arrays(ldict, xlims, ylims, bands=None, decay=False):
    """ Plot the upper limits as a function of DM particle mass and cross section.

    Parameters
    ----------

    ldict      : dict
        A dictionary of strings pointing to pairs of `np.array` objects,
        The keys will be used as labels

    xlims      : tuple
        x-axis limits for the plot

    ylims      : tupel
        y-axis limits for the plot

    bands      : dict
        Dictionary with the expected limit bands

    decay : bool
        Plot limits for decay instead of annihilation
    

    Returns
    -------

    fig : `matplotlib.Figure`
        The figure

    axis : `matplotlib.Axes`
        The plot axes

    leg : `matplotlib.Legend`
        The figure legend

    """
    import matplotlib.pyplot as plt

    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.set_xlabel(MASS_AXIS_LABEL)
    if decay:
        axis.set_ylabel(TAU_AXIS_LABEL)
    else:
        axis.set_ylabel(SIGMAV_AXIS_LABEL)

    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlim((xlims[0], xlims[1]))
    axis.set_ylim((ylims[0], ylims[1]))

    if bands is not None:
        plot_expected_limit_bands(axis, bands)

    for key, val in ldict.items():
        xvals = val[0]
        yvals = val[1]
        if key.lower() == "stacked":
            axis.plot(xvals, yvals, lw=3, label=key)
        else:
            axis.plot(xvals, yvals, label=key)

    leg = axis.legend(loc="upper left")  # ,fontsize=10,ncol=2)
    return fig, axis, leg


def plot_expected_limit_bands(axis, bands):
    """Add the expected limit bands to a plot

    Parameters
    ----------

    axis : `matplotlib.Axes`
        The axes we are adding the bands to

    bands : dict
        Dictionary with the bands

    """
    masses = bands['masses']

    axis.fill_between(masses, bands['q02'], bands['q97'], color='yellow')
    axis.fill_between(masses, bands['q16'], bands['q84'], color='green')
    axis.plot(masses, bands['median'], color='gray')


def plot_mc_truth(axis, mc_model, decay=False):
    """Add a marker showing the true injected signal to the plot

    Parameters
    ----------

    axis : `matplotlib.Axes`
        The axes we are adding the bands to

    mc_model : dict
        Dictionary with truth

    decay : bool
        Plot value for decay instead of annihilation
    
    """
    if decay:
        norm = mc_model['tau']['value']
    else:
        norm = mc_model['sigmav']['value']
    mass = mc_model['mass']['value']
    axis.scatter([mass], [norm])


def plot_limits(sdict, xlims, ylims, alpha=0.05, decay=False):
    """ Plot the upper limits as a function of DM particle mass and cross section.

    Parameters
    ----------

    sdict      : dict
        A dictionary of `DMCastroData` objects with the log-likelihood v. normalization for each energy bin

    xlims      : tuple
        x-axis limits

    ylims      : tuple
        y-axis limits

    alpha      : float
        Confidence level to use in setting limits = 1 - alpha

    decay : bool
        Plot limits for decay instead of annihilation

    """

    ldict = {}
    for key, val in sdict.items():
        ldict[key] = (val.masses, val.getLimits(alpha))
    return plot_limits_from_arrays(ldict, xlims, ylims, decay=decay)


def compare_limits(sdict, xlims, ylims, alpha=0.05, decay=False):
    """ Plot the upper limits as a functino of DM particle mass and cross section.

    Paramters
    ---------

    sdict      : dict
        Dictionary with limits and keys

    xlims      : tuple
        x-axis limits

    ylims      : tuple
        y-axis limits

    alpha      : float 
        Confidence level to use in setting limits = 1 - alpha
  
    decay : bool
        Plot limits for decay instead of annihilation


    Returns
    -------

    fig : `matplotlib.Figure`
        The figure

    axis : `matplotlib.Axes`
        The plot axes

    leg : `matplotlib.Legend`
        The figure legend

    """
    import matplotlib.pyplot as plt

    fig = plt.figure()
    axis = fig.add_subplot(111)

    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlim((xlims[0], xlims[1]))
    axis.set_ylim((ylims[0], ylims[1]))
    axis.set_xlabel(MASS_AXIS_LABEL)
    if decay:
        axis.set_ylabel(TAU_AXIS_LABEL)
    else:
        axis.set_ylabel(SIGMAV_AXIS_LABEL)

    for key, val in sdict.items():
        xvals = val.masses
        yvals = val.getLimits(alpha)
        axis.plot(xvals, yvals, label=key)

    leg = axis.legend(loc="upper left", fontsize=10, ncol=2)
    return fig, axis, leg


def plot_limit(dm_castro_data, ylims, alpha=0.05):
    """Plot the limit curve for a given DMCastroData object
 
    Parameters
    ----------

    dm_castro_data :  `DMCastroData`
        Object with the log-likelihood v. normalization for each mass 

    ylims      : tuple
        Y-axis limits for the plot

    alpha      : float
        Confidence level to use in setting limits = 1 - alpha


    Returns
    -------

    fig : `matplotlib.Figure`
        The figure

    axis : `matplotlib.Axes`
        The plot axes

    """
    import matplotlib.pyplot as plt

    xbins = dm_castro_data.masses
    xmin = xbins[0]
    xmax = xbins[-1]

    fig = plt.figure()
    axis = fig.add_subplot(111)

    axis.set_xscale('log')
    axis.set_yscale('log')
    axis.set_xlabel(MASS_AXIS_LABEL)
    if dm_castro_data.decay:
        axis.set_ylabel(TAU_AXIS_LABEL)
    else:
        axis.set_ylabel(SIGMAV_AXIS_LABEL)
    axis.set_xlim((xmin, xmax))

    if ylims is not None:
        axis.set_ylim((ylims[0], ylims[1]))

    if decay:
        yvals = dm_castro_data.getLimits(1.0 - alpha)
    else:
        yvals = dm_castro_data.getLimits(alpha)
    if yvals.shape[0] == xbins.shape[0]:
        xvals = xbins
    else:
        xvals = np.sqrt(xbins[0:-1] * xbins[1:])
    axis.plot(xvals, yvals)

    return fig, axis

