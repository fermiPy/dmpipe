# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function

import os
import argparse

import numpy as np
from astropy.table import Table, Column

from fermipy.castro import CastroData
from fermipy.spectrum import DMFitFunction
from dmpipe.dm_spectral import DMSpecTable, DMCastroData

import dmsky.roster
import dmsky.skymap


def load_sed_from_txt(sedfile):
    """Load an SED from a txt file in the format used for online materials
    associated with dsph papers.."""

    # FIXME: Discover energy binning
    sed = np.loadtxt(sedfile, unpack=True)

    emin = np.unique(sed[0])
    emax = np.unique(sed[1])
    nbin = len(emin)
    sed = sed.reshape((4, nbin, -1))

    ref_spec = ReferenceSpec(emin, emax,
                             np.ones(nbin), np.ones(nbin),
                             np.ones(nbin), np.ones(nbin))
    cd = CastroData(sed[2], -1.0 * sed[3], ref_spec, 'eflux')
    return cd


def compute_limits(sedfile, roster, chan, masses, alpha=0.05, apply_prior=False, outprefix='lims2'):
    """Generate DM limits from an SED file.

    Parameters
    ----------
    sedfile : str
        Path to SED file.  

    """

    dm_tab = None
    names = []
    cd_stack = []
    lims_out = []
    ts_out = []

    # Loop over targets in roster
    for name, target in roster.items():

        names += [name]
        if len(roster) > 1:
            cd = CastroData.create_from_sedfile(sedfile, target=name)
        else:
            cd = CastroData.create_from_sedfile(sedfile)

        jfactor = target.params['j_integ'].value
        jsigma = target.params['j_sigma'].value
        if dm_tab is None:
            dm_tab = DMSpecTable.create(
                cd.refSpec.emin, cd.refSpec.emax, [], masses)

        chan_code = DMFitFunction.channel_rev_map[chan]
        if apply_prior:
            jfactor = {'j_value': jfactor,
                       'functype': 'lgauss_like', 'mu': 1.0, 'sigma': jsigma}
        cd_dm = dm_tab.convert_castro_data(
            cd, chan_code, 'eflux', jfactor=jfactor)

        print(name, jfactor)
        lims_out += [cd_dm.getLimits(alpha)]
        ts_out += [np.clip(cd_dm.ts_vals(),0,None)]
        cd_stack += [cd_dm]

    cd_dm_comb = DMCastroData.create_from_stack(cd_stack)
    lims_out += [cd_dm_comb.getLimits(alpha)]
    ts_out += [cd_dm_comb.ts_vals()]
    names += ['combined']
    masses = np.vstack([masses for i in range(len(names))])
    lims_out = np.vstack(lims_out)
    ts_out = np.vstack(ts_out)

    header = 'Channel: %s\n' % chan
    header += 'CL: %.3f\n' % (1.0 - alpha)
    header += 'Column 00: Mass (GeV)\n'
    for i, name in enumerate(names):
        header += 'Column %02i: %s sigma-v UL (cm^3 s^-1)\n' % (i + 1, name)
    header += 'Column %02i: combined sigma-v UL (cm^3 s^-1)\n' % (
        len(names) + 1)

    cols = [Column(name='name', data=names),
            Column(name='mass', data=masses, unit='GeV'),
            Column(name='sigmav_ul', data=lims_out, unit='cm^3 / s'),
            Column(name='ts', data=ts_out),
            ]

    text_out = np.vstack([masses[0]] + [lims_out])

    tab = Table(cols)
    tab.write('%s_%s.fits' % (outprefix, chan), overwrite=True)
    np.savetxt('%s_%s.txt' % (outprefix, chan),
               text_out.T, fmt='%12.5g', header=header)


def main(args=None):

    usage = "usage: %(prog)s [options]"
    description = ('Generate DM limits from a SED file.')
    parser = argparse.ArgumentParser(usage=usage, description=description)

    parser.add_argument('--roster', default=None,
                        help='Set the target roster from dmsky.')

    parser.add_argument('--target', default=None,
                        help='Set the name of the target.')

    parser.add_argument('--channel', default=None,
                        help='Set the annihilation channel.')

    parser.add_argument('--jfactor', default=None,
                        help='Set the J-factor of the target.')

    parser.add_argument('--jsigma', default=None,
                        help='Set the J-factor uncertainty of the target.')

    parser.add_argument('--alpha', default=0.05,
                        help='Set the confidence level to be used for the limit calculation.')

    parser.add_argument('--outprefix', default='limits',
                        help='Set the prefix for output files.')

    parser.add_argument('sedfile',
                        help='Set the path to a file containing SEDs for one or more targets.  '
                        'For a multi-object SED file the target name should be specified by a name '
                        'column.  If the file contains only a single target then its name must be '
                        'specified with the target option.')

    args = parser.parse_args(args)

    library = dmsky.roster.RosterLibrary()
    #dsphs = library.create_roster('ackermann2015_dsphs')
    #roster = library.create_roster('albert2017_dsphs')
    channels = ['bb', 'ee', 'mumu', 'tautau', 'uu', 'ww']
    if args.channel is not None:
        channels = args.channel.split(',')
    
    masses_high = [100.00, 158.10, 250.00, 353.60, 500.00, 707.00,
                   1000.00, 1581.00, 2500.00, 3536.00, 5000.00, 7070.00, 10000.00]

    masses_table = dict(
        uu=[2.00, 3.16, 5.00, 7.07, 10.00, 15.81, 25.00, 35.36, 50.00, 70.70, ] + masses_high,
        dd=[2.00, 3.16, 5.00, 7.07, 10.00, 15.81, 25.00, 35.36, 50.00, 70.70, ] + masses_high,
        cc=[2.00, 3.16, 5.00, 7.07, 10.00, 15.81, 25.00, 35.36, 50.00, 70.70, ] + masses_high,
        ss=[2.00, 3.16, 5.00, 7.07, 10.00, 15.81, 25.00, 35.36, 50.00, 70.70, ] + masses_high,
        bb=[6.00, 7.746, 10.00, 15.81, 25.00, 35.36, 50.00, 70.70, ] + masses_high,
        ww=[81.00, ] + masses_high,
        ee=[2.00, 3.16, 5.00, 7.07, 10.00, 15.81, 25.00, 35.36, 50.00, 70.70, ] + masses_high,
        mumu=[2.00, 3.16, 5.00, 7.07, 10.00, 15.81, 25.00, 35.36, 50.00, 70.70, ] + masses_high,
        tautau=[2.00, 3.16, 5.00, 7.07, 10.00, 15.81, 25.00, 35.36, 50.00, 70.70, ] + masses_high,
    )

    print (args.jfactor, args.jsigma)
    if args.roster is None:
        target = dmsky.targets.Dwarf(j_integ=float(args.jfactor), j_sigma=float(args.jsigma),
                                     name=args.target,
                                     profile={'type': 'nfw'},
                                     abbr=args.target, title=args.target)
        roster = dmsky.roster.Roster([target])
    else:
        roster = library.create_roster(args.roster)
        if args.target is not None:
            # Remove all objects in roster but the target
            roster = dmsky.roster.Roster([roster[args.target]])

    for chan in channels:
        masses = masses_table[chan]
        compute_limits(args.sedfile, roster, chan, masses, apply_prior=True,
                       outprefix=args.outprefix, alpha=args.alpha)


if __name__ == "__main__":
    main()
