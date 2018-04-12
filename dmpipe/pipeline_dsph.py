#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Scripts to run the all-sky diffuse analysis
"""
from __future__ import absolute_import, division, print_function

import sys
import os
import argparse


from fermipy.jobs.job_archive import JobArchive
from fermipy.jobs.lsf_impl import check_log
from fermipy.jobs.chain import Link, Chain

from dmpipe.dm_spectral import create_link_spec_table_builder,\
    create_link_stack_likelihood, create_sg_castro_convertor

from dmpipe.target_analysis import create_link_prepare_targets, create_sg_roi_analysis,\
    create_sg_sed_analysis


class Pipeline_dsph(Chain):
    """Small class to chain together the steps of the dSphs pipeline
    """
    default_options = dict(topdir=('dsph_flight', 'Top-level analysis directory.', str),
                           baseconfig=('config_baseline.yaml',
                                       'Template analysis configuration.', str),
                           spec_table=('dm_spectra.fits', 'FITS file with DM spectra', str),
                           targetlist=('target_list.yaml', 'Yaml file with list of targets', str),
                           rosterlist=('roster_list.yaml', 'Yaml file with list of rosters', str),
                           jprior=(None, 'Type of Prior on J-factor', str),
                           dry_run=(False, 'Dry run only', bool))

    def __init__(self, linkname, **kwargs):
        """C'tor
        """
        job_archive = kwargs.get('job_archive', None)
        sg_roi_analysis = create_sg_roi_analysis(linkname="%s.roi-analysis" % linkname,
                                                 mapping={'action': 'action_roi'},
                                                 job_archive=job_archive)
        sg_sed_analysis = create_sg_sed_analysis(linkname="%s.sed-analysis" % linkname,
                                                 mapping={'action': 'action_sed'},
                                                 job_archive=job_archive)
        sg_castro_conv = create_sg_castro_convertor(linkname="%s.dm-castro" % linkname,
                                                    mapping={'action': 'action_castro'},
                                                    job_archive=job_archive)
        link_stack_likelihood = create_link_stack_likelihood(linkname="%s.stack-likelihood" % linkname,
                                                             job_archive=job_archive)

        parser = argparse.ArgumentParser(usage='dmpipe-dsph-chain',
                                         description="Run dSphs analysis chain")
        Chain.__init__(self, linkname,
                       links=[sg_roi_analysis, sg_sed_analysis, sg_castro_conv,
                              link_stack_likelihood],
                       appname='dmpipe-dsph-chain',
                       options=Pipeline_dsph.default_options.copy(),
                       argmapper=self._map_arguments,
                       parser=parser,
                       **kwargs)

    def _map_arguments(self, input_dict):
        """Map from the top-level arguments to the arguments provided to
        the indiviudal links """
        output_dict = input_dict.copy()
        if input_dict.get('dry_run', False):
            action = 'check_status'
        else:
            action = 'run'
        output_dict['action_roi'] = action
        output_dict['action_sed'] = action
        output_dict['action_castro'] = action
        return output_dict

    def run_argparser(self, argv):
        """Initialize a link with a set of arguments using argparser
        """
        args = Link.run_argparser(self, argv)
        for link in self._links.values():
            link.run_link(stream=sys.stdout, dry_run=True)
        return args


def create_chain_dsph_pipeline(**kwargs):
    """Build and return a `Pipeline_dsph` object """
    ret_chain = Pipeline_dsph(linkname=kwargs.pop('linkname', 'Pipeline_dsph'))
    return ret_chain


def main_chain():
    """Energy point for running the entire Cosmic-ray analysis """

    job_archive = JobArchive.build_archive(job_archive_table='job_archive_dSphs.fits',
                                           file_archive_table='file_archive_dSphs.fits',
                                           base_path=os.path.abspath('.') + '/')

    the_chain = Pipeline_dsph('dsphs', job_archive=job_archive)
    args = the_chain.run_argparser(sys.argv[1:])
    logfile = "log_%s_top.log" % the_chain.linkname
    the_chain.archive_self(logfile)
    if args.dry_run:
        outstr = sys.stdout
    else:
        outstr = open(logfile, 'append')
    the_chain.run_chain(outstr, args.dry_run)
    if not args.dry_run:
        outstr.close()
    the_chain.finalize(args.dry_run)
    job_archive.update_job_status(check_log)


if __name__ == '__main__':
    main_chain()
