#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Scripts to run the all-sky diffuse analysis
"""
from __future__ import absolute_import, division, print_function

import sys
import os
import argparse

import yaml

from fermipy.jobs.job_archive import JobArchive
from fermipy.jobs.lsf_impl import check_log

from dmpipe.dm_spectral import create_link_castro_convertor, create_link_spec_table_builder,\
    create_link_stack_likelihood, create_sg_castro_convertor

from dmpipe.target_analysis import create_link_prepare_targets, create_sg_roi_analysis, create_sg_sed_analysis


BUILDER_DICT = {'create_link_castro_convertor':create_link_castro_convertor,
                'create_link_spec_table_builder':create_link_spec_table_builder,
                'create_link_stack_likelihood':create_link_stack_likelihood,
                'create_sg_castro_convertor':create_sg_castro_convertor,
                'create_link_prepare_targets':create_link_prepare_targets,
                'create_sg_roi_analysis':create_sg_roi_analysis,
                'create_sg_sed_analysis':create_sg_sed_analysis}


def build_analysis_link(linktype, **kwargs):
    """Build and return a `fermipy.jobs.Link` object to run a
    part of the analysis"""

    builder_name = 'create_%s'%linktype
    try:
        builder_func = BUILDER_DICT[builder_name]
    except KeyError:
        raise KeyError("Could not build an analysis link using a creator function %s"%builder_name)
    return builder_func(**kwargs)


if __name__ == '__main__':

    JOB_ARCHIVE = JobArchive.build_archive(job_archive_table='job_archive_temp2.fits',
                                           file_archive_table='file_archive_temp2.fits',
                                           base_path=os.path.abspath('.')+'/')

    PARSER = argparse.ArgumentParser(usage="diffuse_analysis.py [options] analyses",
                                     description="Run a high level analysis")
    PARSER.add_argument('--config', type=str, default=None, help="Yaml configuration file")
    PARSER.add_argument('--dry_run', action='store_true', help="Dry run only")
    PARSER.add_argument('analyses', nargs='+', type=str, help="Analysis steps to run")

    ARGS = PARSER.parse_args()

    CONFIG = yaml.load(open(ARGS.config))

    ANALYSES = ARGS.analyses
    if 'ALL' in ANALYSES:
        ANALYSES = ['spec_table', 'prepare_targets',
                    'roi_analysis', 'sed_analysis', 'castro_convertor',
                    'stack_likelihood']  

    for ANALYSIS in ANALYSES:
        ANALYSIS_CONFIG = CONFIG[ANALYSIS]
        LINKTYPE = ANALYSIS_CONFIG.pop('linktype')
        #LINK_KWARGS = dict(linkname=LINKTYPE)
        #LINK = build_analysis_link(LINKTYPE, **LINK_KWARGS)
        LINK = build_analysis_link(LINKTYPE)
        LINK.update_args(ANALYSIS_CONFIG)
        logfile = logfile="log_%s_top.log"%ANALYSIS
        LINK.register_self(logfile=logfile)
        JOB_ARCHIVE.register_jobs(LINK.get_jobs())
        if ARGS.dry_run:
            outstr = sys.stdout
        else:
            outstr = open(logfile, 'append')
        print ("Logfile = %s"%logfile)
        LINK.run(outstr, ARGS.dry_run)
        if not ARGS.dry_run:
            outstr.close()

    JOB_ARCHIVE.file_archive.update_file_status()
    JOB_ARCHIVE.write_table_file()
    JOB_ARCHIVE.update_job_status(check_log)
