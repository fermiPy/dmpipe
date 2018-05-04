#!/usr/bin/env python
#

"""
Interface to search targets and target lists
"""

import argparse

from astropy import table

from dmsky.roster import RosterLibrary
from dmsky.file_io import table as dm_table


class DMTargetFactory(object):
    """ A class to build a set of DM target objects from a dmsky Roster object
    """
    rlib = RosterLibrary()
    colNames = ['name', 'ra', 'dec', 'glat', 'glon', 'j_integ',
                'd_integ', 'j_sigma', 'd_sigma', 'distance', 'proftype']

    def __init__(self):
        """ C'tor, trivial
        """

    @staticmethod
    def create_roster(filepaths):
        """ Create a roster from input yaml files
        """
        roster = None
        for fname in filepaths:
            if roster is None:
                roster = DMTargetFactory.rlib.create_roster(fname)
            else:
                roster += DMTargetFactory.rlib.create_roster(fname)
        return roster

    @staticmethod
    def make_table(filepaths):
        """ Builld a dmsky dm_table object from input yaml files
        """
        roster = DMTargetFactory.create_roster(filepaths)
        tab = dm_table.make_table_for_roster(DMTargetFactory.colNames, roster)
        return tab, roster

    @staticmethod
    def read_table(filepath):
        """ Read a FITS table
        """
        tab = table.Table.read(filepath)
        return tab

    @staticmethod
    def get_sky_crds(tab):
        """ Extract coordinates from a FITS table
        """
        return (tab['glon'], tab['glat'])


if __name__ == "__main__":

    # Argument defintion
    USAGE = "usage: %(prog)s [input]"
    DESCRIPTION = "Make a FITS file with a target table"

    PARSER = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)
    PARSER.add_argument('--output', '-o', default=None, help='Output file.')
    PARSER.add_argument(
        '--clobber',
        action='store_true',
        help='Overwrite output file.')
    PARSER.add_argument('input', nargs="*", help='Input roster.')

    ARGS = PARSER.parse_args()

    TAB, RO = DMTargetFactory.make_table(ARGS.input)

    if ARGS.output:
        TAB.write(ARGS.output)
