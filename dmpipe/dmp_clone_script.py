#!/usr/bin/env python
#

# Description
"""
A simple function to clone a script for each directory in an ROI set
"""

import os
import argparse


from fermipy import utils
from dmpipe import dmp_roi


def clone_script(basedir, roi_set, script_in, options):
    """ Clone a bash script for every ROI in a DMRoiSet
    """
    bash_script = """
cat $0
python {script} {options}
"""

    print(bash_script.format(script=script_in, options=options))

    scriptdir = os.path.join(basedir, 'scripts')
    utils.mkdir(scriptdir)
    os.system('cp %s %s' % (script_in, scriptdir))

    for name in roi_set.roi_dict.keys():

        dirname = os.path.join(basedir, name)
        utils.mkdir(dirname)

        print (dirname)

        script = os.path.basename(script_in)
        scriptpath = os.path.abspath(os.path.join(dirname, script))
        os.system('ln -sf %s %s' % (os.path.abspath(os.path.join(scriptdir, script)),
                                    scriptpath))
        runscript = os.path.abspath(os.path.join(dirname, os.path.splitext(script)[0] + '.sh'))
        with open(os.path.join(runscript), 'wt') as fout:
            fout.write(bash_script.format(script=scriptpath, options=options).format(roi=name))


if __name__ == "__main__":

    USAGE = "usage: %(prog)s [options]"
    DESCRIPTION = "Clone configuration files and make directory structure"

    PARSER = argparse.ArgumentParser(usage=USAGE, description=DESCRIPTION)

    PARSER.add_argument('--basedir', '-b', default=None, required=True)
    PARSER.add_argument('--roi_set', '-r', default=None, required=True,
                        help='YAML with description of the ROI set')
    PARSER.add_argument('--script', '-s', default=None, required=True,
                        help='A python script to clone')
    PARSER.add_argument('--args', '-a', default=None, required=True,
                        help='Arguments provided to script')

    ARGS = PARSER.parse_args()

    ROI_SET, ROI_FILE = dmp_roi.DMROISet.create_from_yaml(ARGS.roi_set)
    clone_script(ARGS.basedir, ROI_SET, ARGS.script, ARGS.args)
