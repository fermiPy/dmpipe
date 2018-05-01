# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function


def test_analysis_link_classes():
    from dmpipe.target_analysis import PrepareTargets, AnalyzeROI, AnalyzeSED
    PrepareTargets.create()
    AnalyzeROI.create()
    AnalyzeSED.create()

def test_analysis_sg_classes():
    from dmpipe.target_analysis import AnalyzeROI_SG, AnalyzeSED_SG
    AnalyzeROI_SG.create()
    AnalyzeSED_SG.create()


if __name__ == '__main__':
    test_analysis_link_classes()
    test_analysis_sg_classes()

