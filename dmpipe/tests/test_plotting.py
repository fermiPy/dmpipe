# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function


def test_plotting_link_classes():
    from dmpipe.dm_plotting import PlotDMSpectra, PlotCastro, PlotLimits, PlotDM
    PlotDMSpectra.create()
    PlotCastro.create()
    PlotLimits.create()
    PlotDM.create()

def test_plotting_sg_classes():
    from dmpipe.dm_plotting import PlotCastro_SG, PlotLimits_SG, PlotStackedLimits_SG, PlotDM_SG, PlotStackedDM_SG
    PlotCastro_SG.create()
    PlotLimits_SG.create()
    PlotStackedLimits_SG.create()
    PlotDM_SG.create()
    PlotStackedDM_SG.create()

if __name__ == '__main__':
    test_plotting_link_classes()
    test_plotting_sg_classes()
