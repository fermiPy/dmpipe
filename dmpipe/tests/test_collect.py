# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function


def test_collect_link_classes():
    from dmpipe.dm_collect import CollectSED, CollectLimits
    CollectSED.create()
    CollectLimits.create()

def test_collect_sg_classes():
    from dmpipe.dm_collect import CollectSED_SG, CollectLimits_SG, CollectStackedLimits_SG
    CollectSED_SG.create()
    CollectLimits_SG.create()
    CollectStackedLimits_SG.create()


if __name__ == '__main__':
    test_collect_link_classes()
    test_collect_sg_classes()

