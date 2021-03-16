# Licensed under a 3-clause BSD style license - see LICENSE.rst



def test_collect_link_classes():
    from dmpipe.dm_collect import CollectLimits
    CollectLimits.create()

def test_collect_sg_classes():
    from dmpipe.dm_collect import CollectLimits_SG, CollectStackedLimits_SG
    CollectLimits_SG.create()
    CollectStackedLimits_SG.create()


if __name__ == '__main__':
    test_collect_link_classes()
    test_collect_sg_classes()

