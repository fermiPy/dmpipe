# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function
 
 
def test_spectral_link_classes():
    from dmpipe.dm_spectral import ConvertCastro, SpecTable, StackLikelihood
    ConvertCastro.create()
    SpecTable.create()
    StackLikelihood.create()

def test_spectral_sg_classes():
    from dmpipe.dm_spectral import ConvertCastro_SG, StackLikelihood_SG
    ConvertCastro_SG.create()
    StackLikelihood_SG.create()

 
if __name__ == '__main__':
    test_spectral_link_classes()
    test_spectral_sg_classes()
