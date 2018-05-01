# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import absolute_import, division, print_function


def test_pipeline_chains():
    from dmpipe.pipeline import PipelineData, PipelineSim, PipelineRandom, Pipeline
    PipelineData.create()
    PipelineSim.create()
    PipelineRandom.create()
    Pipeline.create()


if __name__ == '__main__':
    test_pipeline_chains()

