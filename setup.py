#!/usr/bin/env python
from setuptools import setup, find_packages
import versioneer

setup(
    name='dmpipe',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author='Eric Charles',
    author_email='echarles@slac.stanford.edu',
    description='Pipeline Scripts for LAT DM Analysis',
    license='BSD',
    packages=find_packages(),
    include_package_data=True,
    url="https://github.com/fermiPy/dmpipe",
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Astronomy',
        'Development Status :: 4 - Beta',
    ],
    scripts=[],
    entry_points={'console_scripts': [
        'dmpipe-prepare-targets = dmpipe.dm_prepare:PrepareTargets.main',
        'dmpipe-convert-castro = dmpipe.dm_spectral:ConvertCastro.main',
        'dmpipe-convert-castro-sg = dmpipe.dm_spectral:ConvertCastro_SG.main',
        'dmpipe-stack-likelihood = dmpipe.dm_spectral:StackLikelihood.main',
        'dmpipe-stack-likelihood-sg = dmpipe.dm_spectral:StackLikelihood_SG.main',
        'dmpipe-collect-limits = dmpipe.dm_collect:CollectLimits.main',
        'dmpipe-collect-limits-sg = dmpipe.dm_collect:CollectLimits_SG.main',
        'dmpipe-collect-stacked-limits-sg = dmpipe.dm_collect:CollectStackedLimits_SG.main',
        'dmpipe-spec-table = dmpipe.dm_spectral:SpecTable.main',
        'dmpipe-pipeline-dsph = dmpipe.pipeline_dsph:main_chain',
        'dmpipe-plot-dm-spectra = dmpipe.dm_plotting:PlotDMSpectra.main', 
        'dmpipe-plot-dm = dmpipe.dm_plotting:PlotDM.main', 
        'dmpipe-plot-dm-sg = dmpipe.dm_plotting:PlotDM_SG.main', 
        'dmpipe-plot-stacked-dm-sg = dmpipe.dm_plotting:PlotStackedDM_SG.main', 
        'dmpipe-plot-limits = dmpipe.dm_plotting:PlotLimits.main', 
        'dmpipe-plot-mles = dmpipe.dm_plotting:PlotMLEs.main', 
        'dmpipe-plot-limits-sg = dmpipe.dm_plotting:PlotLimits_SG.main', 
        'dmpipe-plot-stacked-limits-sg = dmpipe.dm_plotting:PlotStackedLimits_SG.main', 
        'dmpipe-plot-control-limits-sg = dmpipe.dm_plotting:PlotControlLimits_SG.main', 
        'dmpipe-plot-control-mles-sg = dmpipe.dm_plotting:PlotControlMLEs_SG.main', 
        'dmpipe-plot-final-limits-sg = dmpipe.dm_plotting:PlotFinalLimits_SG.main', 
        'dmpipe-calc-sed-limits = dmpipe.scripts.calc_sed_limits:main', 
    ]},
    install_requires=[
        'numpy',
        'pyyaml'
        'astropy',
        'matplotlib',
        'scipy',
        'fermipy',
        'healpy',
        'dmsky'
    ],
    extras_require=dict(
        all=[],
    ),
)
