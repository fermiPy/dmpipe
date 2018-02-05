#!/usr/bin/env python
from setuptools import setup, find_packages
import versioneer

setup(
    name='dmpipe',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    author='Matthew Wood',
    author_email='mdwood@slac.stanford.edu',
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
        'dmpipe-analyze-roi = dmpipe.target_analysis:main_roi_single',
        'dmpipe-analyze-sed = dmpipe.target_analysis:main_sed_single',
        'dmpipe-prepare-targets = dmpipe.target_analysis:main_prepare_targets',
        'dmpipe-analyze-roi-sg = dmpipe.target_analysis:main_roi_batch',
        'dmpipe-analyze-sed-sg = dmpipe.target_analysis:main_sed_batch',
        'dmpipe-convert-castro = dmpipe.dm_spectral:main_convert_single',
        'dmpipe-convert-castro-sg = dmpipe.dm_spectral:main_convert_batch',
        'dmpipe-stack-likelihood = dmpipe.dm_spectral:main_stack_likelihood',
        'dmpipe-spec-table = dmpipe.dm_spectral:main_spec_table',
        'dmpipe-pipeline-dsph = dmpipe.pipeline_dsph:main_chain',
        'dmpipe-plot-dm = dmpipe.scripts.plot_castro_dm:main_single', 
        'dmpipe-plot-dm-sg = dmpipe.scripts.plot_castro_dm:main_batch', 
        'dmpipe-plot-castro = dmpipe.scripts.plot_castro:main_single', 
        'dmpipe-plot-castro-sg = dmpipe.scripts.plot_castro:main_batch',
        'dmpipe-calc-sed-limits = dmpipe.scripts.calc_sed_limits:main', 
    ]},
    install_requires=[
        'numpy >= 1.6.1',
        'astropy >= 1.2.1',
        'matplotlib >= 1.5.0',
        'scipy >= 0.14',
        'fermipy >= 0.14.0',
        'pyyaml',
        'healpy',
        'dmsky'
    ],
    extras_require=dict(
        all=[],
    ),
)
