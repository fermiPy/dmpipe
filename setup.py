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
        'dmpipe-analyze-roi-sg = dmpipe.target_analysis:main_roi_batch',
        'dmpipe-analyze-sed-sg = dmpipe.target_analysis:main_sed_batch',
        'dmpipe-prepare-targets = dmpipe.scripts:main',
        'dmpipe-convert-castro = dmpipe.dm_spectral:main_single',
        'dmpipe-convert-castro-sg = dmpipe.dm_spectral:main_batch',
        'dmpipe-stack-likelihood = dmpipe.scripts.stack_dm_likelihood:main',
    ]},
    install_requires=[
        'numpy >= 1.6.1',
        'astropy >= 1.2.1',
        'matplotlib >= 1.5.0',
        'scipy >= 0.14',
        #'fermipy >= 0.13.0',
        'fermipy == 0.13.5+8.g2849a.dirty',
        'pyyaml',
        'healpy',
        'dmsky'
    ],
    extras_require=dict(
        all=[],
    ),
)
