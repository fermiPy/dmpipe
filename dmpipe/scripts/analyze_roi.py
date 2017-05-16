import os
import sys
import copy
import numpy as np
import itertools
import argparse

from fermipy.utils import init_matplotlib_backend

init_matplotlib_backend()

from fermipy.gtanalysis import GTAnalysis
from fermipy.catalog import Catalog3FGL

def main():
        
    usage = "usage: %(prog)s [config file]"
    description = "Run fermipy analysis chain."
    parser = argparse.ArgumentParser(usage=usage,description=description)

    parser.add_argument('--config', default = 'sample_config.yaml')
    parser.add_argument('--source', default = None)

    args = parser.parse_args()
    gta = GTAnalysis(args.config,logging={'verbosity' : 3},
                     fileio={'workdir_regex' : '\.xml$|\.npy$'})


    model0 = { 'SpatialModel' : 'PointSource', 'Index' : 1.5 }
    model1 = { 'SpatialModel' : 'PointSource', 'Index' : 2.0 }
    model2 = { 'SpatialModel' : 'PointSource', 'Index' : 2.7 }

    src_name = gta.config['selection']['target']
    
    gta.setup(overwrite=True)
    gta.free_sources(False)
    gta.print_roi()
    gta.optimize()
    gta.print_roi()

    exclude = []
    
    # Localize all point sources
    for s in sorted(gta.roi.sources, key=lambda t: t['ts'],reverse=True):
#    for s in gta.roi.sources:

        if not s['SpatialModel'] == 'PointSource':
            continue
        if s['offset_roi_edge'] > -0.1:
            continue

        if s.name in exclude:
            continue
        if not '3FGL' in s.name:
            continue
        if s.name == src_name:
            continue
        
        gta.localize(s.name,nstep=5,dtheta_max=0.5,update=True,
                     prefix='base', make_plots=True)

    gta.optimize()
    gta.print_roi()
    
    gta.write_roi('base_roi',make_plots=True)

    exclude = [src_name]
    if not 'carina_2' in exclude:
        exclude += ['carina_2']
    if not 'carina_3' in exclude:
        exclude += ['carina_3']
    
    gta.tsmap('base',model=model0,make_plots=True,exclude=exclude)
    gta.residmap('base',model=model0,make_plots=True,exclude=exclude)
    gta.tsmap('base',model=model1,make_plots=True,exclude=exclude)
    gta.residmap('base',model=model1,make_plots=True,exclude=exclude)
    gta.tsmap('base',model=model2,make_plots=True,exclude=exclude)
    gta.residmap('base',model=model2,make_plots=True,exclude=exclude)

    gta.find_sources(sqrt_ts_threshold=5.0)
    gta.optimize()
    gta.print_roi()
    gta.print_params()
        
    gta.free_sources(skydir=gta.roi.skydir,distance=1.0, pars='norm')
    gta.fit()
    gta.print_roi()
    gta.print_params()
    
    gta.write_roi('fit0_roi',make_plots=True)
    
    m = gta.tsmap('fit0',model=model0,make_plots=True,exclude=exclude)
    gta.plotter.make_tsmap_plots(m, gta.roi,
                                 zoom=2,suffix='tsmap_zoom')    
    gta.residmap('fit0',model=model0,make_plots=True,exclude=exclude)
    gta.tsmap('fit0',model=model1,make_plots=True,exclude=exclude)
    gta.plotter.make_tsmap_plots(m, gta.roi,
                                 zoom=2,suffix='tsmap_zoom')
    gta.residmap('fit0',model=model1,make_plots=True,exclude=exclude)
    gta.tsmap('fit0',model=model2,make_plots=True,exclude=exclude)
    gta.plotter.make_tsmap_plots(m, gta.roi,
                                 zoom=2,suffix='tsmap_zoom')
    gta.residmap('fit0',model=model2,make_plots=True,exclude=exclude)

    gta.sed(src_name, prefix='fit0', make_plots=True, free_radius=1.0)

    gta.free_source(src_name)
    gta.fit(reoptimize=True)
    gta.print_roi()
    gta.print_params()
    
    gta.write_roi('fit1_roi',make_plots=True)
    

if __name__ == '__main__':

    main()
