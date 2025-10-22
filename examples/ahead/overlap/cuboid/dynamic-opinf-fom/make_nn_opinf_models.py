import normaopinf
import normaopinf.opinf
import nnopinf
import nnopinf.training
import os
import numpy as np

if __name__ == '__main__':
    settings = {}
    settings['fom-yaml-file'] = 'cuboid-2.yaml'
    settings['training-data-directories'] = ['/Users/ejparis/Research/ahead/tmp/Norma.jl/examples/ahead/overlap/cuboid/dynamic-opinf-fom']
    settings['model-type'] = 'neural-network'
    settings['stop-training-time'] = 100000.0
    settings['training-skip-steps'] = 1
    settings['forcing'] = False
    settings['truncation-type'] = 'energy'
    settings['boundary-truncation-type'] = 'energy'
    settings['regularization-parameter'] = [0.0005, 0.005, 0.05]
    settings['model-name'] = 'opinf-operator'
    settings['truncation-value'] = 1. - 1.e-5
    settings['boundary-truncation-value'] = 1. - 1.e-5
    settings['trial-space-splitting-type'] = 'split'
    settings['acceleration-computation-type'] = 'finite-difference'
    settings['neural-network-training-settings'] = {'model-name': 'opinf-operator', 'output-path': 'ml-models', 'num-epochs': 15000, 'batch-size': 500, 'learning-rate': 0.001, 'weight-decay': 1e-08, 'lr-decay': 0.9999, 'print-training-output': True, 'epoch': 5000, 'resume': False}
    settings['ensemble-size'] = 5
    snapshots_dict = normaopinf.opinf.get_processed_snapshots(settings)
    normaopinf.opinf.make_opinf_model_from_snapshots_dict(snapshots_dict, settings)
