import normaopinf
import normaopinf.opinf
import nnopinf
import nnopinf.training
import os
import numpy as np

if __name__ == '__main__':
    settings = {}
    settings['fom-yaml-file'] = 'cuboid-2.yaml'
    settings['training-data-directories'] = [os.getcwd()]
    settings['model-type'] = 'neural-network'
    settings['stop-training-time'] = 100000.0
    settings['training-skip-steps'] = 1
    settings['forcing'] = False
    settings['truncation-type'] = 'energy'
    settings['boundary-truncation-type'] = 'energy'
    settings['regularization-parameter'] = [0.0005, 0.005, 0.05]
    settings['model-name'] = 'opinf-operator'
    settings['truncation-value'] = 0.999999
    settings['boundary-truncation-value'] = 0.999999
    settings['trial-space-splitting-type'] = 'combined'
    settings['acceleration-computation-type'] = 'finite-difference'
    settings['neural-network-training-settings'] = {'model-name': 'opinf-operator', 'output-path': 'ml-models', 'print-training-output': True, 'resume': False, 'GN-final-layer': False, 'GN-num-layers': 0, 'GN-final-layer-epoch-frequency': 0, 'GN-final-layer-damping': 0.0, 'GN-verbose': False, 'optimizer': {'method': 'ADAM', 'num-epochs': 10000, 'batch-size': 500, 'learning-rate': 0.001, 'weight-decay': 1e-08, 'lr-decay': 0.9999, 'LBFGS-acceleration': {'enabled': True, 'acceleration-epoch-frequency': 1000, 'acceleration-iterations': 50}}}
    settings['architecture'] = {'model-structure': 'PsdLagrangianOperator', 'n-hidden-layers': 2, 'n-neurons-per-layer': 'auto'}
    settings['ensemble-size'] = 3
    snapshots_dict = normaopinf.opinf.get_processed_snapshots(settings)
    normaopinf.opinf.make_opinf_model_from_snapshots_dict(snapshots_dict, settings)
