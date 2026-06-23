import normaopinf
import normaopinf.opinf
import nnopinf
import nnopinf.training
import os
import numpy as np

if __name__ == '__main__':
  settings = {}
  settings['fom-yaml-file'] = "notched-cylinder-1.yaml"
  settings['training-data-directories'] = [os.getcwd()]
  settings['model-type'] = 'linear'
  settings['stop-training-time'] = 'end'
  settings['training-skip-steps'] = 1
  settings['forcing'] =  False
  settings['truncation-type'] = 'size'
  settings['truncation-value'] = 3
  settings['boundary-truncation-value'] = 99.999
  settings['boundary-truncation-type'] =  'energy'
  settings['regularization-parameter'] =  [1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 1.0e-9, 1.0e-10]
  settings['trial-space-splitting-type'] = 'combined'
  settings['acceleration-computation-type'] = 'acceleration-snapshots'
  snapshots_dict = normaopinf.opinf.get_processed_snapshots(settings)
  settings['model-name'] = 'linear-opinf-operator-1'
  normaopinf.opinf.make_opinf_model_from_snapshots_dict(snapshots_dict,settings)

