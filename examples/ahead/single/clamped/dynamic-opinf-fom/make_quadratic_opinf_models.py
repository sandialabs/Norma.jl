import normaopinf
import normaopinf.opinf
import nnopinf
import nnopinf.training
import os
import numpy as np

if __name__ == '__main__':
  settings = {}
  settings['fom-yaml-file'] = "clamped_rounded_square.yaml"
  settings['training-data-directories'] = [os.getcwd()]
  settings['model-type'] = 'quadratic'
  settings['stop-training-time'] = 'end'
  settings['training-skip-steps'] = 1
  settings['forcing'] =  False
  settings['truncation-type'] = 'size'
  settings['boundary-truncation-type'] =  'size'
  settings['truncation-value'] = 30 
  settings['boundary-truncation-value'] = 4
  settings['regularization-parameter'] =  1e-11
  settings['trial-space-splitting-type'] = 'split'
  settings['acceleration-computation-type'] = 'acceleration-snapshots'
  snapshots_dict = normaopinf.opinf.get_processed_snapshots(settings)
  settings['model-name'] = 'quadratic-opinf-operator'
  normaopinf.opinf.make_opinf_model_from_snapshots_dict(snapshots_dict,settings)

