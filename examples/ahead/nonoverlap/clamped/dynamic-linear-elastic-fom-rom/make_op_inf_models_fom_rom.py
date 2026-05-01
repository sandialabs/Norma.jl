import normaopinf
import normaopinf.opinf
import os
import numpy as np

if __name__ == '__main__':
  
  settings = {}
  settings['fom-yaml-file'] = "../fom-fom/clamped-2.yaml"
  settings['training-data-directories'] = ["../fom-fom"]
  settings['model-type'] = 'linear'
  settings['stop-training-time'] = 'end'
  settings['training-skip-steps'] = 2
  settings['forcing'] =  False

  settings['input-scale'] = 'rms'
  #settings['input-scale'] = 'none'

  settings['truncation-type'] = 'size'
  settings['truncation-value'] = 30
  settings['boundary-truncation-type'] = 'size'
  settings['boundary-truncation-value'] = 4
  
  #settings['truncation-type'] = 'energy'
  #settings['truncation-value'] = 1 - 1e-7
  #settings['boundary-truncation-type'] =  'energy'
  #settings['boundary-truncation-value'] = 1 - 1e-7
  
  #settings['regularization-parameter'] =  1.0e-4
  settings['regularization-parameter'] =  1.0e-5
  settings['trial-space-splitting-type'] = 'split'
  #settings['trial-space-splitting-type'] = 'combined'

  #settings['acceleration-computation-type'] = 'finite-difference'
  settings['acceleration-computation-type'] = 'acceleration-snapshots'

  snapshots_dict = normaopinf.opinf.get_processed_snapshots(settings)
  settings['model-name'] = 'opinf-operator-2'
  normaopinf.opinf.make_opinf_model_from_snapshots_dict(snapshots_dict,settings)

