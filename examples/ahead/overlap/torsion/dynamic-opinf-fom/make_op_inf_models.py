import normaopinf
import normaopinf.opinf
import os
import numpy as np

if __name__ == '__main__':
  settings = {}
  settings['fom-yaml-file'] = "torsion.yaml"
  settings['training-data-directories'] = [os.getcwd()]
  settings['solution-id'] = 1
  settings['model-type'] = 'quadratic'
  settings['forcing'] =  False
  settings['truncation-type'] = 'size'
  settings['boundary-truncation-type'] =  'energy' 
  settings['regularization-parameter'] =  'automatic'
  settings['trial-space-splitting-type'] = 'split'
  settings['acceleration-computation-type'] = 'finite-difference'
  sizes = np.array([2,4,6,8,10,12,14,16,18,20,40,60],dtype=int)
  for i,size in enumerate(sizes):
    settings['truncation-value'] = int(size)
    settings['boundary-truncation-value'] = 1. - 1.e-5 
    settings['operator-name'] = 'quadratic-opinf-operator-rom-dim-' + str(size)
    normaopinf.opinf.make_opinf_model(settings)
