import numpy as np
  
import copy

LOG_ERROR=1
LOG_INFO=2
LOG_WARNING=3
LOG_ALLINFO=4


class loggable:
  """
  Simple class allowing for logging messages through it's debug() method.
  """
  log_level=LOG_WARNING
  _last_message = None
  _last_count = 0
  log_names = ['===', 'EE:', 'II:', 'WW:', 'AA:']
  
  @classmethod
  def debug(cls,text,level):
    if loggable.log_level >= level:
      if text == cls._last_message:
        cls._last_count += 1
      elif cls._last_count > 1:
        print('repeated ({} times)'.format(cls._last_count))
      else:
        print(cls.log_names[level], text)
        cls._last_message = text
        cls._last_count = 1


def metro_dist( v1, v2 ):
  """
  Helper function used in fixing "jumped ions" and problematic polarizations.
  """
  return abs(v1 - v2).max()

def mat2str(m, mode=' % 7.4f'):
  if len(m.shape) == 1:
    m = np.array([m])
  if len(m.shape) > 2:
    m = m.reshape( m.shape[0], -1 )
      
  return '\n'.join([' '.join([mode%col for col in row]) for row in np.array(m) ] )+'\n'

def iterate_all_indices(shape,shape_min=None):
  """
  Generator function for looping through all possible combinations of the values in given ranges.
  """
  if shape_min is None:
    shape_min = [0]*len(shape)

  for i in range(int(shape_min[0]),int(shape[0])):
    if len(shape) > 1:
      for k in iterate_all_indices(shape[1:], shape_min[1:]):
        yield [i]+k
    else:
      yield [i]    

def iterate_shifts():
  """
  Generator function for iterating over all neighbours of 3D cube.
  """
  for dx,dy,dz in iterate_all_indices([3,3,3]):
    yield (dx-1,dy-1,dz-1)
