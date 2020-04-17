"""Module containing some usefull mathematical recepies."""
import numpy as np
import math

def find_nearest(array,value):
    """Return the nearest array entry to the given value"""
    idx = np.searchsorted(array, value)
    if math.fabs(value - array[idx-1]) < math.fabs(value - array[idx]):
        result = array[idx-1]
    else:
        result = array[idx]
    
    return result

