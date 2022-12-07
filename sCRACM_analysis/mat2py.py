import os
import re
import numpy as np
from scipy import io

"""
for pulling info out of matlab .xsg, .mat, and .m file and reading into python dicts 
Functions: 
-helpers: 
    _check_keys: 
    -todict: recursive 
-Execution functions: 
    - loadmat: loads matlab structures as py dicts, improvement upon scipy.io.loadmat. 
    adapted from: https://stackoverflow.com/questions/11955000/how-to-preserve-matlab-struct-when-accessing-in-python
    - mstruct2pydict: Convert a Matlab *.m file __that defines a Matlab struct()__
        into a Python dict.  Each struct.FIELD corresponds to a dict[KEY]. Contributed by 
        Mike Muniak (muniak@ohsu.edu)

"""

def _check_keys( dict):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dict:
        if isinstance(dict[key], io.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, io.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def loadmat(filename):
    """
    this function should be called instead of direct scipy.io .loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    data = io.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


def mstruct2pydict(path):
    """ Convert a Matlab *.m file __that defines a Matlab struct()__
        into a Python dict.  Each struct.FIELD corresponds to a dict[KEY].
        
        True/False strings are interpreted as boolean.
        All valid numbers are interpreted as floats.
        Anything else is a string.
        1D and 2D matrices are interpreted as NumPy ndarrays.
        
        Things that have not been dealt with yet (but could be).
        3D (and beyond) matrices.
        Matlab cells {}.
        Nested struct fields.
        
        v.2021.01.20
        muniak@ohsu.edu
    """
    
    # Does file exist?
    assert os.path.isfile(path), 'Cannot find the file at %s!' % path
    
    # Open/read file contents.
    with open(path) as f:
        mfile = f.read()
    
    # Replace Matlab line continuations (...).
    mfile = re.sub('\.\.\.\s*\n', ' ', mfile)
    
    # Get the base struct name used in the m-file script.
    root = re.match('^[\s\n]*function\s*(\w+)\s*\=.*$', mfile, re.M).group(1)
    
    # Use regex to pull out individual script lines defining struct fields.  The trick here is that 
    # matrices will span multiple lines (e.g., contain line returns), so have to catch that.
    rx_lines = re.compile('^[\s\n]*ROOT\.([^%](.*(\[[^\]]*\]))?[^;\n]*)[\;\s\n]*$'.replace('ROOT', root), re.M)
    
    # Recursive fxn to determine val type: Array, Boolean, Float, or String.
    # Note: if you require some numbers to be integers instead of floats, you will have to convert afterwards.
    def cast_val(v):
        # If value is encased in [] then must be an array/matrix.
        if v.startswith('[') and v.endswith(']'):
            # Semicolons [;] not encased by quotes [''] or [""] are interpreted as new lines/rows.
            v = re.sub(r'\;(?=(?:(?:[^"\']*["\']){2})*[^"\']*\Z)', '\n', v[1:-1])
            # Commas [,] not encased by quotes [''] or [""] are interpreted as delimiters.
            v = re.sub(r'\,(?=(?:(?:[^"\']*["\']){2})*[^"\']*\Z)', ' ', v)
            # Newline chars mean we have a 2D matrix -- process each entry in each line.
            if '\n' in v:
                return np.array([[cast_val(s2) for s2 in s1.split()] for s1 in v.strip().split('\n')])
            # Otherwise we have a 1D matrix -- process each entry.
            else:
                return np.array([cast_val(s) for s in v.split()])
        # Or is it a boolean?
        elif v.lower() == 'true':
            return True
        elif v.lower() == 'false':
            return False
        # Not a boolean...
        else:
            # Is it a number?  Remove surrounding commas if present.
            try:
                return float(v.strip(','))
            # Must be a string.  Remove surrounding quotes if present.
            except ValueError:
                return v.strip(r'\'\"')
    
    out = dict()
    
    # Loop through entries and process.
    for line in rx_lines.finditer(mfile):
        # Split by first '=' to get key/val pairs.  ('=' could also be present later in value string.)
        key, val = line.group(1).split('=', 1)
        # Strip outer whitespace and store key/value pair with correct casting.
        out[key.strip()] = cast_val(val.strip())
    
    # Return dict.
    return out