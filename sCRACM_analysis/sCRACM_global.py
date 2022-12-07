"""
This module contains all variables used as globals 

For use in other modules first import:
from sCRACM_analysis import sCRACM_global
and refer in code as sCRACM_global.path_db
"""

import configparser
import os 

config = configparser.ConfigParser()
ini_file = os.path.join(os.path.dirname(__file__), 'sCRACM.ini')
config.read(ini_file)

path_database = config.get('DATABASE_PATHS', 'path_database')
path_ephys_data = config.get('DATABASE_PATHS', 'path_ephys_data')

filename_cell_db = config.get('DATABASE_FILENAMES', 'filename_cell_db')
