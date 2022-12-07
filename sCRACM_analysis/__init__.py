name = "sCRACM_analysis"
__version__ = '0.1.1'  # Don't forget to update in setup.py
import matplotlib as mpl
import configparser
import os

mpl.use('TkAgg')

config = configparser.ConfigParser()
ini_file = os.path.join(os.path.dirname(__file__), 'sCRACM.ini')


def generate_ini_file(ini_file_path):
    """Generate a new initiation file.

    Args:
        ini_file_path: full path to initiation file.

    """
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    from tkinter.filedialog import askdirectory
    import pandas as pd

    paths_dict = {'path_database':
                      {'prompt': 'path to main folder containing all databases',
                       'ui_path': None,
                       'template_path': '/Users/username/cloudservice/'
                                        'sCRACM_analysis'},
                  'path_ephys_data':
                      {'prompt': 'path to folder containing analyzed map data. Should'
                                 ' contain .csv files for each recorded neuron',
                       'ui_path': None,
                       'template_path': '/Users/username/cloudservice/'
                                        'insularDatabases/xsg'}}
                  

    files_dict = {
        'filename_cell_db': {'prompt': 'filename of the database; '
                                       ' should be stored in '
                                       'path_database!',
                             'template_file':
                                 'template_key_sCRACM_db.csv'}}

    # create ini file
    cfgfile = open(ini_file, 'w')

    # add settings based on the dicts
    config.add_section('DATABASE_PATHS')
    for key, val_dict in paths_dict.items():
        print('\nPlease select the ' + val_dict['prompt'])
        print('if not present generate the folder and select in gui, '
              'template : ' + val_dict['template_path'])
        Tk().withdraw()
        if key != 'path_database':
            ui_path = askdirectory(initialdir=paths_dict['path_database']
            ['ui_path'])
        else:
            ui_path = askdirectory()
        # store ui in path_dict and initialization file
        paths_dict[key]['ui_path'] = ui_path
        config.set('DATABASE_PATHS', key, ui_path)

        if 'required_folder' in paths_dict[key].keys():
            # test whether all required folder are present otherwise generate
            # that folder
            for req_folder in paths_dict[key]['required_folders']:
                required_folder = os.path.join(paths_dict[key]['ui_path'],
                                               req_folder)
                if not os.path.exists(required_folder):
                    os.mkdir(required_folder, )

    config.add_section('DATABASE_FILENAMES')
    for key, val_dict in files_dict.items():
        existing_db = input('\nIs a ' + key.split('filename_')[-1] + ' present?'
                            + 'Yes: press Return, No: enter any character')
        if existing_db == '':
            print('\nPlease select the ' + val_dict['prompt'])
            Tk().withdraw()
            ui_file_path = askopenfilename(initialdir=paths_dict
            ['path_database']['ui_path'], filetypes=(('csv files', '*.csv'),))
            ui_file = os.path.basename(ui_file_path)
        else:
            print('\nNo database present yet: '
                  'Will create a new database .csv file now')
            db_filename = input('Provide a name for the ' +
                                key.split('filename_')[-1] + 'do not include'
                                                             'extensionname, '
                                                             'eg.: ".csv"')
            print('Please remember to update the descriptors of the key file '
                  'for : ' + db_filename + ' stored in database path : ' +
                  paths_dict['path_database']['ui_path'])
            template_db_path = os.path.join(os.path.dirname(__file__),
                                            val_dict['template_file'])
            template = pd.read_csv(template_db_path)
            db = pd.DataFrame(columns=template.columnName.tolist())
            db_filepath = os.path.join(paths_dict['path_database']['ui_path'],
                                       (db_filename + '.csv'))
            key_db_filepath = os.path.join(
                paths_dict['path_database']['ui_path'],
                ('key_' + db_filename + '.csv'))
            db.to_csv(db_filepath, index=False)
            template.to_csv(key_db_filepath, index=False)
            ui_file = (db_filename + '.csv')
        config.set('DATABASE_FILENAMES', key, ui_file)

    
    config.write(cfgfile)
    cfgfile.close()
    pass


if not os.path.isfile(ini_file):
    generate_ini_file(ini_file_path=ini_file)

#add imports for global functions 
from sCRACM_analysis.sCRACM_global import *
from sCRACM_analysis.mat2py import *
#add imports for helper functions 
