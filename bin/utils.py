'''
Utilities for input and output handling
'''

import os


def get_basename(file_name):
    '''
    Retrieves the basename of the file to allow flexible output names
    '''
    # retrives the basename of the input file
    if "/" in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    
    return(basename)


def check_output_dir(output_dir):
    """
    checks if output dir exists, if not create output dir
    checks that output dir is formatted correctly with '/' on the end - if not, adds '/'
    """
    # checks if output directory exists or creates one
    if not (os.path.exists(output_dir)):
        os.mkdir(output_dir)
        print('Creating output directory')

    if output_dir.endswith("/"):
        output_dir = output_dir
    else:
        output_dir = output_dir + "/"
    # checks if output directory exists after attempting to create one
    if not (os.path.exists(output_dir)):
        print('Output dir could not be created - check permissions & create output dir manually if necessary')
    
    return(output_dir)