#!/usr/bin/env python2.7
# ''' Input/Output helper scripts

# This file provided various python functions that use the os and sys
# modules in order to list files, folders etc.

# '''

import os


def list_folder(infolder):
    '''Returns a list of all files in a folder with the full path'''
    return [os.path.join(infolder, f) for f in os.listdir(infolder)]


def readable_dir(infolder):
    '''Tests if the provided infolder is a readable directory'''
    if not os.path.isdir(infolder):
        raise Exception("readable_dir:{0} is not a valid path".format(infolder))
    if os.access(infolder, os.R_OK):
        return infolder
    else:
        raise Exception("readable_dir:{0} is not a readable dir".format(infolder))


def list_files(current_dir):
    '''Recursivly returns all files in a given folder'''
    file_list = []
    for path, subdirs, files in os.walk(current_dir):  # Walk directory tree
        for name in files:
            f = os.path.join(path, name)
            file_list.append(f)
    return file_list


def splitpath(path):
    '''Returns the root and extension of a path'''
    return os.path.splitext(path)