#!/usr/bin/env python
'''
    make.py
    
    Builds the model code after substituting parameters into template files.
    Run `make.py -h' to see documentation of command-line options.
'''

import sys
if sys.version_info < (2,7):
    print('This script requires Python 2.7 or later.')
    sys.exit(1)

import os
import shutil
import importlib
import json
import argparse
import subprocess

try:
    basestring = basestring
except NameError:
    basestring = str

script_dir = os.path.dirname(__file__)

def main():
    '''Runs everything.'''
    args = parse_arguments()
    
    if os.path.exists(args.d):
        if not os.path.isdir(args.d):
            print('Destination path {} exists but is not directory.'.format(args.d))
            sys.exit(1)
    else:
        os.makedirs(args.d)
    
    params = load_parameters(args.p)
    
    copy_sources(params, args.d)

def parse_arguments():
    '''Parses command-line arguments.'''
    
    parser = argparse.ArgumentParser(
        description = '''Builds the model in the specified directory using provided parameters.
            
            Copies source code into the `src' subdirectory, with 
        ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-d', metavar = '<destination>',
        default = 'build',
        help = 'Destination directory for built model.'
    )
    parser.add_argument(
        '-p', metavar = '<params-file>',
        default = os.path.join(script_dir, 'parameters-example.py'),
        help = '''
            Filename for model parameters, either a Python module (.py) or JSON (.json).
            These parameters are used to generate `parameters.hpp' from `parameters.hpp.template'
            (or any other source files ending in `.template'),
            filling in all entries of the form {{variable_name}} with the corresponding value defined in the file.
            Use parameter names as variable names in a Python module, or as dictionary keys in a JSON file.
        '''
    )
    parser.add_argument('-c', metavar = '<compiler>', default = 'c++', help = 'C++ compiler.')
    parser.add_argument('-f', metavar = '<flags>', default = '-O2', help = 'Compiler flags.')
    
    return parser.parse_args()

def load_parameters(params_filename):
    '''Loads parameters from Python module or JSON file.'''
    
    base, ext = os.path.splitext(params_filename)
    if ext == '.py':
        return load_parameters_python(params_filename)
    elif ext == '.json':
        return load_parameters_json(params_filename)
    
    print('{}: unknown file type'.format(params_filename))
    sys.exit(1)

def load_parameters_python(params_filename):
    '''Loads parameters from Python module using importlib.'''
    
    old_sys_path = sys.path
    sys.path.insert(0, os.path.dirname(params_filename))
    module_name = os.path.splitext(os.path.basename(params_filename))[0]
    params = importlib.import_module(module_name)
    sys.path = old_sys_path
    
    return params
    
def load_parameters_json(params_filename):
    '''Loads parameters from JSON file into object attributes.'''
    
    class JSONObject(object):
        def __init__(self, pairs):
            for key, value in pairs:
                setattr(self, key, value)
    
    with open(params_filename) as f:
        return json.load(f, object_pairs_hook = JSONObject)

def copy_sources(params, build_dir):
    '''Copies source files, inserting parameter values into templates along the way.'''
    
    code_dir = os.path.join(script_dir, 'src')
    build_code_dir = os.path.join(build_dir, 'src')
    
    if os.path.exists(build_code_dir):
        print('{} already exists.'.format(build_code_dir))
        sys.exit(1)
    os.makedirs(build_code_dir)
    
    for dirpath, dirnames, filenames in os.walk(code_dir):
        build_dirpath = os.path.join(build_code_dir, os.path.relpath(dirpath, code_dir))
        try:
            os.makedirs(build_dirpath)
        except:
            pass
        
        for filename in filenames:
            if filename.startswith('.'):
                continue
            elif filename.endswith('.template'):
                process_template(params, dirpath, build_dirpath, filename)
            else:
                src_filename = os.path.join(dirpath, filename)
                dst_filename = os.path.join(build_dirpath, filename)
                shutil.copy2(src_filename, dst_filename)

def process_template(params, dirpath, build_dirpath, filename):
    '''Replaces `{{varname}}' with `varname = value', where value is taken from params.
    
    Value is formatted appropriately using the format_value function.
    '''
    src_filename = os.path.join(dirpath, filename)
    dst_filename = os.path.join(build_dirpath, os.path.splitext(filename)[0])
    
    with open(src_filename) as sf:
        with open(dst_filename, 'w') as df:
            for line in sf:
                pieces1 = line.split('{{')
                if len(pieces1) == 2:
                    pieces2 = pieces1[1].split('}}')
                    if len(pieces2) == 2:
                        varname = pieces2[0]
                        
                        if not hasattr(params, pieces2[0]):
                            value = None
                        else:
                            value = getattr(params, pieces2[0])
                        line_subs = '{}{} = {}{}'.format(
                            pieces1[0],
                            varname, format_value(varname, value),
                            pieces2[1]
                        )
                        df.write(line_subs)
                    else:
                        df.write(line)
                else:
                    df.write(line)

def format_value(varname, value):
    '''Formats a parameter value for substitution into a template.
    
    The Python type is used to format the value for C++.
    For simplicity: can find only a single instance of a variable substitution per line;
    this is good enough for us.
    '''
    if isinstance(value, list):
        print(value)
        values = [format_value(varname, val) for val in value]
        print(values)
        print(', '.join(values))
        values_formatted = '{{{}}}'.format(', '.join(values))
        print(values_formatted)
        return values_formatted
    elif value is None:
        return 'NULL'
    elif value is False:
        return 'false'
    elif value is True:
        return 'true'
    elif isinstance(value, int):
        return str(value)
    elif isinstance(value, float):
        return json.dumps(value) # It's possible this will break in weird situations
    elif isinstance(value, basestring):
        return '"{}"'.format(value)
    print('Invalid type {} for value {} for parameter {}'.format(type(value), value, varname))
    sys.exit(1)

if __name__ == '__main__':
    main()
