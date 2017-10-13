#!/usr/bin/env python
'''
    make.py
    
    Builds the model code.
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

script_dir = os.path.dirname(__file__)

try:
    basestring = basestring
except NameError:
    basestring = str

def main():
    '''Runs everything.'''
    args = parse_arguments()
    
    copy_sources(args.dest_dir)
    generate_managers(args.dest_dir)
    generate_parameters(args.dest_dir, args.params_file)
    build(args.dest_dir, args.compiler, args.flags)

def parse_arguments():
    '''Parses command-line arguments.'''
    
    parser = argparse.ArgumentParser(
        description = '''Builds the model in the specified directory using provided parameters.
            
            Copies source code into the `src' subdirectory, with 
        ''',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-p', '-P', '--params-file', metavar = '<params-file>',
        default = os.path.join(script_dir, 'parameters-example.py'),
        help = 'Path to parameters file in .py or .json format.'
    )
    parser.add_argument(
        '-d', '-D', '--dest-dir', metavar = '<dest-dir>',
        default = os.path.join(script_dir, 'build'),
        help = 'Path to destination directory for built model.'
    )
    parser.add_argument('-c', '--compiler', metavar = '<compiler>', default = 'c++', help = 'C++ compiler.')
    parser.add_argument('-f', '--flags', metavar = '<flags>', default = '-O2 -g', help = 'Compiler flags.')
    
    return parser.parse_args()

def copy_sources(dst_dirname):
    os.makedirs(dst_dirname)
    
    shutil.copytree(os.path.join(script_dir, 'src'), os.path.join(dst_dirname, 'src'))

def generate_managers(dst_dirname):
    subprocess.Popen([
        os.path.join(script_dir, 'generate_managers.py'),
        '-d',
        dst_dirname
    ]).wait()

def generate_parameters(dst_dirname, params_filename):
    subprocess.Popen([
        os.path.join(script_dir, 'generate_parameters.py'),
        '-d', dst_dirname,
        '-p', params_filename
    ]).wait()

def build(dst_dirname, compiler_cmd, compiler_flags):
    os.makedirs(os.path.join(dst_dirname, 'bin'))
    
    compile_cmd =  compiler_cmd + \
         ' ' + compiler_flags + \
         ' -std=c++11' + \
         ' -lsqlite3' + \
         ' -o bin/varmodel2' + \
         ' generated/managers/*.cpp' + \
         ' src/*.cpp' + \
         ' src/util/*.cpp' + \
         ' -I src' + \
         ' -I src/datamodel' + \
         ' -I src/managers' + \
         ' -I src/util' + \
         ' -I generated' + \
         ' -I generated/managers'
    print(compile_cmd)
    subprocess.Popen(compile_cmd, cwd = dst_dirname, shell = True).wait()

if __name__ == '__main__':
    main()