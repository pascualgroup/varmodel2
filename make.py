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

script_dir = os.path.dirname(__file__)

try:
    basestring = basestring
except NameError:
    basestring = str

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
    
    os.makedirs(os.path.join(args.d, 'src'))
    os.makedirs(os.path.join(args.d, 'generated'))
    
    copy_sources(args.d)
    process_template(
        params,
        os.path.join(script_dir, 'templates', 'parameters.hpp.template'),
        os.path.join(args.d, 'generated', 'parameters.hpp')
    )
    for object_type in [
        'Strain', 'Gene',
        'Population', 'Host', 'Infection', 'Immunity'
    ]:
        generate_manager(
            object_type,
            os.path.join(args.d, 'src'),
            os.path.join(args.d, 'generated')
        )

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
        default = os.path.join(script_dir, 'build'),
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
    
    return {name: getattr(params, name) for name in params.__dict__.keys() if not name.startswith('_')}
    
def load_parameters_json(params_filename):
    '''Loads parameters from JSON dictionary.'''
    with open(params_filename) as f:
        return json.load(f)

def copy_sources(build_dir):
    '''Copies source files, inserting parameter values into templates along the way.'''
    
    src_dir = os.path.join(script_dir, 'src')
    build_src_dir = os.path.join(build_dir, 'src')
    
    for filename in os.listdir(src_dir):
        if filename.endswith('.hpp') or filename.endswith('.cpp'):
            src_filename = os.path.join(src_dir, filename)
            dst_filename = os.path.join(build_src_dir, filename)
            shutil.copy2(src_filename, dst_filename)

def process_template(params, src_filename, dst_filename):
    '''Replaces `{varname}' with `varname = value', where value is taken from params.
    
    Value is formatted appropriately using the format_value function.
    '''
    
    # Construct map: varname -> `varname = value` to use as a formatting dictionary
    param_value_map = {}
    for name, value in params.items():
        param_value_map[name] = '{} = {}'.format(name, format_value(name, value))
    
    with open(src_filename) as sf:
        with open(dst_filename, 'w') as df:
            df.write(sf.read().format(**param_value_map))

def format_value(varname, value):
    '''Formats a parameter value for substitution into a template.
    
    The Python type is used to format the value for C++.
    '''
    if isinstance(value, list):
        values = [format_value(varname, val) for val in value]
        values_formatted = '{{{}}}'.format(', '.join(values))
        return values_formatted
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
    print('Error: invalid type {} for value {} for parameter {}'.format(type(value), value, varname))
    sys.exit(1)


### DATABASE MANAGER GENERATION FUNCTIONS ###

def parse_type(object_type, input_filename):
    columns = []
    vectors = []
    reflists_IM = []
    reflists_us = []
    
    with open(input_filename) as input_file:
        state = 'START'
        for line in input_file:
            tokens = line.strip().split('//')[0].split()
            print(tokens, state)
            if len(tokens) == 0:
                continue
            if state == 'START':
                if len(tokens) == 3 and tokens[0] == 'struct' and tokens[1] == object_type and  tokens[2] == '{':
                    state = 'SCANNING'
            if state == 'SCANNING':
                if tokens[0] == '};':
                    state = 'END'
                elif len(tokens) == 2:
                    if tokens[0] == 'uint64_t' or tokens[0] == 'int64_t' or tokens[0] == 'double':
                        if tokens[1].endswith(';'):
                            columns.append(('FIELD', tokens[0], tokens[1][:-1]))
                    elif tokens[0] == 'std::vector<uint64_t>' or tokens[0] == 'std::vector<int64_t>' or tokens[0] == 'std::vector<double>':
                        if tokens[1].endswith(';'):
                            field_type = tokens[0][len('std::vector')+1:-1]
                            vectors.append((field_type, tokens[1][:-1]))
                    elif tokens[0].startswith('IndexedSet') \
                        and tokens[0][len('IndexedSet')] == '<' and tokens[0][-1] == '>' \
                        and tokens[1].endswith(';'):
                            map_object_type = tokens[0][len('IndexedSet')+1:-1]
                            reflists_IM.append((map_object_type, tokens[1][:-1]))
                elif len(tokens) == 3:
                    if tokens[1] == '*' and tokens[2].endswith(';'):
                        columns.append(('REF', tokens[0], tokens[2][:-1]))
                    elif tokens[0].startswith('std::unordered_set<') \
                        and tokens[1] == '*>' and tokens[2].endswith(';'):
                        set_object_type = tokens[0][len('std::unordered_set')+1:]
                        reflists_us.append((set_object_type, tokens[2][:-1]))
    
    return columns, vectors, reflists_IM, reflists_us

with open(os.path.join(script_dir, 'templates', 'Manager.hpp.template')) as f:
    manager_header_format = f.read()

with open(os.path.join(script_dir, 'templates', 'Manager.cpp.template')) as f:
    manager_implementation_format = f.read()

with open(os.path.join(script_dir, 'templates', 'create_reflist_block_IndexedSet.cpp.template')) as f:
    create_reflist_block_IndexedSet_format = f.read()

with open(os.path.join(script_dir, 'templates', 'create_reflist_block_unordered_set.cpp.template')) as f:
    create_reflist_block_unordered_set_format = f.read()

with open(os.path.join(script_dir, 'templates', 'create_vector_block.cpp.template')) as f:
    create_vector_block_format = f.read()

resolve_relationships_signature_format = \
    'void {prefix}resolve_relationships(sqlite3 * db{ref_manager_args})'

def generate_manager(object_type, src_dir, dst_dir):
    src_filename = os.path.join(src_dir, '{}.hpp'.format(object_type))
    
    manager_type = object_type + 'Manager'
    print('Generating {}...'.format(manager_type))
    
    columns, vectors, reflists_IM, reflists_us = parse_type(object_type, src_filename)
    
    ref_cols = [c for c in columns if c[0] == 'REF']
    ref_vars = [c[2] for c in columns if c[0] == 'REF']
    reflist_vars = [rl[1] for rl in reflists_IM + reflists_us]
    
    print('Vector variables: {}'.format(json.dumps(vectors)))
    print('Reference variables: {}'.format(json.dumps(ref_vars)))
    print('Reference lists: {}'.format(json.dumps(reflist_vars)))
    
    def format_forward_declarations():
        return '\n'.join([
            'struct {}Manager;'.format(c[1]) for c in ref_cols
        ])
    
    def format_table_name_args(vars):
        if len(vars) == 0:
            return ''
        return ', ' + ', '.join(['char const * const {}_table_name'.format(var) for var in vars])
    
    def format_manager_args():
        if len(ref_cols) == 0:
            return ''
        return ', ' + ', '.join(['{}Manager * {}_manager'.format(c[1], c[2]) for c in ref_cols])
    
    def format_manager_header():
        def format_resolve_relationships_signature():
            return resolve_relationships_signature_format.format(
                prefix = '',
                ref_manager_args = format_manager_args()
            )
        
        return manager_header_format.format(
            object_type = object_type,
            manager_type = manager_type,
            forward_declarations = format_forward_declarations(),
            resolve_relationships_signature = format_resolve_relationships_signature()
        )
    
    manager_header_filename = os.path.join(dst_dir, manager_type + '.hpp')
    if os.path.exists(manager_header_filename):
        print('{} already exists.'.format(manager_header_filename))
        sys.exit(1)
    with open(manager_header_filename, 'w') as header_file:
        header_file.write(format_manager_header())
    
    def sqlite_type_for_type(type_str):
        if type_str.startswith('uint') or type_str.startswith('int'):
            return 'int64'
        elif type_str == 'double':
            return 'double'
        elif type_str == 'std::string':
            return 'text'
        return 'int64'
     
    def sql_type_for_type(type_str):
        if type_str.startswith('uint') or type_str.startswith('int'):
            return 'INTEGER'
        elif type_str == 'double':
            return 'REAL'
        elif type_str == 'std::string':
            return 'TEXT'
        return 'INTEGER'
    
    def format_manager_implementation():
        def format_manager_includes():
            return '\n'.join([
                '#include "{}Manager.hpp"'.format(c[0])
                for c in reflists_IM + reflists_us
            ])
        
        def format_object_includes():
            return '\n'.join([
                '#include "{}.hpp"'.format(c[0])
                for c in reflists_IM + reflists_us
            ] + [
                '#include "{}.hpp"'.format(c[1])
                for c in ref_cols
            ])
        
        def format_resolve_relationships_signature():
            return resolve_relationships_signature_format.format(
                prefix = manager_type + '::',
                ref_manager_args = format_manager_args(),
            )
        
        def format_load_column_statements():
            def format_load_column_statement(col_info, index):
                return 'obj->{name} = sqlite3_column_{sqlite_type}(stmt, {index});'.format(
                    name = col_info[2],
                    sqlite_type = sqlite_type_for_type(col_info[1]),
                    index = index
                )
            
            return '\n        '.join([
                format_load_column_statement(col_info, index + 1) for index, col_info in enumerate(columns)
                if col_info[0] == 'FIELD'
            ])
        
        def format_sql_create_columns():
            if len(columns) == 0:
                return ''
            
            def format_create_column(col_info):
                if col_info[0] == 'FIELD':
                    return '{} {}'.format(col_info[2], sql_type_for_type(col_info[1]))
                else:
                    assert col_info[0] == 'REF'
                    return '{}_id INTEGER'.format(col_info[2])
            
            return ', ' + ', '.join([format_create_column(col_info) for col_info in columns])
        
        def format_sql_insert_qmarks():
            if len(columns) == 0:
                return ''
            
            return ',' + ','.join(['?'] * len(columns))
        
        def format_bind_column_statements():
            def format_bind_column_statement(col_info, index):
                if col_info[0] == 'FIELD':
                    return 'sqlite3_bind_{}(stmt, {}, obj->{});'.format(
                        sqlite_type_for_type(col_info[1]),
                        index,
                        col_info[2]
                    )
                else:
                    assert col_info[0] == 'REF'
                    return 'sqlite3_bind_int64(stmt, {}, obj->{}->id);'.format(
                        index,
                        col_info[2]
                    )
            
            return '\n        '.join([
                format_bind_column_statement(col_info, index + 2) for index, col_info in enumerate(columns)
            ])
        
        def format_create_reflist_blocks_IndexedSet():
            def create_reflist_block(reflist_spec):
                ref_type, reflist_var = reflist_spec
                return create_reflist_block_IndexedSet_format.format(
                    object_type = object_type,
                    ref_type = ref_type,
                    reflist_var = reflist_var
                )
            
            return '\n        '.join([
                create_reflist_block(reflist_spec) for reflist_spec in reflists_IM
            ])
        
        def format_create_reflist_blocks_unordered_set():
            def create_reflist_block(reflist_spec):
                ref_type, reflist_var = reflist_spec
                return create_reflist_block_unordered_set_format.format(
                    object_type = object_type,
                    ref_type = ref_type,
                    reflist_var = reflist_var
                )
            
            return '\n        '.join([
                create_reflist_block(reflist_spec) for reflist_spec in reflists_us
            ])
        
        def format_create_vector_blocks():
            def create_vector_block(spec):
                value_type, vector_var = spec
                return create_vector_block_format.format(
                    object_type = object_type,
                    sql_type = sql_type_for_type(value_type),
                    vector_var = vector_var,
                    sqlite3_bind_type = sqlite_type_for_type(value_type)
                )
            
            return '\n        '.join([
                create_vector_block(spec) for spec in vectors
            ])
        
        return manager_implementation_format.format(
            object_type = object_type,
            manager_type = manager_type,
            manager_includes = format_manager_includes(),
            object_includes = format_object_includes(),
            resolve_relationships_signature = format_resolve_relationships_signature(),
            load_column_statements = format_load_column_statements(),
            sql_create_columns = format_sql_create_columns(),
            sql_insert_qmarks = format_sql_insert_qmarks(),
            bind_column_statements = format_bind_column_statements(),
            create_vector_blocks = format_create_vector_blocks(),
            create_reflist_blocks_IndexedSet = format_create_reflist_blocks_IndexedSet(),
            create_reflist_blocks_unordered_set = format_create_reflist_blocks_unordered_set()
        )
    
    manager_implementation_filename = os.path.join(dst_dir, manager_type + '.cpp')
    if os.path.exists(manager_implementation_filename):
        print('{} already exists.'.format(manager_implementation_filename))
        sys.exit(1)
    with open(manager_implementation_filename, 'w') as implementation_file:
        implementation_file.write(format_manager_implementation())

if __name__ == '__main__':
    main()
