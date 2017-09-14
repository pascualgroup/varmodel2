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
        os.path.join(script_dir, 'parameters.hpp.template'),
        os.path.join(args.d, 'generated', 'parameters.hpp')
    )
    for object_type in ['Population', 'Host']:
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
    '''Replaces `{{varname}}' with `varname = value', where value is taken from params.
    
    Value is formatted appropriately using the format_value function.
    '''
    with open(src_filename) as sf:
        with open(dst_filename, 'w') as df:
            for line in sf:
                pieces1 = line.split('{{')
                if len(pieces1) == 2:
                    pieces2 = pieces1[1].split('}}')
                    if len(pieces2) == 2:
                        varname = pieces2[0]
                        
                        if not hasattr(params, pieces2[0]):
                            print('Error: parameter {} not present'.format(varname))
                            sys.exit(1)
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
    reflists = []
    
    with open(input_filename) as input_file:
        state = 'START'
        for line in input_file:
            tokens = line.strip().split()
            if len(tokens) == 0:
                continue
            if state == 'START':
                if len(tokens) == 3 and tokens[0] == 'struct' and tokens[1] == object_type and  tokens[-2] == '{':
                    state = 'SCANNING'
            if state == 'SCANNING':
                if tokens[0] == '};':
                    state = 'END'
                elif len(tokens) == 2:
                    if tokens[0] == 'uint64_t' or tokens[0] == 'int64_t' or tokens[0] == 'double':
                        if tokens[1].endswith(';'):
                            columns.append(('FIELD', tokens[0], tokens[1][:-1]))
                    elif tokens[0].startswith('IndexedMap') \
                        and tokens[0][len('IndexedMap')] == '<' and tokens[0][-1] == '>' \
                        and tokens[1].endswith(';'):
                            map_object_type = tokens[0][len('IndexedMap')+1:-1]
                            reflists.append((map_object_type, tokens[1][:-1]))
                elif len(tokens) == 3 and tokens[1] == '*' and tokens[2].endswith(';'):
                    columns.append(('REF', tokens[0], tokens[2][:-1]))
    
    return columns, reflists

manager_header_format = '''#ifndef {manager_type}_hpp
#define {manager_type}_hpp

#include "{object_type}.hpp"
#include "IndexedMap.hpp"

namespace varmodel {{

struct {manager_type} {{
    uint64_t next_id;
    IndexedMap<{object_type}> collection;
    
    // Object management
    {object_type} * create();
    {object_type} * create(uint64_t id);
    {object_type} * object_for_id(uint64_t id);
    
    // Checkpointing methods
    void load_from_checkpoint(sqlite3 * db, char const * const table_name);
    {resolve_relationships_signature};
    {save_to_checkpoint_signature};
}};

}} // namespace varmodel

#endif // #define {manager_type}_hpp
'''

manager_implementation_format = '''#include "{manager_type}.hpp"
namespace varmodel {{

{object_type} * {manager_type}::create() {{
    return create(next_id++);
}}

{object_type} * {manager_type}::create(uint64_t id) {{
    {object_type} * obj = new {object_type}(id);
    collection.add(obj);
    return obj;
}}

{object_type} * {manager_type}::object_for_id(uint64_t id) {{
     return collection.object_for_id(id);
}}

void {manager_type}::load_from_checkpoint(sqlite3 * db, char const * const table_name) {{
    char sql[8192];
    int result = snprintf(sql, sizeof(sql), "SELECT * FROM %s;", table_name);
    assert(result >= 0 && result < sizeof(sql));
    
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    while(true) {{
        if(sqlite3_step(stmt) != SQLITE_ROW) {{
            break;
        }}
        {object_type} * obj = create(sqlite3_column_int64(stmt, 0));
        {load_column_statements}
    }}
    sqlite3_finalize(stmt);
}}

{resolve_relationships_signature} {{
    
}}

{save_to_checkpoint_signature} {{
    char create_sql[8192];
    int result = snprintf(create_sql, sizeof(create_sql),
        "CREATE TABLE %s (id INTEGER{sql_create_columns});",
        table_name
    );
    assert(result >= 0 && result < sizeof(create_sql));
    sqlite3_exec(db, create_sql, NULL, NULL, NULL);
    
    char insert_sql[8192];
    result = snprintf(insert_sql, sizeof(insert_sql),
        "INSERT INTO %s VALUES (?{sql_insert_qmarks});",
        table_name
    );
    assert(result >= 0 && result <= sizeof(insert_sql));
    
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, insert_sql, -1, &stmt, NULL);
    for({object_type} * obj : collection.as_vector()) {{
        sqlite3_bind_int64(stmt, 1, obj->id);
        {bind_column_statements}
        sqlite3_step(stmt);
        sqlite3_reset(stmt); 
    }}
    sqlite3_finalize(stmt);
}}

}} // namespace varmodel
'''

resolve_relationships_signature_format = \
    'void {prefix}resolve_relationships(sqlite3 * db{ref_manager_args})'

save_to_checkpoint_signature_format = \
    'void {prefix}save_to_checkpoint(sqlite3 * db, char const * const table_name{reflist_table_name_args})'

def generate_manager(object_type, src_dir, dst_dir):
    src_filename = os.path.join(src_dir, '{}.hpp'.format(object_type))
    
    manager_type = object_type + 'Manager'
    
    columns, reflists = parse_type(object_type, src_filename)
    
    ref_cols = [c for c in columns if c[0] == 'REF']
    ref_vars = [c[2] for c in columns if c[0] == 'REF']
    reflist_vars = [rl[1] for rl in reflists]
    
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
        
        def format_save_to_checkpoint_signature():
            return save_to_checkpoint_signature_format.format(
                prefix = '',
                reflist_table_name_args = format_table_name_args(reflist_vars)
            )
        
        return manager_header_format.format(
            object_type = object_type,
            manager_type = manager_type,
            resolve_relationships_signature = format_resolve_relationships_signature(),
            save_to_checkpoint_signature = format_save_to_checkpoint_signature()
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
        def format_resolve_relationships_signature():
            return resolve_relationships_signature_format.format(
                prefix = manager_type + '::',
                ref_manager_args = format_manager_args(),
            )
        
        def format_save_to_checkpoint_signature():
            return save_to_checkpoint_signature_format.format(
                prefix = manager_type + '::',
                reflist_table_name_args = format_table_name_args(reflist_vars)
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
        
        return manager_implementation_format.format(
            object_type = object_type,
            manager_type = manager_type,
            resolve_relationships_signature = format_resolve_relationships_signature(),
            save_to_checkpoint_signature = format_save_to_checkpoint_signature(),
            load_column_statements = format_load_column_statements(),
            sql_create_columns = format_sql_create_columns(),
            sql_insert_qmarks = format_sql_insert_qmarks(),
            bind_column_statements = format_bind_column_statements()
        )
    
    manager_implementation_filename = os.path.join(dst_dir, manager_type + '.cpp')
    if os.path.exists(manager_implementation_filename):
        print('{} already exists.'.format(manager_implementation_filename))
        sys.exit(1)
    with open(manager_implementation_filename, 'w') as implementation_file:
        implementation_file.write(format_manager_implementation())

if __name__ == '__main__':
    main()
