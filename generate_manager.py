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
    '''Generates a simple database manager for an object type.'''
    args = parse_arguments()
    
    if not os.path.exists(args.type_header_filename):
        print('Header path {} does not exist.'.format(args.header_filename))
        sys.exit(1)
    
    output_dir = os.path.dirname(args.type_header_filename)
    type_header_prefix = os.path.splitext(os.path.basename(args.type_header_filename))[0]
    
    with open(args.type_header_filename) as input_file:
        typeinfo = parse_type(input_file)
    
    print(typeinfo.name)
    print(typeinfo.columns)
    print(typeinfo.reflists)
    
    generate_manager(output_dir, type_header_prefix, typeinfo)

def parse_arguments():
    '''Parses command-line arguments.'''
    
    parser = argparse.ArgumentParser(
        description = 'Generates a database manager from the provided object header file.',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        'type_header_filename', metavar = '<type-header-filename>',
        help = 'Filename of object type header.'
    )
    
    return parser.parse_args()

def parse_type(input_file):
    class TypeInfo:
        def __init__(self, name, columns, reflists):
            self.name = name
            self.columns = columns
            self.reflists = reflists
    
    name = None
    columns = []
    reflists = []
    
    state = 'START'
    for line in input_file:
        tokens = line.strip().split()
        if len(tokens) == 0:
            continue
        if state == 'START':
            if tokens[0] == 'DB_TYPE':
                if len(tokens) != 4 or not tokens[1] == 'struct' or not tokens[3] == '{':
                    print('Cannot parse DB_TYPE line:\n{}'.format(line))
                name = tokens[2]
                state = 'SCANNING'
        if state == 'SCANNING':
            if tokens[0] == '};':
                state = 'END'
            elif tokens[0] == 'DB_FIELD':
                if len(tokens) != 3 or not tokens[2].endswith(';'):
                    print('Cannot parse DB_FIELD line:\n{}'.format(line))
                    sys.exit(1)
                columns.append(('FIELD', tokens[1], tokens[2][:-1]))
            elif tokens[0] == 'DB_REF':
                if len(tokens) != 4 or not tokens[2] == '*' or not tokens[3].endswith(';'):
                    print('Cannot parse DB_REF line:\n{}'.format(line))
                    sys.exit(1)
                columns.append(('REF', tokens[1], tokens[3][:-1]))
            elif tokens[0] == 'DB_REFLIST':
                if len(tokens) != 3 or not tokens[2].endswith(';'):
                    print('Cannot parse DB_REFLIST line:\n{}'.format(line))
                    sys.exit(1)
                reflists.append((tokens[1], tokens[2][:-1]))
    
    return TypeInfo(name, columns, reflists)

manager_header_format = '''#ifndef {manager_type_header_prefix}_hpp
#define {manager_type_header_prefix}_hpp

#include "{object_type_header_prefix}.hpp"
#include "IndexedMap.hpp"

namespace varmodel {{

struct {manager_type} {{
    uint64_t next_id;
    IndexedMap<{object_type}> collection;
    
    // Constructors
    {manager_type}();
    {manager_type}(uint64_t next_id);
    
    // Object management
    {object_type} * create();
    {object_type} * create(uint64_t id);
    {object_type} * object_for_id(uint64_t id);
    
    // Checkpointing methods
    void load_from_checkpoint(sqlite3 * db, char const * const table_name);
    {resolve_relationships_signature};
    {save_to_checkpoint_signature};
}}

}} // namespace varmodel

#endif // #define {manager_type_header_prefix}_hpp
'''

manager_implementation_format = '''#include "{manager_type_header_prefix}.hpp"
namespace varmodel {{

{object_type} * {manager_type}::create() {{
    return create(next_id++);
}}

{object_type} * {manager_type}::create(uint64_t id) {{
    {object_type} * obj = new {object_type}(id);
    populations.add(obj);
    return obj;
}}

{object_type} * {manager_type}::object_for_id(uint64_t id) {{
     return collection.object_for_id(id);
}}

void {manager_type}::load_from_checkpoint(sqlite3 * db, char const * const table_name) {{
    char sql[1024];
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

}}

}}
'''

resolve_relationships_signature_format = \
    'void {prefix}resolve_relationships(sqlite3 * db, char const * const table_name{ref_table_name_args})'

save_to_checkpoint_signature_format = \
    'void {prefix}save_to_checkpoint(sqlite3 * db, char const * const table_name{reflist_table_name_args})'

def generate_manager(output_dir, object_type_header_prefix, typeinfo):
    manager_type_header_prefix = object_type_header_prefix + 'Manager'
    
    object_type = typeinfo.name
    manager_type = object_type + 'Manager'
    
    ref_vars = [c[2] for c in typeinfo.columns if c[0] == 'REF']
    reflist_vars = [rl[1] for rl in typeinfo.reflists]
    
    def format_table_name_args(vars):
        if len(vars) == 0:
            return ''
        return ', ' + ', '.join(['char const * const {}_table_name'.format(var) for var in vars])
    
    def format_manager_header():
        def format_resolve_relationships_signature():
            return resolve_relationships_signature_format.format(
                prefix = '',
                ref_table_name_args = format_table_name_args(ref_vars + reflist_vars)
            )
        
        def format_save_to_checkpoint_signature():
            return save_to_checkpoint_signature_format.format(
                prefix = '',
                reflist_table_name_args = format_table_name_args(reflist_vars)
            )
        
        return manager_header_format.format(
            object_type_header_prefix = object_type_header_prefix,
            object_type = object_type,
            manager_type_header_prefix = manager_type_header_prefix,
            manager_type = manager_type,
            resolve_relationships_signature = format_resolve_relationships_signature(),
            save_to_checkpoint_signature = format_save_to_checkpoint_signature()
        )
    
    manager_header_filename = os.path.join(output_dir, manager_type_header_prefix + '.hpp')
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
    
    def format_manager_implementation():
        def format_resolve_relationships_signature():
            return resolve_relationships_signature_format.format(
                prefix = manager_type + '::',
                ref_table_name_args = format_table_name_args(ref_vars + reflist_vars)
            )
        
        def format_save_to_checkpoint_signature():
            return save_to_checkpoint_signature_format.format(
                prefix = manager_type + '::',
                reflist_table_name_args = format_table_name_args(reflist_vars)
            )
        
        def format_load_column_statements():
            def format_load_column_statement(col_info, index):
                print(col_info, index)
                return 'obj->{name} = sqlite3_column_{sqlite_type}(stmt, {index});'.format(
                    name = col_info[2],
                    sqlite_type = sqlite_type_for_type(col_info[1]),
                    index = index
                )
            
            return '\n        '.join([
                format_load_column_statement(col_info, index + 1) for index, col_info in enumerate(typeinfo.columns)
                if col_info[0] == 'FIELD'
            ])
        
        return manager_implementation_format.format(
            object_type_header_prefix = object_type_header_prefix,
            object_type = object_type,
            manager_type_header_prefix = manager_type_header_prefix,
            manager_type = manager_type,
            resolve_relationships_signature = format_resolve_relationships_signature(),
            save_to_checkpoint_signature = format_save_to_checkpoint_signature(),
            load_column_statements = format_load_column_statements()
        )
    
    manager_implementation_filename = os.path.join(output_dir, manager_type_header_prefix + '.cpp')
    if os.path.exists(manager_implementation_filename):
        print('{} already exists.'.format(manager_implementation_filename))
        sys.exit(1)
    with open(manager_implementation_filename, 'w') as implementation_file:
        implementation_file.write(format_manager_implementation())

if __name__ == '__main__':
    main()
