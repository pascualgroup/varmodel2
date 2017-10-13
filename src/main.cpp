//#include "varmodel.hpp"
#include <sqlite3.h>
#include <assert.h>
#include "stack_traces.hpp"

void print_traceback_and_abort();
void handle_sqlite_error(void * p_arg, int error_code, const char * msg);

int main(int argc, const char * argv[]) {
    program_name = argv[0];
    
    sqlite3_config(SQLITE_CONFIG_LOG, handle_sqlite_error, NULL);
    set_signal_handler();
    
    //sqlite3_exec(nullptr, "CREATE TABLE hello", nullptr, nullptr, nullptr);
    //assert(false);
    cause_calamity();
    //varmodel::validate_and_load_parameters();
    //varmodel::load_checkpoint();
    //varmodel::initialize();
    //varmodel::run();
}

void print_traceback_and_abort() {
    posix_print_stack_trace();
    exit(1);
}

void handle_sqlite_error(void * p_arg, int error_code, const char * msg) {
    fprintf(stderr, "SQLite error occurred! Aborting.\n");
    fprintf(stderr, "(%d) %s\n", error_code, msg);
    print_traceback_and_abort();
}

