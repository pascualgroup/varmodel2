//#include "varmodel.hpp"
#include <sqlite3.h>
#include <assert.h>
#include "error_handling.hpp"

int main(int argc, const char * argv[]) {
    program_name = argv[0];
    
    sqlite3_config(SQLITE_CONFIG_LOG, handle_sqlite_error, NULL);
    set_signal_handler();
    
    //sqlite3_exec(nullptr, "CREATE TABLE hello", nullptr, nullptr, nullptr);
    assert(false);
    //varmodel::validate_and_load_parameters();
    //varmodel::load_checkpoint();
    //varmodel::initialize();
    //varmodel::run();
}
