#include "varmodel.hpp"
#include <sqlite3.h>
#include <assert.h>
#include "error_handling.hpp"

int main(int argc, const char * argv[]) {
    program_name = argv[0];
    
    sqlite3_config(SQLITE_CONFIG_LOG, handle_sqlite_error, NULL);
    register_signal_handler();
    
    varmodel::run();
}
