#include "varmodel.hpp"

int main(int argc, const char * argv[]) {
    varmodel::validate_and_load_parameters();
    varmodel::initialize();
    varmodel::save_checkpoint();
    varmodel::run();
}
