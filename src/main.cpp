#include "varmodel.hpp"

int main(int argc, const char * argv[]) {
    varmodel::validate_and_load_parameters();
    varmodel::load_checkpoint();
    //varmodel::initialize();
    //varmodel::run();
}
