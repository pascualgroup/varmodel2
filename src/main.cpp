#include "varmodel.hpp"

int main(int argc, const char * argv[]) {
    varmodel::initialize();
    varmodel::save_checkpoint();
    varmodel::run();
}
