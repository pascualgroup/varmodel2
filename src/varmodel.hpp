#ifndef varmodel_hpp
#define varmodel_hpp

#include "Population.hpp"

namespace varmodel {

void initialize();

void load_checkpoint();
void save_checkpoint();

void run();

}; // namespace varmodel

#endif // #ifndef varmodel_hpp
