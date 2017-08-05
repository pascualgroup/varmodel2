#ifndef varmodel_hpp
#define varmodel_hpp

#include "Population.hpp"

namespace varmodel {

void initialize();

void load_checkpoint();
void save_checkpoint();

void run();

void do_update_biting_rate_event();
void do_biting_event(Population * population);
void do_immigration_event(Population * population);
void do_death_event(Host * host);
void do_transition_event(Infection * infection);
void do_clearance_event(Infection * infection);

}; // namespace varmodel

#endif // #ifndef varmodel_hpp
