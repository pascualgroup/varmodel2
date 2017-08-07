#ifndef state_hpp
#define state_hpp

#include "parameters.hpp"
#include "random.hpp"
#include "Population.hpp"
#include "Host.hpp"
#include "Event.hpp"

namespace varmodel {

extern rng_t * rng;
extern double current_time;

extern PopulationManager * population_manager;
extern HostManager * host_manager;
extern EventManager * event_manager;

extern CheckpointEvent * checkpoint_event;

}; // namespace varmodel

#endif // #ifndef state_hpp
