#ifndef __state_hpp__
#define __state_hpp__

#include "parameters.hpp"
#include "random.hpp"
#include "Population.hpp"
#include "Host.hpp"
#include "Event.hpp"

namespace varmodel {

extern rng_t * rng;

extern PopulationManager * population_manager;
extern HostManager * host_manager;
extern EventManager * event_manager;

extern CheckpointEvent * checkpoint_event;
extern BitingRateUpdateEvent * biting_rate_update_event;

}; // namespace varmodel

#endif // #ifndef __state_hpp__
