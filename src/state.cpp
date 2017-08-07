#include "state.hpp"

namespace varmodel {

rng_t * rng;
double current_time;

PopulationManager * population_manager;
HostManager * host_manager;
EventManager * event_manager;

CheckpointEvent * checkpoint_event;

} // namespace varmodel
