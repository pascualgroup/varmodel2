#include "state.hpp"

namespace varmodel {

double current_time;
rng_t * rng;

PopulationManager * population_manager;
HostManager * host_manager;
EventManager * event_manager;

CheckpointEvent * checkpoint_event;
BitingRateUpdateEvent * biting_rate_update_event;

} // namespace varmodel
