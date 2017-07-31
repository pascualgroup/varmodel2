#include "Event.hpp"
#include "parameters.hpp"
#include "random.hpp"
#include "varmodel.hpp"
#include "state.hpp"

namespace varmodel {

Event::Event(
    int64_t id, EventType event_type, double time
) : id(id), event_type(event_type) {
    this->time = time;
}


/*** CheckpointEvent ***/

CheckpointEvent::CheckpointEvent(int64_t id, double time) : Event(id, CHECKPOINT_EVENT, time) {
}

void CheckpointEvent::do_event() {
    save_checkpoint();
}


/*** BitingRateUpdateEvent ***/

BitingRateUpdateEvent::BitingRateUpdateEvent(int64_t id, double time) : Event(id, BITING_RATE_UPDATE_EVENT, time) {
    update_biting_rate();
}

void BitingRateUpdateEvent::do_event() {
    update_biting_rate();
}

/*** EventManager implementation ***/

EventManager::EventManager() : EventManager(1) {
}

EventManager::EventManager(int64_t next_id) {
    this->next_id = next_id;
}

CheckpointEvent * EventManager::create_checkpoint_event() {
    return nullptr;
}

BitingRateUpdateEvent * EventManager::create_biting_rate_update_event() {
    return nullptr;
}


} // namespace varmodel
