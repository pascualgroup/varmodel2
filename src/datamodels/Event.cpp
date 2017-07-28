#include "Event.hpp"

namespace varmodel {

Event::Event(
    int64_t id, EventType event_type
) : id(id), event_type(event_type) {
}

PeriodicEvent::PeriodicEvent(
    int64_t id, EventType event_type, double initial_time, double period, int64_t count
) : Event(id, event_type) {
    this->initial_time = initial_time;
    this->period = period;
    this->count = count;
}

double PeriodicEvent::get_next_time() {
    return initial_time + period * count;
}

void PeriodicEvent::iterate() {
    count += 1;
}


/*** EventManager implementation ***/

EventManager::EventManager() : EventManager(1) {
}

EventManager::EventManager(int64_t next_id) {
    this->next_id = next_id;
}


} // namespace varmodel
