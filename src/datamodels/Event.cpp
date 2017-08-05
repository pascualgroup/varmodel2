#include "Event.hpp"
#include "parameters.hpp"
#include "random.hpp"
#include "varmodel.hpp"
#include "state.hpp"

namespace varmodel {

Event::Event(uint64_t id, EventType event_type, double time)
    : id(id), event_type(event_type), time(time) {
}


/*** CheckpointEvent ***/

CheckpointEvent::CheckpointEvent(uint64_t id, double time)
    : Event(id, CHECKPOINT_EVENT, time) {
}

void CheckpointEvent::do_event() {
    save_checkpoint();
}


/*** BitingRateUpdateEvent ***/

BitingRateUpdateEvent::BitingRateUpdateEvent(uint64_t id, double time)
    : Event(id, BITING_RATE_UPDATE_EVENT, time) {
}

void BitingRateUpdateEvent::do_event() {
    do_update_biting_rate_event();
}


/*** BitingEvent ***/

BitingEvent::BitingEvent(uint64_t id, Population * population, double time)
    : Event(id, BITING_EVENT, time), population(population) {
}

void BitingEvent::do_event() {
    do_biting_event(population);
}


/*** ImmigrationEvent ***/

ImmigrationEvent::ImmigrationEvent(uint64_t id, Population * population, double time)
    : Event(id, IMMIGRATION_EVENT, time), population(population) {
}

void ImmigrationEvent::do_event() {
    do_immigration_event(population);
}



/*** DeathEvent ***/

DeathEvent::DeathEvent(uint64_t id, Host * host, double time)
    : Event(id, DEATH_EVENT, time), host(host) {
}

void DeathEvent::do_event() {
    do_death_event(host);
}



/*** TransitionEvent ***/

TransitionEvent::TransitionEvent(uint64_t id, Infection * infection, double time)
    : Event(id, TRANSITION_EVENT, time), infection(infection) {
}

void TransitionEvent::do_event() {
    do_transition_event(infection);
}



/*** ClearanceEvent ***/

ClearanceEvent::ClearanceEvent(uint64_t id, Infection * infection, double time)
    : Event(id, TRANSITION_EVENT, time), infection(infection) {
}

void ClearanceEvent::do_event() {
    do_clearance_event(infection);
}


/*** EventManager implementation ***/

EventManager::EventManager() : EventManager(1) {
}

EventManager::EventManager(uint64_t next_id) {
    this->next_id = next_id;
}

CheckpointEvent * EventManager::create_checkpoint_event() {
    return nullptr;
}

BitingRateUpdateEvent * EventManager::create_biting_rate_update_event() {
    return nullptr;
}


} // namespace varmodel
