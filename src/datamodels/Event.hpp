#ifndef __Event_hpp__
#define __Event_hpp__

#include <sqlite3.h>
#include <vector>
#include "Host.hpp"
#include "IndexedPriorityQueue.hpp"

namespace varmodel {

/*** EventType enum used for introspection ***/

enum EventType {
    CHECKPOINT_EVENT,
    BITING_RATE_UPDATE_EVENT,
    BITING_EVENT,
    IMMIGRATION_EVENT,
    DEATH_EVENT,
    TRANSITION_EVENT,
    CLEARANCE_EVENT,
    IMMUNITY_LOSS_EVENT
};


/*** Abstract event types ***/

struct Event {
    const int64_t id;
    const EventType event_type;
    double time;
    
    Event(int64_t id, EventType event_type, double time);
    virtual void do_event() = 0;
};


/*** Concrete event types ***/

struct CheckpointEvent : Event {
    CheckpointEvent(int64_t id, double time);
    virtual void do_event();
};

struct BitingRateUpdateEvent : Event {
    BitingRateUpdateEvent(int64_t id, double time);
    virtual void do_event();
};

struct BitingEvent : Event { };
struct ImmigrationEvent : Event { };
struct DeathEvent : Event { };
struct TransitionEvent : Event { };
struct ClearanceEvent : Event { };
struct ImmunityLossEvent : Event { };


/*** CompareEvent functor: compares by get_next_time() and then by event type ***/ 

struct CompareEvents
{
	bool operator()(Event * e1, Event * e2)
	{
        if(e1->time == e2->time) {
            return e1->event_type < e2->event_type;
        }
		return e1->time < e2->time;
	}
};


/*** EventManager declaration ***/

struct EventManager {
    int64_t next_id;
    IndexedPriorityQueue<Event *, nullptr, std::hash<Event *>, CompareEvents> events;
    
    EventManager();
    EventManager(int64_t next_id);
    
    // Object management
    CheckpointEvent * create_checkpoint_event();
    BitingRateUpdateEvent * create_biting_rate_update_event();
    
    BitingEvent * create_biting_event();
    ImmigrationEvent * create_immigration_event();
    DeathEvent * create_death_event();
    TransitionEvent * create_transition_event();
    ClearanceEvent * create_clearance_event();
    ImmunityLossEvent * create_immunity_loss_event();
    
    void update_event(Event * event);
    void delete_event(Event * event);
    
    // Database management
    void load_from_db(sqlite3 * db);
    void resolve_relationships_from_db(sqlite3 * db);
    
    void write_to_db(sqlite3 * db);
    void write_population_table(sqlite3 * db);
    void write_population_hosts_table(sqlite3 * db);
};

} // namespace varmodel

#endif // #ifndef __Event_hpp__
