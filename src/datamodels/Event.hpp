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
    RATE_UPDATE_EVENT,
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
    
    Event(int64_t id, EventType event_type);
    virtual double get_next_time() = 0;
    virtual void iterate() = 0;
};

struct PeriodicEvent : Event {
    double initial_time;
    double period;
    int64_t count;
    
    PeriodicEvent(int64_t id, EventType event_type, double initial_time, double period, int64_t count);
    virtual double get_next_time();
    virtual void iterate();
};

struct RateEvent : Event {
    double rate;
    double next_time;
    
    RateEvent(int64_t id, EventType event_type, double rate, double initial_time);
};


/*** Concrete event types ***/

struct CheckpointEvent : PeriodicEvent {
};

struct RateUpdateEvent : Event { };
struct BitingEvent : Event { };
struct ImmigrationEvent : Event { };
struct DeathEvent : Event { };
struct TransitionEvent : Event { };
struct ClearanceEvent : Event { };
struct ImmunityLossEvent : Event { };


/*** CompareEvent functor: compares by get_next_time() ***/ 

struct CompareEvents
{
	bool operator()(Event * ep1, Event * ep2)
	{
		return ep1->get_next_time() < ep2->get_next_time();
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
    RateUpdateEvent * create_rate_update_event();
    BitingEvent * create_biting_event();
    ImmigrationEvent * create_immigration_event();
    DeathEvent * create_death_event();
    TransitionEvent * create_transition_event();
    ClearanceEvent * create_clearance_event();
    ImmunityLossEvent * create_immunity_loss_event();
    
    // Database management
    void load_from_db(sqlite3 * db);
    void resolve_relationships_from_db(sqlite3 * db);
    
    void write_to_db(sqlite3 * db);
    void write_population_table(sqlite3 * db);
    void write_population_hosts_table(sqlite3 * db);
};

} // namespace varmodel

#endif // #ifndef __Event_hpp__
