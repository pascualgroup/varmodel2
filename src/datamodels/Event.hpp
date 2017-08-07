#ifndef Event_hpp
#define Event_hpp

#include <sqlite3.h>
#include <vector>
#include "IndexedPriorityQueue.hpp"

namespace varmodel {

/*** Forward type definitions ***/

struct Population;
struct Host;
struct Infection;

/*** EventType enum used for introspection ***/

enum EventType {
    CHECKPOINT_EVENT,
    BITING_EVENT,
    IMMIGRATION_EVENT,
    DEATH_EVENT,
    TRANSITION_EVENT,
    CLEARANCE_EVENT,
    IMMUNITY_LOSS_EVENT
};


/*** Abstract event superclass ***/

struct Event {
    const uint64_t id;
    const EventType event_type;
    double time;
    
    Event(uint64_t id, EventType event_type, double time);
    virtual void do_event() = 0;
};


/*** Concrete event types ***/

struct CheckpointEvent : Event {
    CheckpointEvent(uint64_t id, double time);
    virtual void do_event();
};

struct BitingEvent : Event {
    BitingEvent(uint64_t id, Population * population, double time);
    virtual void do_event();
    
    Population * population;
};

struct ImmigrationEvent : Event {
    ImmigrationEvent(uint64_t id, Population * population, double time);
    virtual void do_event();
    
    Population * population;
};

struct DeathEvent : Event {
    DeathEvent(uint64_t id, Host * host, double time);
    virtual void do_event();
    
    Host * host;
};

struct TransitionEvent : Event {
    TransitionEvent(uint64_t id, Infection * infection, double time);
    virtual void do_event();
    
    Infection * infection; 
};

struct ClearanceEvent : Event {
    ClearanceEvent(uint64_t id, Infection * infection, double time);
    virtual void do_event();
    
    Infection * infection;
};

struct ImmunityLossEvent : Event {
};


/*** CompareEvent functor: compares by get_next_time() and then by event ID as a tiebreaker ***/ 

struct CompareEvents
{
	bool operator()(Event * e1, Event * e2)
	{
        if(e1->time == e2->time) {
            return e1->id < e2->id;
        }
		return e1->time < e2->time;
	}
};


/*** EventManager declaration ***/

struct EventManager {
    uint64_t next_id;
    IndexedPriorityQueue<Event *, nullptr, std::hash<Event *>, CompareEvents> events;
    
    EventManager();
    EventManager(uint64_t next_id);
    
    // Object management
    CheckpointEvent * create_checkpoint_event();
    BitingEvent * create_biting_event(
        Population * population,
        double time
    );
    ImmigrationEvent * create_immigration_event(
        Population * population,
        double time
    );
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

#endif // #ifndef Event_hpp
