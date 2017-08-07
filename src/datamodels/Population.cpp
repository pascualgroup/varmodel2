#include "state.hpp"

namespace varmodel {

/*** Population class implementation ***/

Population::Population(uint64_t id) : id(id) {
}


/*** PopulationManager clas implementation ***/

PopulationManager::PopulationManager() : PopulationManager(1) { }

PopulationManager::PopulationManager(uint64_t next_id) : next_id(next_id) {
}

Population * PopulationManager::create() {
    return create(next_id++);
}

Population * PopulationManager::create(uint64_t id) {
    Population * pop = new Population(id);
    populations.add(pop);
    return pop;
}

Population * PopulationManager::population_for_id(uint64_t id) {
    return populations.object_for_id(id);
}

void PopulationManager::load_from_db(sqlite3 * db) {
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, "SELECT * FROM population;", -1, &stmt, NULL);
    while(true) {
        if(sqlite3_step(stmt) != SQLITE_ROW) {
            break;
        }
        uint64_t id = sqlite3_column_int64(stmt, 0);
        Population * pop = create(id);
        pop->index = sqlite3_column_int64(stmt, 1);
        pop->transmission_count = sqlite3_column_double(stmt, 2);
    }
    sqlite3_finalize(stmt);
}

void PopulationManager::resolve_relationships_from_db(sqlite3 * db) {
    resolve_population_hosts_table(db);
}

void PopulationManager::resolve_population_hosts_table(sqlite3 * db) {
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, "SELECT * FROM population_hosts;", -1, &stmt, NULL);
    while(true) {
        if(sqlite3_step(stmt) != SQLITE_ROW) {
            break;
        }
        uint64_t pop_id = sqlite3_column_int64(stmt, 0);
        Population * pop = population_for_id(pop_id);
        
        uint64_t host_id = sqlite3_column_int64(stmt, 1);
        Host * host = host_manager->host_for_id(host_id);
        
        pop->hosts.add(host);
    }
    
    sqlite3_finalize(stmt);
}

void PopulationManager::write_to_db(sqlite3 * db) {
    write_population_table(db);
    write_population_hosts_table(db); 
}

void PopulationManager::write_population_table(sqlite3 * db) {
    sqlite3_exec(db,
        "CREATE TABLE population (id INTEGER, index INTEGER, transmission_count INTEGER);",
        NULL, NULL, NULL
    );
    
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, "INSERT INTO population VALUES (?,?,?);", -1, &stmt, NULL);
    for(Population * pop : populations.as_vector()) {
        sqlite3_bind_int64(stmt, 1, pop->id);
        sqlite3_bind_int64(stmt, 2, pop->index);
        sqlite3_bind_int64(stmt, 3, pop->transmission_count);
        sqlite3_step(stmt);
        sqlite3_reset(stmt); 
    }
    sqlite3_finalize(stmt);
}

void PopulationManager::write_population_hosts_table(sqlite3 * db) {
    sqlite3_exec(db,
        "CREATE TABLE population_hosts (population_id INTEGER, host_id INTEGER);",
        NULL, NULL, NULL
    );
    
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, "INSERT INTO population_hosts VALUES (?,?);", -1, &stmt, NULL);
    for(Population * pop : populations.as_vector()) {
        for(Host * host : pop->hosts.as_vector()) {
            sqlite3_bind_int64(stmt, 1, pop->id);
            sqlite3_bind_int64(stmt, 2, host->id);
            sqlite3_step(stmt);
            sqlite3_reset(stmt);
        }
    }
    sqlite3_finalize(stmt);
}

} // namespace varmodel
