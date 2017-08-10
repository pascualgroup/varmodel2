#include "PopulationManager.hpp"
namespace varmodel {

Population * PopulationManager::create() {
    return create(next_id++);
}

Population * PopulationManager::create(uint64_t id) {
    Population * obj = new Population(id);
    populations.add(obj);
    return obj;
}

Population * PopulationManager::object_for_id(uint64_t id) {
     return collection.object_for_id(id);
}

void PopulationManager::load_from_checkpoint(sqlite3 * db, char const * const table_name) {
    char sql[1024];
    int result = snprintf(sql, sizeof(sql), "SELECT * FROM %s;", table_name);
    assert(result >= 0 && result < sizeof(sql));
    
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, sql, -1, &stmt, NULL);
    while(true) {
        if(sqlite3_step(stmt) != SQLITE_ROW) {
            break;
        }
        Population * obj = create(sqlite3_column_int64(stmt, 0));
        obj->index = sqlite3_column_int64(stmt, 1);
        obj->transmission_count = sqlite3_column_int64(stmt, 2);
    }
    sqlite3_finalize(stmt);
}

void PopulationManager::resolve_relationships(sqlite3 * db, BitingEventManager * biting_event_manager, ImmigrationEventManager * immigration_event_manager) {
    
}

void PopulationManager::save_to_checkpoint(sqlite3 * db, char const * const table_name, char const * const hosts_table_name) {
    char create_sql[8192];
    int result = snprintf(create_sql, sizeof(create_sql),
        "CREATE TABLE %s (id INTEGER, index INTEGER, transmission_count INTEGER, biting_event_id INTEGER, immigration_event_id INTEGER);",
        table_name
    );
    assert(result >= 0 && result < sizeof(create_sql));
    sqlite3_exec(db, create_sql, NULL, NULL, NULL);
    
    char insert_sql[8192];
    int result = snprintf(insert_sql, sizeof(insert_sql),
        "INSERT INTO %s VALUES (?,?,?,?,?);",
        table_name
    );
    assert(result >= 0 && result <= sizeof(insert_sql));
    
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, insert_sql, -1, &stmt, NULL);
    for(Population * obj : objects.as_vector()) {
        sqlite3_bind_int64(stmt, 1, obj->id);
        sqlite3_bind_int64(stmt, 2, obj->index);
        sqlite3_bind_int64(stmt, 3, obj->transmission_count);
        sqlite3_bind_int64(stmt, 4, obj->biting_event->id);
        sqlite3_bind_int64(stmt, 5, obj->immigration_event->id);
        sqlite3_step(stmt);
        sqlite3_reset(stmt); 
    }
    sqlite3_finalize(stmt);
}

} // namespace varmodel
