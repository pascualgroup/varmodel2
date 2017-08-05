#include "Host.hpp"

namespace varmodel {

/*** Host class implementation ***/

Host::Host(uint64_t id) : id(id) {
}


/*** HostManager class implementation ***/

HostManager::HostManager() : HostManager(1) { }

HostManager::HostManager(uint64_t next_id) {
    this->next_id = next_id;
}

Host * HostManager::create() {
    return create(next_id++);
}

Host * HostManager::create(uint64_t id) {
    Host * host = new Host(id);
    hosts.add(host);
    return host;
}

Host * HostManager::host_for_id(uint64_t id) {
    return hosts.object_for_id(id);
}

void HostManager::load_from_db(sqlite3 * db) {
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, "SELECT * FROM host;", -1, &stmt, NULL);
    while(true) {
        if(sqlite3_step(stmt) != SQLITE_ROW) {
            break;
        }
        uint64_t id = sqlite3_column_int64(stmt, 0);
        create(id);
    }
    sqlite3_finalize(stmt);
}

void HostManager::resolve_relationships_from_db(sqlite3 * db) {
}

void HostManager::write_to_db(sqlite3 * db) {
    write_host_table(db); 
}

void HostManager::write_host_table(sqlite3 * db) {
    sqlite3_exec(db,
        "CREATE TABLE host (id INTEGER, transmission_count INTEGER);",
        NULL, NULL, NULL
    );
    
    sqlite3_stmt * stmt = NULL;
    sqlite3_prepare_v2(db, "INSERT INTO population VALUES (?,?);", -1, &stmt, NULL);
    for(Host * host : hosts.as_vector()) {
        sqlite3_bind_int64(stmt, 1, host->id);
        sqlite3_step(stmt);
        sqlite3_reset(stmt); 
    }
    sqlite3_finalize(stmt);
}


} // namespace varmodel
