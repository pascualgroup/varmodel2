#ifndef IndexedMap_hpp
#define IndexedMap_hpp

#include <vector>
#include <unordered_map>
#include <cassert>

namespace varmodel {

template<typename T>
struct IndexedMap {
    std::vector<T *> vec;
    std::unordered_map<int64_t, size_t> id_index_map; 
    
    void add(T * obj) {
        assert(id_index_map.find(obj->id) == id_index_map.end());
        vec.push_back(obj);
        id_index_map[obj->id] = vec.size() - 1;
    }
    
    void remove(T * obj) {
        size_t index = id_index_map[obj->id];
        id_index_map.erase(obj->id);
        if(index < vec.size() - 1) {
            vec[index] = vec.back();
            id_index_map[vec[index]->id] = index;
        }
        vec.pop_back();
    }
    
    T * object_for_id(int64_t id) {
        return object_at_index(id_index_map[id]);
    }
    
    T * object_at_index(size_t index) {
        return vec[index];
    }
    
    std::vector<T *> & as_vector() {
        return vec;
    }
};

}; // namespace varmodel

#endif // #ifndef IndexedMap_hpp
