#ifndef util_hpp
#define util_hpp

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

namespace varmodel {

template<typename T>
struct HashVector
{
	size_t operator()(std::vector<T> const & vec) const
	{
		if(vec.size() == 0) {
			return 0;
		}
		size_t hash_val = vec.size();
		for(size_t i = 0; i < vec.size(); i++) {
			hash_val ^= _hash(vec[i]) + 0x9e3779b9 + (hash_val << 6) + (hash_val >> 2);
		}
		return hash_val;
	}
    
	std::hash<T> _hash;
};

std::vector<size_t> makeRange(size_t from, size_t to);
std::vector<size_t> makeRange(size_t size);

double add(std::vector<double> const & vec);
std::vector<double> addCumulative(std::vector<double> const & vec);

bool file_exists(std::string const & filename);

} // namespace varmodel

#endif // #ifndef util_hpp
