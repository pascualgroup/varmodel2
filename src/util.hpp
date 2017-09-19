#ifndef util_hpp
#define util_hpp

#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>

namespace varmodel {

std::vector<size_t> makeRange(size_t from, size_t to);
std::vector<size_t> makeRange(size_t size);

double add(std::vector<double> const & vec);
std::vector<double> addCumulative(std::vector<double> const & vec);

bool file_exists(std::string const & filename);

} // namespace varmodel

#endif // #ifndef util_hpp
