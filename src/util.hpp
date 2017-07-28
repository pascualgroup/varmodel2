#ifndef __util_hpp__
#define __util_hpp__

#include <vector>
#include <unordered_set>
#include <unordered_map>

namespace varmodel {

std::vector<size_t> makeRange(size_t from, size_t to);
std::vector<size_t> makeRange(size_t size);

double add(std::vector<double> const & vec);
std::vector<double> addCumulative(std::vector<double> const & vec);

} // namespace varmodel

#endif
