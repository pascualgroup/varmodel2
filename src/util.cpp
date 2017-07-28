#include <cassert>
#include "util.hpp"

using namespace std;

namespace varmodel {

std::vector<size_t> makeRange(size_t size)
{
	return makeRange(0, size);
}

std::vector<size_t> makeRange(size_t from, size_t to)
{
	assert(to >= from);
	size_t size = to - from;
	vector<size_t> range(size);
	for(size_t i = 0; i < size; i++)
	{
		range[i] = from + i;
	}
	return range;
}

double add(std::vector<double> const & vec)
{
	double sum = 0.0;
	for(double val : vec) {
		sum += val;
	}
	return sum;
}

std::vector<double> addCumulative(std::vector<double> const & vec)
{
	vector<double> cumSum(vec);
	for(size_t i = 1; i < cumSum.size(); i++) {
		cumSum[i] += cumSum[i-1];
	}
	return cumSum;
}

} // namespace varmodel
