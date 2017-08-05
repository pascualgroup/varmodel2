#ifndef IndexedPriorityQueue_hpp
#define IndexedPriorityQueue_hpp

#include <vector>
#include <unordered_map>
#include <functional>
#include <cassert>
#include <iostream>

namespace varmodel {

template<typename T, T nullValue, typename Hash, typename Compare>
class IndexedPriorityQueue;

template<typename T, T nullValue, typename Hash, typename Compare=std::less<T>>
class IndexedPriorityQueueProbe
{
public:
	IndexedPriorityQueueProbe(std::vector<T> const & initHeap) :
		queue(initHeap, true)
	{
	}
	
	IndexedPriorityQueue<T, nullValue, Hash, Compare> queue;
	
	void update(T obj) { queue.update(obj); }
	
	bool add(T const & obj) { return queue.add(obj); }
	bool remove(T const & obj) { return queue.remove(obj); }
	
	size_t getHead() { return queue.getHead(); }
	
	void buildHeap() { queue.buildHeap(); }
	bool verifyHeap() { return queue.verifyHeap(); }
	bool heapifyUp(size_t i) { return queue.heapifyUp(i); }
	bool heapifyDown(size_t i) { return queue.heapifyDown(i); }
	void swap(size_t i1, size_t i2) { queue.swap(i1, i2); }
	
	std::vector<T> getHeap() { return queue.heap; }
};


template<typename T, T nullValue, typename Hash=std::hash<T>, typename Compare=std::less<T>>
class IndexedPriorityQueue
{
friend class IndexedPriorityQueueProbe<T, nullValue, Compare>;
public:
	IndexedPriorityQueue():
		heap(0),
		indexes(0)
	{
	}
	
	IndexedPriorityQueue(std::vector<T> const & initVec):
		heap(initVec),
		indexes(makeIndexes())
	{
		buildHeap();
	}
	
	void update(T const & obj)
	{
		assert(indexes.find(obj) != indexes.end());
		updateAtHeapIndex(indexes[obj]);
	}
	
	bool add(T const & value)
	{
		if(indexes.find(value) != indexes.end()) {
			return false;
		}
		
		heap.push_back(value);
		size_t index = heap.size() - 1;
		indexes[value] = index;
		
		assert(indexes.size() == heap.size());
		heapifyUp(index);
		
		return true;
	}
	
	bool remove(T const & value)
	{
		auto itr = indexes.find(value);
		if(itr == indexes.end()) {
			return false;
		}
		
		size_t index = itr->second;
		if(index == heap.size() - 1) {
			heap.pop_back();
			indexes.erase(itr);
		}
		else {
			T obj = heap.back();
			heap.pop_back();
			heap[index] = obj;
			indexes.erase(itr);
			indexes[obj] = index;
			updateAtHeapIndex(index);
		}
		assert(indexes.size() == heap.size());
		
		return true;
	}
	
	T getHead()
	{
		if(heap.size() == 0) {
			return nullValue;
		}
		
		return heap[0];
	}
	
	size_t getMaxParentIndex()
	{
		assert(heap.size() > 1);
		return (heap.size() - 2)/2;
	}
	
	bool verifyHeap()
	{
		if(heap.size() <= 1) {
			return true;
		}
		
		size_t maxParentIndex = getMaxParentIndex();
		for(size_t i = 0; i <= maxParentIndex; i++)
		{
			size_t left = 2*i + 1;
			assert(left < heap.size());
			if(compare(heap[left], heap[i])) {
                std::cerr << "A" << std::endl;
				return false;
			}
			
			size_t right = left + 1;
			if(right < heap.size()) {
				if(compare(heap[right], heap[i])) {
                    std::cerr << "B" << std::endl;
					return false;
				}
			}
		}
		assert(2 * (maxParentIndex + 1) + 1 >= heap.size());
		return true;
	}
	
	size_t size()
	{
		return heap.size();
	}
	
	bool contains(T const & value)
	{
		return indexes.find(value) != indexes.end();
	}
	
private:
	Compare compare;
	std::vector<T> heap;
	std::unordered_map<T, size_t, Hash> indexes;
	
	IndexedPriorityQueue(std::vector<T> const & initVec, bool useAsInitialHeap):
		heap(initVec),
		indexes(makeIndexes())
	{
		if(!useAsInitialHeap) {
			buildHeap();
		}
	}
	
	std::unordered_map<T, size_t, Hash> makeIndexes()
	{
		std::unordered_map<T, size_t, Hash> tmpIndexes(heap.size());
		for(size_t i = 0; i < heap.size(); i++)
		{
			assert(tmpIndexes.find(heap[i]) == tmpIndexes.end());
			indexes[heap[i]] = i;
		}
		
		return tmpIndexes;
	}
	
	void buildHeap()
	{
		if(heap.size() < 2) {
			return;
		}
		
		// Build heap
		size_t maxParentIndex = getMaxParentIndex();
		for(size_t i = maxParentIndex + 1; i > 0; i--) {
			heapifyDown(i - 1);
		}
	}
	
	bool heapifyUp(size_t i)
	{
		assert(i  < heap.size());
		
		bool moved = false;
		while(i > 0) {
			size_t parent = (i-1)/2;
			if(compare(heap[i], heap[parent])) {
				swap(i, parent);
				i = parent;
				moved = true;
			}
			else break;
		}
		return moved;
	}
	
	bool heapifyDown(size_t i)
	{
		return heapifyDown(i, false);
	}
	
	bool heapifyDown(size_t i, bool removing)
	{
		assert(i < heap.size());
		if(heap.size() < 2) {
			return false;
		}
		
		bool moved = false;
		
		size_t size = heap.size();
		size_t maxParentIndex = (size - 2) / 2;
		while(i <= maxParentIndex)
		{
			T nVal = heap[i];
			
			size_t left = 2*i + 1;
			if(left >= size) break;
			T leftVal = heap[left];
			
			size_t right = left + 1;
			T rightVal = right < size ? heap[right] : nullValue;
			
			size_t min = i;
			
			if(removing || compare(leftVal, nVal))
			{
				if(rightVal != nullValue && compare(rightVal, leftVal)) {
					min = right;
				}
				else {
					min = left;
				}
			}
			else if(rightVal != nullValue && compare(rightVal, nVal)) {
				min = right;
			}
			else {
				break;
			}
			
			moved = true;
			swap(i, min);
			i = min;
		}
		return moved;
	}
	
	void updateAtHeapIndex(size_t index)
	{
		if(!heapifyUp(index))
		{
			heapifyDown(index);
		}
	}
	
	void swap(size_t i1, size_t i2)
	{
		assert(i1 < heap.size());
		assert(i2 < heap.size());
		assert(i1 != i2);
		
		T e1 = heap[i2];
		T e2 = heap[i1];
		heap[i1] = e1;
		heap[i2] = e2;
		
		indexes[e1] = i1;
		indexes[e2] = i2;
	}
};

} // namespace varmodel

#endif // #ifndef IndexedPriorityQueue_hpp
