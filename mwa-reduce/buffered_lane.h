#ifndef BUFFER_LANE_H
#define BUFFER_LANE_H

#include <vector>

#include "lane.h"

template<typename Tp>
class lane_write_buffer 
{
public:
	typedef typename lane<Tp>::size_type size_type;
	typedef typename lane<Tp>::value_type value_type;
	
	lane_write_buffer(lane<Tp>& lane, size_type buffer_size) : _lane(lane), _buffer_size(buffer_size)
	{
		_buffer.reserve(buffer_size);
	}
	
	~lane_write_buffer()
	{
		flush();
	}
	
	void clear()
	{
		lane<Tp>::clear();
		_buffer.clear();
	}
		
	void write(const value_type &element)
	{
		_buffer.push_back(element);
		if(_buffer.size() == _buffer_size)
			flush();
	}
	
	void write_end()
	{
		flush();
		_lane.write_end();
	}
		
	void flush()
	{
		lane<Tp>::write(&_buffer[0], _buffer.size());
		_buffer.clear();
	}
private:
	size_type _buffer_size;
	std::vector<value_type> _buffer;
	lane<Tp>& _lane;
};

#endif
