#ifndef AO_LANE_03_H
#define AO_LANE_03_H

#include <cstring>
#include <deque>

#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

namespace ao
{

template<typename T>
class lane
{
	public:
		typedef std::size_t size_type;
		typedef T value_type;
		
		lane() :
			_buffer(0),
			_capacity(0),
			_write_position(0),
			_free_write_space(0),
			_status(status_normal)
		{
		}
		
		explicit lane(size_t capacity) :
			_buffer(new T[capacity]),
			_capacity(capacity),
			_write_position(0),
			_free_write_space(_capacity),
			_status(status_normal)
		{
		}
		
		~lane()
		{
			delete[] _buffer;
		}
		
		void clear()
		{
			boost::mutex::scoped_lock lock(_mutex);
			_write_position = 0;
			_free_write_space = _capacity;
			_status = status_normal;
		}
		
		void write(const value_type &element)
		{
			write(&element, 1);
		}
		
		void write(const value_type *elements, size_t n)
		{
			boost::mutex::scoped_lock lock(_mutex);
			
			if(_status == status_normal)
			{
				size_t write_size = _free_write_space > n ? n : _free_write_space;
				immediate_write(elements, write_size);
				n -= write_size;
				
				while(n != 0) {
					elements += write_size;
				
					do {
						_has_changed_condition.wait(lock);
					} while(_free_write_space == 0 && _status == status_normal);
					
					write_size = _free_write_space > n ? n : _free_write_space;
					immediate_write(elements, write_size);
					n -= write_size;
				} while(n != 0);
			}
		}
		
		bool read(value_type &destination)
		{
			return read(&destination, 1) == 1;
		}
		
		size_t read(value_type *destinations, size_t n)
		{
			size_t n_left = n;
			
			boost::mutex::scoped_lock lock(_mutex);
			
			size_t free_space = free_read_space();
			size_t read_size = free_space > n ? n : free_space;
			immediate_read(destinations, read_size);
			n_left -= read_size;
			
			while(n_left != 0 && _status == status_normal)
			{
				destinations += read_size;
				
				do {
					_has_changed_condition.wait(lock);
				} while(free_read_space() == 0 && _status == status_normal);
				
				free_space = free_read_space();
				read_size = free_space > n_left ? n_left : free_space;
				immediate_read(destinations, read_size);
				n_left -= read_size;
			}
			return n - n_left;
		}
		
		void write_end()
		{
			boost::mutex::scoped_lock lock(_mutex);
			_status = status_end;
			_has_changed_condition.notify_all();
		}
		
		size_t capacity() const
		{
			return _capacity;
		}
		
		size_t size() const
		{
			boost::mutex::scoped_lock lock(_mutex);
			return _capacity - _free_write_space;
		}
		
		bool empty() const
		{
			boost::mutex::scoped_lock lock(_mutex);
			return _capacity == _free_write_space;
		}
		
		/**
		 * Change the capacity of the lane. This will erase all data in the lane.
		 */
		void resize(size_t new_capacity)
		{
			T *new_buffer = new T[new_capacity];
			delete[] _buffer;
			_buffer = new_buffer;
			_capacity = new_capacity;
			_write_position = 0;
			_free_write_space = new_capacity;
			_status = status_normal;
		}
	private:
		T *_buffer;
		
		size_t _capacity;
		
		size_t _write_position;
		
		size_t _free_write_space;
		
		enum { status_normal, status_end } _status;
		
		mutable boost::mutex _mutex;
		
		boost::condition _has_changed_condition;
		
		lane(const lane<T>& source) { }
		
		void operator=(const lane<T>& source) { }
		
		size_t read_position() const
		{
			return (_write_position + _free_write_space) % _capacity;
		}
		
		size_t free_read_space() const
		{
			return _capacity - _free_write_space;
		}
		
		void immediate_write(const value_type *elements, size_t n)
		{
			// Split the writing in two ranges if needed. The first range fits in
			// [_write_position, _capacity), the second range in [0, end). By doing
			// so, we only have to calculate the modulo in the write position once.
			if(n > 0)
			{
				size_t nPart;
				if(_write_position + n > _capacity)
				{
					nPart = _capacity - _write_position;
				} else {
					nPart = n;
				}
				for(size_t i = 0; i < nPart ; ++i, ++_write_position)
				{
					_buffer[_write_position] = elements[i];
				}
				
				_write_position = _write_position % _capacity;
				
				for(size_t i = nPart; i < n ; ++i, ++_write_position)
				{
					_buffer[_write_position] = elements[i];
				}
				
				_free_write_space -= n;
				
				// Now that there is less free write space, there is more free read
				// space and thus readers can possibly continue.
				_has_changed_condition.notify_all();
			}
		}
		
		void immediate_read(value_type *elements, size_t n)
		{
			// As with write, split in two ranges if needed. The first range fits in
			// [read_position(), _capacity), the second range in [0, end).
			if(n > 0)
			{
				size_t nPart;
				size_t position = read_position();
				if(position + n > _capacity)
				{
					nPart = _capacity - position;
				} else {
					nPart = n;
				}
				for(size_t i = 0; i < nPart ; ++i, ++position)
				{
					elements[i] = _buffer[position];
				}
				
				position = position % _capacity;
				
				for(size_t i = nPart; i < n ; ++i, ++position)
				{
					elements[i] = _buffer[position];
				}
				
				_free_write_space += n;
				
				// Now that there is more free write space, writers can possibly continue.
				_has_changed_condition.notify_all();
			}
		}
};

}

#endif // AO_LANE03_H
