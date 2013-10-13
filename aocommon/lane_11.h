#ifndef AO_LANE_11_H
#define AO_LANE_11_H

#include <cstring>
#include <deque>
#include <mutex>
#include <condition_variable>

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
		
		lane(const lane<T>& source) = delete;
		
		lane(lane<T>&& source)
		{
			swap(source);
		}
		
		~lane()
		{
			delete[] _buffer;
		}
		
		lane<T>& operator=(const lane<T>& source) = delete;
		
		lane<T>& operator=(lane<T>&& source)
		{
			swap(source);
			return *this;
		}
		
		void swap(lane<T>& other)
		{
			std::swap(_buffer, other._buffer);
			std::swap(_capacity, other._capacity);
			std::swap(_write_position, other._write_position);
			std::swap(_free_write_space, other._free_write_space);
			std::swap(_status, other._status);
		}
		
		void clear()
		{
			std::lock_guard<std::mutex> lock(_mutex);
			_write_position = 0;
			_free_write_space = _capacity;
			_status = status_normal;
		}
		
		void write(const value_type& element)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			
			if(_status == status_normal)
			{
				while(_free_write_space == 0)
					_writing_possible_condition.wait(lock);
				
				_buffer[_write_position] = element;
				_write_position = (_write_position+1) % _capacity;
				--_free_write_space;
				// Now that there is less free write space, there is more free read
				// space and thus readers can possibly continue.
				_reading_possible_condition.notify_all();
			}
		}
		
		void write(value_type&& element)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			
			if(_status == status_normal)
			{
				while(_free_write_space == 0)
					_writing_possible_condition.wait(lock);
				
				_buffer[_write_position] = std::move(element);
				_write_position = (_write_position+1) % _capacity;
				// Now that there is less free write space, there is more free read
				// space and thus readers can possibly continue.
				_reading_possible_condition.notify_all();
			}
		}
		
		void write(const value_type* elements, size_t n)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			
			if(_status == status_normal)
			{
				size_t write_size = _free_write_space > n ? n : _free_write_space;
				immediate_write(elements, write_size);
				n -= write_size;
				
				while(n != 0) {
					elements += write_size;
				
					do {
						_writing_possible_condition.wait(lock);
					} while(_free_write_space == 0 && _status == status_normal);
					
					write_size = _free_write_space > n ? n : _free_write_space;
					immediate_write(elements, write_size);
					n -= write_size;
				} while(n != 0);
			}
		}
		
		bool read(value_type& destination)
		{
			std::unique_lock<std::mutex> lock(_mutex);
			while(free_read_space() == 0 && _status == status_normal)
				_reading_possible_condition.wait(lock);
			if(free_read_space() == 0)
				return false;
			else
			{
				destination = _buffer[read_position()];
				++_free_write_space;
				// Now that there is more free write space, writers can possibly continue.
				_writing_possible_condition.notify_all();
				return true;
			}
		}
		
		size_t read(value_type* destinations, size_t n)
		{
			size_t n_left = n;
			
			std::unique_lock<std::mutex> lock(_mutex);
			
			size_t free_space = free_read_space();
			size_t read_size = free_space > n ? n : free_space;
			immediate_read(destinations, read_size);
			n_left -= read_size;
			
			while(n_left != 0 && _status == status_normal)
			{
				destinations += read_size;
				
				do {
					_reading_possible_condition.wait(lock);
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
			std::lock_guard<std::mutex> lock(_mutex);
			_status = status_end;
			_writing_possible_condition.notify_all();
			_reading_possible_condition.notify_all();
		}
		
		size_t capacity() const
		{
			return _capacity;
		}
		
		size_t size() const
		{
			std::lock_guard<std::mutex> lock(_mutex);
			return _capacity - _free_write_space;
		}
		
		bool empty() const
		{
			std::lock_guard<std::mutex> lock(_mutex);
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
		T* _buffer;
		
		size_t _capacity;
		
		size_t _write_position;
		
		size_t _free_write_space;
		
		enum { status_normal, status_end } _status;
		
		mutable std::mutex _mutex;
		
		std::condition_variable _writing_possible_condition, _reading_possible_condition;
		
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
				_reading_possible_condition.notify_all();
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
				_writing_possible_condition.notify_all();
			}
		}
};

template<typename T>
void swap(ao::lane<T>& first, ao::lane<T>& second)
{
	first.swap(second);
}

} // end of namespace

#endif // AO_LANE11_H
