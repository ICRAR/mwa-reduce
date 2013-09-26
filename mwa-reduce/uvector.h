#ifndef UVECTOR_WRAPPER_H
#define UVECTOR_WRAPPER_H

#ifdef HAVE_UVECTOR

#include "aocommon/uvector.h"

#else

#warning "uvector does not work on this platform: using slower std::vector instead"

#include <vector>

namespace ao
{
	template<typename Tp, typename Alloc = std::allocator<Tp> >
	class uvector : public std::vector<Tp> { };
}

#endif

#endif
