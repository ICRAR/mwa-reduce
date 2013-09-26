#ifndef UVECTOR_WRAPPER_H
#define UVECTOR_WRAPPER_H

//#ifdef HAVE_UVECTOR

#include "aocommon/uvector.h"

//#else

#warning "uvector does not work on this platform: using slower std::vector instead"

#include <vector>
#include <memory>

template<typename Tp, typename Alloc = std::allocator<Tp> >
using uvector = std::vector<Tp, Alloc>;

//#endif

#endif
