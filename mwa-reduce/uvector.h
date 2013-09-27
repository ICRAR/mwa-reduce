#ifndef UVECTOR_WRAPPER_H
#define UVECTOR_WRAPPER_H

#ifdef HAVE_UVECTOR

#include "aocommon/uvector.h"

#else

#warning "uvector does not work on this platform: using slower std::vector instead"

#include <vector>

namespace ao {
	using std::vector;
}
#define uvector vector

#endif

#endif
