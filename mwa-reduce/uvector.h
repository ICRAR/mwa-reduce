#ifndef UVECTOR_WRAPPER_H
#define UVECTOR_WRAPPER_H

#ifdef HAVE_UVECTOR
#include "bits/uvector.h"
#else
#warning "uvector does not work on this platform: using slower std::vector instead"
#include <vector>
#define uvector std::vector
#endif

#endif
