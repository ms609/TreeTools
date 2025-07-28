#ifndef TreeTools_types_
#define TreeTools_types_

#include <cstdint> /* for uint_fast32_t, int_fast32_t, INTPTR_MAX */

using int16 = int_fast16_t;
using int32 = int_fast32_t;
using intx = int_fast32_t;
using uintx = uint_fast32_t;

namespace TreeTools {
  const intx INTX_MAX = intx(INT_FAST32_MAX);
  const uintx UINTX_MAX = uintx(UINT_FAST32_MAX);
}

#endif
