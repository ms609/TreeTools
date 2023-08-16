#ifndef TreeTools_types_
#define TreeTools_types_

#include <cstdint> /* for uint_fast32_t, int_fast32_t, INTPTR_MAX */

typedef int_fast16_t int16;
typedef int_fast32_t int32;
typedef int_fast32_t intx;
typedef uint_fast32_t uintx;

namespace TreeTools {
  const intx INTX_MAX = intx(INT_FAST32_MAX);
  const uintx UINTX_MAX = uintx(UINT_FAST32_MAX);
}

#endif
