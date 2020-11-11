#ifndef TreeTools_types_
#define TreeTools_types_

#include <cstdint> /* for uint_fast16_t, int_fast16_t, INTPTR_MAX */

typedef int_fast16_t intx;
typedef int_fast16_t int16;
typedef uint_fast16_t uintx;
typedef int_fast32_t int32;

const intx INTX_MAX = intx(INT_FAST16_MAX);
const intx INTX_CONSERVATIVE_MAX = 0x7FFF;
const uintx UINTX_MAX = uintx(UINT_FAST16_MAX);

#endif
