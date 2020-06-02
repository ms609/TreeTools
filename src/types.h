#include <cstdint> /* for uint_fast16_t, int_fast16_t, INTPTR_MAX */

typedef int_fast16_t intx;
typedef uint_fast16_t uintx;

const intx INTX_MAX = intx(INT_FAST16_MAX);
const intx INTX_CONSERVATIVE_MAX = 0x7FFF;
const uintx UINTX_MAX = uintx(UINT_FAST16_MAX);
