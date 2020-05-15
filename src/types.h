#include <cstdint> /* for uint16_t */

typedef int16_t intx;
 /* Must be unsigned so we can use quicksort */
 /* If we chose signed, we'd have to impose a limit on n_children, which
  * would exclude star trees */
const intx INTX_MAX = INT16_MAX;
