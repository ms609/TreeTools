#include <cstdint> /* for uint16_t */

typedef int32_t xint;
 /* Must be unsigned so we can use quicksort */
 /* If we chose signed, we'd have to impose a limit on n_children, which
  * would exclude star trees */
const xint XINT_MAX = INT32_MAX;
