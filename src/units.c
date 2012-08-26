
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "units.h"

#define _UNITS_APPLY_SCALING(U,X) (U)->characteristic_scale * (X)
#define _UNITS_APPLY_INVERSE_SCALING(U,X) (X) / (U)->characteristic_scale


