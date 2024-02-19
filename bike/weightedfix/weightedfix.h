#pragma once

#include "types.h"

ret_t decode_weightedfix(OUT split_e_t *e,
                        IN const syndrome_t *original_s,
                        IN const ct_t *ct,
                        IN const sk_t *sk);

void hello_world();
