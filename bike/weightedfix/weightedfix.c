#include "decode.h"
#include "weightedfix.h"

#define UNUSED(x) (void)x;

ret_t decode_weightedfix(OUT split_e_t *e,
                        IN const syndrome_t *original_s,
                        IN const ct_t *ct,
                        IN const sk_t *sk) {
    UNUSED(e);
    UNUSED(original_s);
    UNUSED(ct);
    UNUSED(sk);

    
    hello_world();
    
    DMSG("This decoder is not implemented yet.");
    BIKE_ERROR(E_DECODING_FAILURE);
}
