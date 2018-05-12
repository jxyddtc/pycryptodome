/*
 *  galois.c: arithmetic in Galois Fields
 *
 * ===================================================================
 * The contents of this file are dedicated to the public domain.  To
 * the extent that dedication to the public domain is not available,
 * everyone is granted a worldwide, perpetual, royalty-free,
 * non-exclusive license to exercise all rights associated with the
 * contents of this file for any purpose whatsoever.
 * No rights are reserved.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * ===================================================================
 */

#include "common.h"

FAKE_INIT(galois)

#define ALIGNMENT 32

/**
 * A V table is a 4096 bytes table that will contain the expanded
 * GHASH key (H). It is used to speed up the GF(128) multiplication Z = X*H.
 *
 * The table contains 128 entries, one for each bit of X.
 * Each entry takes 32 bytes and can fit into the cache line of a modern
 * processor. If we assume that access to memory mapped to the same
 * cache line is somewhat constant, we can make GHASH robust again
 * cache timing attacks.
 */
typedef uint64_t t_v_tables[128][2][2];

/**
 * Create a V table. V[i] is the value H*x^i (i=0..127).
 * \param h         The 16 byte GHASH key
 * \param tables    A pointer to an allocated V table
 */
static void make_v_tables(const uint8_t h[16], t_v_tables *tables)
{
    uint64_t (*cur)[2];
    int i;

    memset(tables, 0, sizeof(t_v_tables));

    cur = &((*tables)[0][1]);

    (*cur)[0] = LOAD_U64_BIG(&h[0]);
    (*cur)[1] = LOAD_U64_BIG(&h[8]);

    for (i=1; i<128; i++) {
        uint64_t c;
        uint64_t (*next)[2];

        next = &((*tables)[i][1]);

        /** v = (v&1)*0xE1000000000000000000000000000000L ^ (v>>1) **/
        c = (*cur)[1]&1 ? 0xE100000000000000 : 0;
        (*next)[1] = (*cur)[1]>>1 | (*cur)[0]<<63;
        (*next)[0] = (*cur)[0]>>1 ^ c;

        cur = next;
    }
}

/**
 * Multiply two elements of GF(2**128) using the reducing polynomial
 * (x^128 + x^7 + x^2 + x + 1).
 *
 * \param   out         The 16 byte buffer that will receive the result
 * \param   key_tables  One factor, expanded into a V table
 * \param   x           The other factor (16 bytes)
 */
static void gcm_mult2(uint8_t out[16], const t_v_tables *key_tables, const uint8_t x[16])
{
    int i, bit_scan_128;
    uint64_t z[2];

    z[0] = z[1] = 0;
    bit_scan_128 = 0;
    for (i=0; i<16; i++) {
        unsigned xi;
        int j;

        xi = x[i];
        for (j=0; j<8; j++) {
            unsigned bit;

            bit = xi>>7 & 1; /** Constant time */
            z[0] ^= (*key_tables)[bit_scan_128][bit][0];
            z[1] ^= (*key_tables)[bit_scan_128][bit][1];

            xi <<= 1;
            bit_scan_128++;
        }
    }
    
    STORE_U64_BIG(out,   z[0]);
    STORE_U64_BIG(out+8, z[1]);
}

/**
 * Compute the GHASH of a piece of data given an arbitrary Y_0,
 * as specified in NIST SP 800 38D.
 *
 * \param y_out      The resulting GHASH (16 bytes).
 * \param block_data Pointer to the data to hash.
 * \param len        Length of the data to hash (multiple of 16).
 * \param y_in       The initial Y (Y_0, 16 bytes).
 * \param exp_key    The expanded hash key (16*256*16 bytes + alignment).
 *
 * y_out and y_int can point to the same buffer.
 */
EXPORT_SYM int ghash(
        uint8_t y_out[16],
        const uint8_t block_data[],
        size_t len,
        const uint8_t y_in[16],
        const t_v_tables *v_tables
        )
{
    unsigned i;

    if (NULL==y_out || NULL==block_data || NULL==y_in || NULL==v_tables)
        return ERR_NULL;

    if (len % 16)
        return ERR_NOT_ENOUGH_DATA;
     
    memcpy(y_out, y_in, 16);
    for (i=0; i<len; i+=16) {
        unsigned j;
        uint8_t x[16];

        for (j=0; j<16; j++) {
            x[j] = y_out[j] ^ block_data[i+j];
        }
        gcm_mult2(y_out, v_tables, x);
    }

    return 0;
}

/**
 * Expand the AES key into a Python (byte) string object.
 */ 
EXPORT_SYM int ghash_expand(const uint8_t h[16], t_v_tables **ghash_tables)
{
    if (NULL==h || NULL==ghash_tables)
        return ERR_NULL;

    *ghash_tables = align_alloc(sizeof(t_v_tables), ALIGNMENT);
    if (NULL == *ghash_tables)
        return ERR_MEMORY;
    
    make_v_tables(h, *ghash_tables);
    
    return 0;
}

EXPORT_SYM int ghash_destroy(t_v_tables *ghash_tables)
{
    free(ghash_tables);
    return 0;
}
