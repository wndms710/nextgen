#include "tomcrypt.hpp"
//DAP src/hashes/helper/hash_memory.c
/**
  Hash a block of memory and store the digest.
  @param hash   The index of the hash you wish to use
  @param in     The data you wish to hash
  @param inlen  The length of the data to hash (octets)
  @param out    [out] Where to store the digest
  @param outlen [in/out] Max size and resulting size of the digest
  @return CRYPT_OK if successful
*/
int hash_memory(int hash, const unsigned char *in, unsigned long inlen, unsigned char *out, unsigned long *outlen)
{
    hash_state *md;
    int err;

    LTC_ARGCHK(in     != NULL);
    LTC_ARGCHK(out    != NULL);
    LTC_ARGCHK(outlen != NULL);

    /* Hard coded now, DAP
    if ((err = hash_is_valid(hash)) != CRYPT_OK) {
        return err;
    }
    */

    if (*outlen < hash_descriptor[hash].hashsize) {
       *outlen = hash_descriptor[hash].hashsize;
       return CRYPT_BUFFER_OVERFLOW;
    }

    md = static_cast<hash_state*>(malloc(sizeof(hash_state)));
    if (md == NULL) {
       return CRYPT_MEM;
    }

    if ((err = hash_descriptor[hash].init(md)) != CRYPT_OK) {
       goto LBL_ERR;
    }
    if ((err = hash_descriptor[hash].process(md, in, inlen)) != CRYPT_OK) {
       goto LBL_ERR;
    }
    err = hash_descriptor[hash].done(md, out);
    *outlen = hash_descriptor[hash].hashsize;
LBL_ERR:
    XFREE(md);

    return err;
}

//DAP src/hashes/sha1.c
struct ltc_hash_descriptor hash_descriptor[TAB_SIZE] = {
  { NULL, 0, 0, 0, { 0 }, 0, NULL, NULL, NULL }
};

const struct ltc_hash_descriptor sha1_desc =
{
    "sha1",
    2,
    20,
    64,

    /* OID */
   { 1, 3, 14, 3, 2, 26,  },
   6,

    &sha1_init,
    &sha1_process,
    &sha1_done
};

#define F0(x,y,z)  (z ^ (x & (y ^ z)))
#define F1(x,y,z)  (x ^ y ^ z)
#define F2(x,y,z)  ((x & y) | (z & (x | y)))
#define F3(x,y,z)  (x ^ y ^ z)

static int  sha1_compress(hash_state *md, unsigned char *buf)
{
    ulong32 a,b,c,d,e,W[80],i;

    /* copy the state into 512-bits into W[0..15] */
    for (i = 0; i < 16; i++) {
        LOAD32H(W[i], buf + (4*i));
    }

    /* copy state */
    a = md->sha1.state[0];
    b = md->sha1.state[1];
    c = md->sha1.state[2];
    d = md->sha1.state[3];
    e = md->sha1.state[4];

    /* expand it */
    for (i = 16; i < 80; i++) {
        W[i] = ROL(W[i-3] ^ W[i-8] ^ W[i-14] ^ W[i-16], 1); 
    }

    /* compress */
    /* round one */
    #define FF0(a,b,c,d,e,i) e = (ROLc(a, 5) + F0(b,c,d) + e + W[i] + 0x5a827999UL); b = ROLc(b, 30);
    #define FF1(a,b,c,d,e,i) e = (ROLc(a, 5) + F1(b,c,d) + e + W[i] + 0x6ed9eba1UL); b = ROLc(b, 30);
    #define FF2(a,b,c,d,e,i) e = (ROLc(a, 5) + F2(b,c,d) + e + W[i] + 0x8f1bbcdcUL); b = ROLc(b, 30);
    #define FF3(a,b,c,d,e,i) e = (ROLc(a, 5) + F3(b,c,d) + e + W[i] + 0xca62c1d6UL); b = ROLc(b, 30);
 
    for (i = 0; i < 20; ) {
       FF0(a,b,c,d,e,i++);
       FF0(e,a,b,c,d,i++);
       FF0(d,e,a,b,c,i++);
       FF0(c,d,e,a,b,i++);
       FF0(b,c,d,e,a,i++);
    }

    /* round two */
    for (; i < 40; )  { 
       FF1(a,b,c,d,e,i++);
       FF1(e,a,b,c,d,i++);
       FF1(d,e,a,b,c,i++);
       FF1(c,d,e,a,b,i++);
       FF1(b,c,d,e,a,i++);
    }

    /* round three */
    for (; i < 60; )  { 
       FF2(a,b,c,d,e,i++);
       FF2(e,a,b,c,d,i++);
       FF2(d,e,a,b,c,i++);
       FF2(c,d,e,a,b,i++);
       FF2(b,c,d,e,a,i++);
    }

    /* round four */
    for (; i < 80; )  { 
       FF3(a,b,c,d,e,i++);
       FF3(e,a,b,c,d,i++);
       FF3(d,e,a,b,c,i++);
       FF3(c,d,e,a,b,i++);
       FF3(b,c,d,e,a,i++);
    }

    #undef FF0
    #undef FF1
    #undef FF2
    #undef FF3

    /* store */
    md->sha1.state[0] = md->sha1.state[0] + a;
    md->sha1.state[1] = md->sha1.state[1] + b;
    md->sha1.state[2] = md->sha1.state[2] + c;
    md->sha1.state[3] = md->sha1.state[3] + d;
    md->sha1.state[4] = md->sha1.state[4] + e;

    return CRYPT_OK;
}

/**
   Initialize the hash state
   @param md   The hash state you wish to initialize
   @return CRYPT_OK if successful
*/
int sha1_init(hash_state * md)
{
   LTC_ARGCHK(md != NULL);
   md->sha1.state[0] = 0x67452301UL;
   md->sha1.state[1] = 0xefcdab89UL;
   md->sha1.state[2] = 0x98badcfeUL;
   md->sha1.state[3] = 0x10325476UL;
   md->sha1.state[4] = 0xc3d2e1f0UL;
   md->sha1.curlen = 0;
   md->sha1.length = 0;
   return CRYPT_OK;
}

/**
   Process a block of memory though the hash
   @param md     The hash state
   @param in     The data to hash
   @param inlen  The length of the data (octets)
   @return CRYPT_OK if successful
*/
HASH_PROCESS(sha1_process, sha1_compress, sha1, 64)

/**
   Terminate the hash to get the digest
   @param md  The hash state
   @param out [out] The destination of the hash (20 bytes)
   @return CRYPT_OK if successful
*/
int sha1_done(hash_state * md, unsigned char *out)
{
    int i;

    LTC_ARGCHK(md  != NULL);
    LTC_ARGCHK(out != NULL);

    if (md->sha1.curlen >= sizeof(md->sha1.buf)) {
       return CRYPT_INVALID_ARG;
    }

    /* increase the length of the message */
    md->sha1.length += md->sha1.curlen * 8;

    /* append the '1' bit */
    md->sha1.buf[md->sha1.curlen++] = (unsigned char)0x80;

    /* if the length is currently above 56 bytes we append zeros
     * then compress.  Then we can fall back to padding zeros and length
     * encoding like normal.
     */
    if (md->sha1.curlen > 56) {
        while (md->sha1.curlen < 64) {
            md->sha1.buf[md->sha1.curlen++] = (unsigned char)0;
        }
        sha1_compress(md, md->sha1.buf);
        md->sha1.curlen = 0;
    }

    /* pad upto 56 bytes of zeroes */
    while (md->sha1.curlen < 56) {
        md->sha1.buf[md->sha1.curlen++] = (unsigned char)0;
    }

    /* store length */
    STORE64H(md->sha1.length, md->sha1.buf+56);
    sha1_compress(md, md->sha1.buf);

    /* copy output */
    for (i = 0; i < 5; i++) {
        STORE32H(md->sha1.state[i], out+(4*i));
    }
    return CRYPT_OK;
}
//DAP src/math/ltm_desc.c

//DAP src/math/multi.c
#include <stdarg.h>

int ltc_init_multi(void **a, ...)
{
   void    **cur = a;
   int       np  = 0;
   va_list   args;

   va_start(args, a);
   while (cur != NULL) {
       if (mp_init(cur) != CRYPT_OK) {
          /* failed */
          va_list clean_list;

          va_start(clean_list, a);
          cur = a;
          while (np--) {
              mp_clear(*cur);
              cur = va_arg(clean_list, void**);
          }
          va_end(clean_list);
          return CRYPT_MEM;
       }
       ++np;
       cur = va_arg(args, void**);
   }
   va_end(args);
   return CRYPT_OK;   
}

void ltc_deinit_multi(void *a, ...)
{
   void     *cur = a;
   va_list   args;

   va_start(args, a);
   while (cur != NULL) {
       mp_clear(cur);
       cur = va_arg(args, void *);
   }
   va_end(args);
}

//DAP src/misc/base64/base64_decode.c
static const unsigned char map[256] = {
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255,  62, 255, 255, 255,  63,
 52,  53,  54,  55,  56,  57,  58,  59,  60,  61, 255, 255,
255, 254, 255, 255, 255,   0,   1,   2,   3,   4,   5,   6,
  7,   8,   9,  10,  11,  12,  13,  14,  15,  16,  17,  18,
 19,  20,  21,  22,  23,  24,  25, 255, 255, 255, 255, 255,
255,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,
 37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,
 49,  50,  51, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
255, 255, 255, 255 };

/**
   base64 decode a block of memory
   @param in       The base64 data to decode
   @param inlen    The length of the base64 data
   @param out      [out] The destination of the binary decoded data
   @param outlen   [in/out] The max size and resulting size of the decoded data
   @return CRYPT_OK if successful
*/
int base64_decode(const unsigned char *in,  unsigned long inlen, 
                        unsigned char *out, unsigned long *outlen)
{
   unsigned long t, x, y, z;
   unsigned char c;
   int           g;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);

   g = 3;
   for (x = y = z = t = 0; x < inlen; x++) {
       c = map[in[x]&0xFF];
       if (c == 255) continue;
       /* the final = symbols are read and used to trim the remaining bytes */
       if (c == 254) { 
          c = 0; 
          /* prevent g < 0 which would potentially allow an overflow later */
          if (--g < 0) {
             return CRYPT_INVALID_PACKET;
          }
       } else if (g != 3) {
          /* we only allow = to be at the end */
          return CRYPT_INVALID_PACKET;
       }

       t = (t<<6)|c;

       if (++y == 4) {
          if (z + g > *outlen) { 
             return CRYPT_BUFFER_OVERFLOW; 
          }
          out[z++] = (unsigned char)((t>>16)&255);
          if (g > 1) out[z++] = (unsigned char)((t>>8)&255);
          if (g > 2) out[z++] = (unsigned char)(t&255);
          y = t = 0;
       }
   }
   if (y != 0) {
       return CRYPT_INVALID_PACKET;
   }
   *outlen = z;
   return CRYPT_OK;
}

//DAP src/misc/base64/base64_encode.c
static const char *codes = 
"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/**
   base64 Encode a buffer (NUL terminated)
   @param in      The input buffer to encode
   @param inlen   The length of the input buffer
   @param out     [out] The destination of the base64 encoded data
   @param outlen  [in/out] The max size and resulting size
   @return CRYPT_OK if successful
*/
int base64_encode(const unsigned char *in,  unsigned long inlen, 
                        unsigned char *out, unsigned long *outlen)
{
   unsigned long i, len2, leven;
   unsigned char *p;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);

   /* valid output size ? */
   len2 = 4 * ((inlen + 2) / 3);
   if (*outlen < len2 + 1) {
      *outlen = len2 + 1;
      return CRYPT_BUFFER_OVERFLOW;
   }
   p = out;
   leven = 3*(inlen / 3);
   for (i = 0; i < leven; i += 3) {
       *p++ = codes[(in[0] >> 2) & 0x3F];
       *p++ = codes[(((in[0] & 3) << 4) + (in[1] >> 4)) & 0x3F];
       *p++ = codes[(((in[1] & 0xf) << 2) + (in[2] >> 6)) & 0x3F];
       *p++ = codes[in[2] & 0x3F];
       in += 3;
   }
   /* Pad it if necessary...  */
   if (i < inlen) {
       unsigned a = in[0];
       unsigned b = (i+1 < inlen) ? in[1] : 0;

       *p++ = codes[(a >> 2) & 0x3F];
       *p++ = codes[(((a & 3) << 4) + (b >> 4)) & 0x3F];
       *p++ = (i+1 < inlen) ? codes[(((b & 0xf) << 2)) & 0x3F] : '=';
       *p++ = '=';
   }

   /* append a NULL byte */
   *p = '\0';

   /* return ok */
   *outlen = p - out;
   return CRYPT_OK;
}

//DAP src/misc/crypt/crypt_argchk.c
#include <signal.h>

/**
  @file crypt_argchk.c
  Perform argument checking, Tom St Denis
*/  

void crypt_argchk(const char *v, const char *s, int d)
{
 fprintf(stderr, "LTC_ARGCHK '%s' failure on line %d of file %s\n",
         v, d, s);
 (void)raise(SIGABRT);
}

//DAP src/misc/crypt/crypt_ltc_mp_descriptor.c
ltc_math_descriptor ltc_mp;

//DAP src/misc/error_to_string.c
/**
  @file error_to_string.c
  Convert error codes to ASCII strings, Tom St Denis
*/

static const char *err_2_str[] =
{
   "CRYPT_OK",
   "CRYPT_ERROR",
   "Non-fatal 'no-operation' requested.",

   "Buffer overflow.",
   "Invalid input packet.",

   "Error reading the PRNG.",   

   "Invalid hash specified.",

   "Out of memory.",
   "A private PK key is required.",

   "Invalid argument provided.",

   "Invalid PK type.",
   "Invalid sized parameter.",

   "Invalid Padding."

};

/**
   Convert an LTC error code to ASCII
   @param err    The error code
   @return A pointer to the ASCII NUL terminated string for the error or "Invalid error code." if the err code was not valid.
*/
const char *error_to_string(int err)
{
   if (err < 0 || err >= (int)(sizeof(err_2_str)/sizeof(err_2_str[0]))) {
      return "Invalid error code.";
   } else {
      return err_2_str[err];
   }   
}

//DAP src/misc/zeromem.c
/**
   Zero a block of memory
   @param out    The destination of the area to zero
   @param outlen The length of the area to zero (octets)
*/
void zeromem(void *out, size_t outlen)
{
  unsigned char *mem = static_cast<unsigned char*>(out);
   LTC_ARGCHKVD(out != NULL);
   while (outlen-- > 0) {
      *mem++ = 0;
   }
}

//DAP src/pk/asn1/der/bit/der_decode_bit_string.c
/**
  Store a BIT STRING
  @param in      The DER encoded BIT STRING
  @param inlen   The size of the DER BIT STRING
  @param out     [out] The array of bits stored (one per char)
  @param outlen  [in/out] The number of bits stored
  @return CRYPT_OK if successful
*/
int der_decode_bit_string(const unsigned char *in,  unsigned long inlen,
                                unsigned char *out, unsigned long *outlen)
{
   unsigned long dlen, blen, x, y;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);

   /* packet must be at least 4 bytes */
   if (inlen < 4) {
       return CRYPT_INVALID_ARG;
   }

   /* check for 0x03 */
   if ((in[0]&0x1F) != 0x03) {
      return CRYPT_INVALID_PACKET;
   }

    /* offset in the data */
    x = 1;

   /* get the length of the data */
   if (in[x] & 0x80) {
      /* long format get number of length bytes */
      y = in[x++] & 0x7F;

      /* invalid if 0 or > 2 */
      if (y == 0 || y > 2) {
         return CRYPT_INVALID_PACKET;
      }

      /* read the data len */
      dlen = 0;
      while (y--) {
         dlen = (dlen << 8) | (unsigned long)in[x++];
      }
   } else {
      /* short format */
      dlen = in[x++] & 0x7F;
   }
  
   /* is the data len too long or too short? */
   if ((dlen == 0) || (dlen + x > inlen)) {
       return CRYPT_INVALID_PACKET;
   }

   /* get padding count */
   blen = ((dlen - 1) << 3) - (in[x++] & 7);

   /* too many bits? */
   if (blen > *outlen) {
      *outlen = blen;
      return CRYPT_BUFFER_OVERFLOW;
   }

   /* decode/store the bits */
   for (y = 0; y < blen; y++) {
       out[y] = (in[x] & (1 << (7 - (y & 7)))) ? 1 : 0;
       if ((y & 7) == 7) {
          ++x;
       }
   }

   /* we done */
   *outlen = blen;
   return CRYPT_OK;
}

//DAP src/pk/asn1/der/bit/der_length_bit_string.c
/**
  Gets length of DER encoding of BIT STRING 
  @param nbits  The number of bits in the string to encode
  @param outlen [out] The length of the DER encoding for the given string
  @return CRYPT_OK if successful
*/
int der_length_bit_string(unsigned long nbits, unsigned long *outlen)
{
   unsigned long nbytes;
   LTC_ARGCHK(outlen != NULL);

   /* get the number of the bytes */
   nbytes = (nbits >> 3) + ((nbits & 7) ? 1 : 0) + 1;
 
   if (nbytes < 128) {
      /* 03 LL PP DD DD DD ... */
      *outlen = 2 + nbytes;
   } else if (nbytes < 256) {
      /* 03 81 LL PP DD DD DD ... */
      *outlen = 3 + nbytes;
   } else if (nbytes < 65536) {
      /* 03 82 LL LL PP DD DD DD ... */
      *outlen = 4 + nbytes;
   } else {
      return CRYPT_INVALID_ARG;
   }

   return CRYPT_OK;
}

//DAP src/pk/asn1/der/boolean/der_decode_boolean.c
/**
  Read a BOOLEAN
  @param in      The destination for the DER encoded BOOLEAN
  @param inlen   The size of the DER BOOLEAN
  @param out     [out]  The boolean to decode
  @return CRYPT_OK if successful
*/
int der_decode_boolean(const unsigned char *in, unsigned long inlen,
                                       int *out)
{
   LTC_ARGCHK(in  != NULL);
   LTC_ARGCHK(out != NULL);
   
   if (inlen != 3 || in[0] != 0x01 || in[1] != 0x01 || (in[2] != 0x00 && in[2] != 0xFF)) {
      return CRYPT_INVALID_ARG;
   }
   
   *out = (in[2]==0xFF) ? 1 : 0;
   
   return CRYPT_OK;
}

//DAP src/pk/asn1/der/boolean/der_length_boolean.c
/**
  Gets length of DER encoding of a BOOLEAN 
  @param outlen [out] The length of the DER encoding
  @return CRYPT_OK if successful
*/
int der_length_boolean(unsigned long *outlen)
{
   LTC_ARGCHK(outlen != NULL);
   *outlen = 3;
   return CRYPT_OK;
}

//DAP src/pk/asn1/der/choice/der_decode_choice.c
/**
   Decode a CHOICE
   @param in       The DER encoded input
   @param inlen    [in/out] The size of the input and resulting size of read type
   @param list     The list of items to decode
   @param outlen   The number of items in the list
   @return CRYPT_OK on success
*/
int der_decode_choice(const unsigned char *in,   unsigned long *inlen,
                            ltc_asn1_list *list, unsigned long  outlen)
{
   unsigned long size, x, z;
   void          *data;

   LTC_ARGCHK(in    != NULL);
   LTC_ARGCHK(inlen != NULL);
   LTC_ARGCHK(list  != NULL);

   /* get blk size */
   if (*inlen < 2) {
      return CRYPT_INVALID_PACKET;
   }

   /* set all of the "used" flags to zero */
   for (x = 0; x < outlen; x++) {
       list[x].used = 0;
   }

   /* now scan until we have a winner */
   for (x = 0; x < outlen; x++) {
       size = list[x].size;
       data = list[x].data;

       switch (list[x].type) {
           case LTC_ASN1_INTEGER:
               if (der_decode_integer(in, *inlen, data) == CRYPT_OK) {
                  if (der_length_integer(data, &z) == CRYPT_OK) {
                      list[x].used = 1;
                      *inlen       = z;
                      return CRYPT_OK;
                  }
               }
               break;

           case LTC_ASN1_SHORT_INTEGER:
	     if (der_decode_short_integer(in, *inlen, static_cast<long unsigned int*>(data)) == CRYPT_OK) {
                  if (der_length_short_integer(size, &z) == CRYPT_OK) {
                      list[x].used = 1;
                      *inlen       = z;
                      return CRYPT_OK;
                  }
               }
               break;

           case LTC_ASN1_BIT_STRING:
	     if (der_decode_bit_string(in, *inlen, static_cast<unsigned char*>(data), &size) == CRYPT_OK) {
                  if (der_length_bit_string(size, &z) == CRYPT_OK) {
                     list[x].used = 1;
                     list[x].size = size;
                     *inlen       = z;
                     return CRYPT_OK;
                  }
               }
               break;

           case LTC_ASN1_OCTET_STRING:
	     if (der_decode_octet_string(in, *inlen, static_cast<unsigned char*>(data), &size) == CRYPT_OK) {
                  if (der_length_octet_string(size, &z) == CRYPT_OK) {
                     list[x].used = 1;
                     list[x].size = size;
                     *inlen       = z;
                     return CRYPT_OK;
                  }
               }
               break;

           case LTC_ASN1_NULL:
               if (*inlen == 2 && in[x] == 0x05 && in[x+1] == 0x00) {
                  *inlen = 2;
                  list[x].used   = 1;
                  return CRYPT_OK;
               }
               break;
                  
           case LTC_ASN1_OBJECT_IDENTIFIER:
	     if (der_decode_object_identifier(in, *inlen, static_cast<long unsigned int*>(data), &size) == CRYPT_OK) {
	       if (der_length_object_identifier(static_cast<long unsigned int*>(data), size, &z) == CRYPT_OK) {
                     list[x].used = 1;
                     list[x].size = size;
                     *inlen       = z;
                     return CRYPT_OK;
                  }
               }
               break;

           case LTC_ASN1_IA5_STRING:
	     if (der_decode_ia5_string(in, *inlen, static_cast<unsigned char*>(data), &size) == CRYPT_OK) {
	       if (der_length_ia5_string(static_cast<const unsigned char*>(data), size, &z) == CRYPT_OK) {
                     list[x].used = 1;
                     list[x].size = size;
                     *inlen       = z;
                     return CRYPT_OK;
                  }
               }
               break;


           case LTC_ASN1_PRINTABLE_STRING:
	     if (der_decode_printable_string(in, *inlen, static_cast<unsigned char*>(data), &size) == CRYPT_OK) {
	       if (der_length_printable_string(static_cast<const unsigned char*>(data), size, &z) == CRYPT_OK) {
                     list[x].used = 1;
                     list[x].size = size;
                     *inlen       = z;
                     return CRYPT_OK;
                  }
               }
               break;

           case LTC_ASN1_UTF8_STRING:
	     if (der_decode_utf8_string(in, *inlen, static_cast<wchar_t*>(data), &size) == CRYPT_OK) {
	       if (der_length_utf8_string(static_cast<const wchar_t*>(data), size, &z) == CRYPT_OK) {
                     list[x].used = 1;
                     list[x].size = size;
                     *inlen       = z;
                     return CRYPT_OK;
                  }
               }
               break;

           case LTC_ASN1_UTCTIME:
               z = *inlen;
               if (der_decode_utctime(in, &z, static_cast<ltc_utctime*>(data)) == CRYPT_OK) {
                  list[x].used = 1;
                  *inlen       = z;
                  return CRYPT_OK;
               }
               break;

           case LTC_ASN1_SET:
           case LTC_ASN1_SETOF:
           case LTC_ASN1_SEQUENCE:
	     if (der_decode_sequence(in, *inlen, static_cast<ltc_asn1_list*>(data), size) == CRYPT_OK) {
	       if (der_length_sequence(static_cast<ltc_asn1_list*>(data), size, &z) == CRYPT_OK) {
                     list[x].used = 1;
                     *inlen       = z;
                     return CRYPT_OK;
                  }
               }
               break;

           default:
               return CRYPT_INVALID_ARG;
       }
   }

   return CRYPT_INVALID_PACKET;
}

//DAP src/pk/asn1/der/ia5/der_decode_ia5_string.c
/**
  Store a IA5 STRING
  @param in      The DER encoded IA5 STRING
  @param inlen   The size of the DER IA5 STRING
  @param out     [out] The array of octets stored (one per char)
  @param outlen  [in/out] The number of octets stored
  @return CRYPT_OK if successful
*/
int der_decode_ia5_string(const unsigned char *in, unsigned long inlen,
                                unsigned char *out, unsigned long *outlen)
{
   unsigned long x, y, len;
   int           t;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);

   /* must have header at least */
   if (inlen < 2) {
      return CRYPT_INVALID_PACKET;
   }

   /* check for 0x16 */
   if ((in[0] & 0x1F) != 0x16) {
      return CRYPT_INVALID_PACKET;
   }
   x = 1;

   /* decode the length */
   if (in[x] & 0x80) {
      /* valid # of bytes in length are 1,2,3 */
      y = in[x] & 0x7F;
      if ((y == 0) || (y > 3) || ((x + y) > inlen)) {
         return CRYPT_INVALID_PACKET;
      }

      /* read the length in */
      len = 0;
      ++x;
      while (y--) {
         len = (len << 8) | in[x++];
      }
   } else {
      len = in[x++] & 0x7F;
   }

   /* is it too long? */
   if (len > *outlen) {
      *outlen = len;
      return CRYPT_BUFFER_OVERFLOW;
   }

   if (len + x > inlen) {
      return CRYPT_INVALID_PACKET;
   }

   /* read the data */
   for (y = 0; y < len; y++) {
       t = der_ia5_value_decode(in[x++]);
       if (t == -1) {
           return CRYPT_INVALID_ARG;
       }
       out[y] = t;
   }

   *outlen = y;

   return CRYPT_OK;
}
 
//DAP src/pk/asn1/der/ia5/der_length_ia5_string.c
static const struct {
   int code, value;
} ia5_table[] = {
{ '\0', 0 },
{ '\a', 7 }, 
{ '\b', 8 }, 
{ '\t', 9 }, 
{ '\n', 10 }, 
{ '\f', 12 }, 
{ '\r', 13 }, 
{ ' ', 32 }, 
{ '!', 33 }, 
{ '"', 34 }, 
{ '#', 35 }, 
{ '$', 36 }, 
{ '%', 37 }, 
{ '&', 38 }, 
{ '\'', 39 }, 
{ '(', 40 }, 
{ ')', 41 }, 
{ '*', 42 }, 
{ '+', 43 }, 
{ ',', 44 }, 
{ '-', 45 }, 
{ '.', 46 }, 
{ '/', 47 }, 
{ '0', 48 }, 
{ '1', 49 }, 
{ '2', 50 }, 
{ '3', 51 }, 
{ '4', 52 }, 
{ '5', 53 }, 
{ '6', 54 }, 
{ '7', 55 }, 
{ '8', 56 }, 
{ '9', 57 }, 
{ ':', 58 }, 
{ ';', 59 }, 
{ '<', 60 }, 
{ '=', 61 }, 
{ '>', 62 }, 
{ '?', 63 }, 
{ '@', 64 }, 
{ 'A', 65 }, 
{ 'B', 66 }, 
{ 'C', 67 }, 
{ 'D', 68 }, 
{ 'E', 69 }, 
{ 'F', 70 }, 
{ 'G', 71 }, 
{ 'H', 72 }, 
{ 'I', 73 }, 
{ 'J', 74 }, 
{ 'K', 75 }, 
{ 'L', 76 }, 
{ 'M', 77 }, 
{ 'N', 78 }, 
{ 'O', 79 }, 
{ 'P', 80 }, 
{ 'Q', 81 }, 
{ 'R', 82 }, 
{ 'S', 83 }, 
{ 'T', 84 }, 
{ 'U', 85 }, 
{ 'V', 86 }, 
{ 'W', 87 }, 
{ 'X', 88 }, 
{ 'Y', 89 }, 
{ 'Z', 90 }, 
{ '[', 91 }, 
{ '\\', 92 }, 
{ ']', 93 }, 
{ '^', 94 }, 
{ '_', 95 }, 
{ '`', 96 }, 
{ 'a', 97 }, 
{ 'b', 98 }, 
{ 'c', 99 }, 
{ 'd', 100 }, 
{ 'e', 101 }, 
{ 'f', 102 }, 
{ 'g', 103 }, 
{ 'h', 104 }, 
{ 'i', 105 }, 
{ 'j', 106 }, 
{ 'k', 107 }, 
{ 'l', 108 }, 
{ 'm', 109 }, 
{ 'n', 110 }, 
{ 'o', 111 }, 
{ 'p', 112 }, 
{ 'q', 113 }, 
{ 'r', 114 }, 
{ 's', 115 }, 
{ 't', 116 }, 
{ 'u', 117 }, 
{ 'v', 118 }, 
{ 'w', 119 }, 
{ 'x', 120 }, 
{ 'y', 121 }, 
{ 'z', 122 }, 
{ '{', 123 }, 
{ '|', 124 }, 
{ '}', 125 }, 
{ '~', 126 }
};

int der_ia5_char_encode(int c)
{
   int x;
   for (x = 0; x < (int)(sizeof(ia5_table)/sizeof(ia5_table[0])); x++) {
       if (ia5_table[x].code == c) {
          return ia5_table[x].value;
       }
   }
   return -1;
}

int der_ia5_value_decode(int v)
{
   int x;
   for (x = 0; x < (int)(sizeof(ia5_table)/sizeof(ia5_table[0])); x++) {
       if (ia5_table[x].value == v) {
          return ia5_table[x].code;
       }
   }
   return -1;
}
   
/**
  Gets length of DER encoding of IA5 STRING 
  @param octets   The values you want to encode 
  @param noctets  The number of octets in the string to encode
  @param outlen   [out] The length of the DER encoding for the given string
  @return CRYPT_OK if successful
*/
int der_length_ia5_string(const unsigned char *octets, unsigned long noctets, unsigned long *outlen)
{
   unsigned long x;

   LTC_ARGCHK(outlen != NULL);
   LTC_ARGCHK(octets != NULL);

   /* scan string for validity */
   for (x = 0; x < noctets; x++) {
       if (der_ia5_char_encode(octets[x]) == -1) {
          return CRYPT_INVALID_ARG;
       }
   }

   if (noctets < 128) {
      /* 16 LL DD DD DD ... */
      *outlen = 2 + noctets;
   } else if (noctets < 256) {
      /* 16 81 LL DD DD DD ... */
      *outlen = 3 + noctets;
   } else if (noctets < 65536UL) {
      /* 16 82 LL LL DD DD DD ... */
      *outlen = 4 + noctets;
   } else if (noctets < 16777216UL) {
      /* 16 83 LL LL LL DD DD DD ... */
      *outlen = 5 + noctets;
   } else {
      return CRYPT_INVALID_ARG;
   }

   return CRYPT_OK;
}

//DAP src/pk/asn1/der/integer/der_decode_integer.c
/**
  Read a mp_int integer
  @param in       The DER encoded data
  @param inlen    Size of DER encoded data
  @param num      The first mp_int to decode
  @return CRYPT_OK if successful
*/
int der_decode_integer(const unsigned char *in, unsigned long inlen, void *num)
{
   unsigned long x, y, z;
   int           err;

   LTC_ARGCHK(num    != NULL);
   LTC_ARGCHK(in     != NULL);

   /* min DER INTEGER is 0x02 01 00 == 0 */
   if (inlen < (1 + 1 + 1)) {
      return CRYPT_INVALID_PACKET;
   }

   /* ok expect 0x02 when we AND with 0001 1111 [1F] */
   x = 0;
   if ((in[x++] & 0x1F) != 0x02) {
      return CRYPT_INVALID_PACKET;
   }

   /* now decode the len stuff */
   z = in[x++];

   if ((z & 0x80) == 0x00) {
      /* short form */

      /* will it overflow? */
      if (x + z > inlen) {
         return CRYPT_INVALID_PACKET;
      }
     
      /* no so read it */
      if ((err = mp_read_unsigned_bin(num, (unsigned char *)in + x, z)) != CRYPT_OK) {
         return err;
      }
   } else {
      /* long form */
      z &= 0x7F;
      
      /* will number of length bytes overflow? (or > 4) */
      if (((x + z) > inlen) || (z > 4) || (z == 0)) {
         return CRYPT_INVALID_PACKET;
      }

      /* now read it in */
      y = 0;
      while (z--) {
         y = ((unsigned long)(in[x++])) | (y << 8);
      }

      /* now will reading y bytes overrun? */
      if ((x + y) > inlen) {
         return CRYPT_INVALID_PACKET;
      }

      /* no so read it */
      if ((err = mp_read_unsigned_bin(num, (unsigned char *)in + x, y)) != CRYPT_OK) {
         return err;
      }
   }

   /* see if it's negative */
   if (in[x] & 0x80) {
      void *tmp;
      if (mp_init(&tmp) != CRYPT_OK) {
         return CRYPT_MEM;
      }

      if (mp_2expt(tmp, mp_count_bits(num)) != CRYPT_OK || mp_sub(num, tmp, num) != CRYPT_OK) {
         mp_clear(tmp);
         return CRYPT_MEM;
      }
      mp_clear(tmp);
   } 

   return CRYPT_OK;

}

//DAP src/pk/asn1/der/integer/der_length_integer.c
/**
  Gets length of DER encoding of num 
  @param num    The int to get the size of 
  @param outlen [out] The length of the DER encoding for the given integer
  @return CRYPT_OK if successful
*/
int der_length_integer(void *num, unsigned long *outlen)
{
   unsigned long z, len;
   int           leading_zero;

   LTC_ARGCHK(num     != NULL);
   LTC_ARGCHK(outlen  != NULL);

   if (mp_cmp_d(num, 0) != LTC_MP_LT) {
      /* positive */

      /* we only need a leading zero if the msb of the first byte is one */
      if ((mp_count_bits(num) & 7) == 0 || mp_iszero(num) == LTC_MP_YES) {
         leading_zero = 1;
      } else {
         leading_zero = 0;
      }

      /* size for bignum */
      z = len = leading_zero + mp_unsigned_bin_size(num);
   } else {
      /* it's negative */
      /* find power of 2 that is a multiple of eight and greater than count bits */
      leading_zero = 0;
      z = mp_count_bits(num);
      z = z + (8 - (z & 7));
      if (((mp_cnt_lsb(num)+1)==mp_count_bits(num)) && ((mp_count_bits(num)&7)==0)) --z;
      len = z = z >> 3;
   }

   /* now we need a length */
   if (z < 128) {
      /* short form */
      ++len;
   } else {
      /* long form (relies on z != 0), assumes length bytes < 128 */
      ++len;

      while (z) {
         ++len;
         z >>= 8;
      }
   }

   /* we need a 0x02 to indicate it's INTEGER */
   ++len;

   /* return length */
   *outlen = len; 
   return CRYPT_OK;
}

//DAP src/pk/asn1/der/object_identifier/der_decode_object_identifier.c
/**
  Decode OID data and store the array of integers in words
  @param in      The OID DER encoded data
  @param inlen   The length of the OID data
  @param words   [out] The destination of the OID words
  @param outlen  [in/out] The number of OID words
  @return CRYPT_OK if successful
*/
int der_decode_object_identifier(const unsigned char *in,    unsigned long  inlen,
                                       unsigned long *words, unsigned long *outlen)
{
   unsigned long x, y, t, len;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(words  != NULL);
   LTC_ARGCHK(outlen != NULL);

   /* header is at least 3 bytes */
   if (inlen < 3) {
      return CRYPT_INVALID_PACKET;
   }

   /* must be room for at least two words */
   if (*outlen < 2) {
      return CRYPT_BUFFER_OVERFLOW;
   }

   /* decode the packet header */
   x = 0;
   if ((in[x++] & 0x1F) != 0x06) {
      return CRYPT_INVALID_PACKET;
   }
   
   /* get the length */
   if (in[x] < 128) {
      len = in[x++]; 
   } else {
       if (in[x] < 0x81 || in[x] > 0x82) {
          return CRYPT_INVALID_PACKET;
       }
       y   = in[x++] & 0x7F;
       len = 0;
       while (y--) {
          len = (len << 8) | (unsigned long)in[x++];
       }
   }

   if (len < 1 || (len + x) > inlen) {
      return CRYPT_INVALID_PACKET;
   }

   /* decode words */
   y = 0;
   t = 0;
   while (len--) {
       t = (t << 7) | (in[x] & 0x7F);
       if (!(in[x++] & 0x80)) {
           /* store t */
           if (y >= *outlen) {
              return CRYPT_BUFFER_OVERFLOW;
           }
      if (y == 0) {
         words[0] = t / 40;
         words[1] = t % 40;
         y = 2;
      } else {
              words[y++] = t;
      }
           t          = 0;
       }
   }
       
   *outlen = y;
   return CRYPT_OK;
}

//DAP src/pk/asn1/der/object_identifier/der_length_object_identifier.c
unsigned long der_object_identifier_bits(unsigned long x)
{
   unsigned long c;
   x &= 0xFFFFFFFF;
   c  = 0;
   while (x) {
     ++c;
     x >>= 1;
   }
   return c;
}


/**
  Gets length of DER encoding of Object Identifier
  @param nwords   The number of OID words 
  @param words    The actual OID words to get the size of
  @param outlen   [out] The length of the DER encoding for the given string
  @return CRYPT_OK if successful
*/
int der_length_object_identifier(unsigned long *words, unsigned long nwords, unsigned long *outlen)
{
   unsigned long y, z, t, wordbuf;   

   LTC_ARGCHK(words  != NULL);
   LTC_ARGCHK(outlen != NULL);


   /* must be >= 2 words */
   if (nwords < 2) {
      return CRYPT_INVALID_ARG;
   }

   /* word1 = 0,1,2,3 and word2 0..39 */
   if (words[0] > 3 || (words[0] < 2 && words[1] > 39)) {
      return CRYPT_INVALID_ARG;
   }

   /* leading word is the first two */
   z = 0;
   wordbuf = words[0] * 40 + words[1];
   for (y = 1; y < nwords; y++) {
       t = der_object_identifier_bits(wordbuf);
       z += t/7 + ((t%7) ? 1 : 0) + (wordbuf == 0 ? 1 : 0);
       if (y < nwords - 1) {
          /* grab next word */
          wordbuf = words[y+1];
       }
   }

   /* now depending on the length our length encoding changes */
   if (z < 128) {
      z += 2;
   } else if (z < 256) {
      z += 3;
   } else if (z < 65536UL) {
      z += 4;
   } else {
      return CRYPT_INVALID_ARG;
   }

   *outlen = z;
   return CRYPT_OK;
}

//DAP src/pk/asn1/der/octet/der_decode_octet_string.c
/**
  Store a OCTET STRING
  @param in      The DER encoded OCTET STRING
  @param inlen   The size of the DER OCTET STRING
  @param out     [out] The array of octets stored (one per char)
  @param outlen  [in/out] The number of octets stored
  @return CRYPT_OK if successful
*/
int der_decode_octet_string(const unsigned char *in, unsigned long inlen,
                                  unsigned char *out, unsigned long *outlen)
{
   unsigned long x, y, len;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);

   /* must have header at least */
   if (inlen < 2) {
      return CRYPT_INVALID_PACKET;
   }

   /* check for 0x04 */
   if ((in[0] & 0x1F) != 0x04) {
      return CRYPT_INVALID_PACKET;
   }
   x = 1;

   /* decode the length */
   if (in[x] & 0x80) {
      /* valid # of bytes in length are 1,2,3 */
      y = in[x] & 0x7F;
      if ((y == 0) || (y > 3) || ((x + y) > inlen)) {
         return CRYPT_INVALID_PACKET;
      }

      /* read the length in */
      len = 0;
      ++x;
      while (y--) {
         len = (len << 8) | in[x++];
      }
   } else {
      len = in[x++] & 0x7F;
   }

   /* is it too long? */
   if (len > *outlen) {
      *outlen = len;
      return CRYPT_BUFFER_OVERFLOW;
   }

   if (len + x > inlen) {
      return CRYPT_INVALID_PACKET;
   }

   /* read the data */
   for (y = 0; y < len; y++) {
       out[y] = in[x++];
   }

   *outlen = y;

   return CRYPT_OK;
}

//DAP src/pk/asn1/der/octet/der_length_octet_string.c
/**
  Gets length of DER encoding of OCTET STRING 
  @param noctets  The number of octets in the string to encode
  @param outlen   [out] The length of the DER encoding for the given string
  @return CRYPT_OK if successful
*/
int der_length_octet_string(unsigned long noctets, unsigned long *outlen)
{
   LTC_ARGCHK(outlen != NULL);

   if (noctets < 128) {
      /* 04 LL DD DD DD ... */
      *outlen = 2 + noctets;
   } else if (noctets < 256) {
      /* 04 81 LL DD DD DD ... */
      *outlen = 3 + noctets;
   } else if (noctets < 65536UL) {
      /* 04 82 LL LL DD DD DD ... */
      *outlen = 4 + noctets;
   } else if (noctets < 16777216UL) {
      /* 04 83 LL LL LL DD DD DD ... */
      *outlen = 5 + noctets;
   } else {
      return CRYPT_INVALID_ARG;
   }

   return CRYPT_OK;
}

//DAP src/pk/asn1/der/printable_string/der_decode_printable_string.c
/**
  Store a printable STRING
  @param in      The DER encoded printable STRING
  @param inlen   The size of the DER printable STRING
  @param out     [out] The array of octets stored (one per char)
  @param outlen  [in/out] The number of octets stored
  @return CRYPT_OK if successful
*/
int der_decode_printable_string(const unsigned char *in, unsigned long inlen,
                                unsigned char *out, unsigned long *outlen)
{
   unsigned long x, y, len;
   int           t;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);

   /* must have header at least */
   if (inlen < 2) {
      return CRYPT_INVALID_PACKET;
   }

   /* check for 0x13 */
   if ((in[0] & 0x1F) != 0x13) {
      return CRYPT_INVALID_PACKET;
   }
   x = 1;

   /* decode the length */
   if (in[x] & 0x80) {
      /* valid # of bytes in length are 1,2,3 */
      y = in[x] & 0x7F;
      if ((y == 0) || (y > 3) || ((x + y) > inlen)) {
         return CRYPT_INVALID_PACKET;
      }

      /* read the length in */
      len = 0;
      ++x;
      while (y--) {
         len = (len << 8) | in[x++];
      }
   } else {
      len = in[x++] & 0x7F;
   }

   /* is it too long? */
   if (len > *outlen) {
      *outlen = len;
      return CRYPT_BUFFER_OVERFLOW;
   }

   if (len + x > inlen) {
      return CRYPT_INVALID_PACKET;
   }

   /* read the data */
   for (y = 0; y < len; y++) {
       t = der_printable_value_decode(in[x++]);
       if (t == -1) {
           return CRYPT_INVALID_ARG;
       }
       out[y] = t;
   }

   *outlen = y;

   return CRYPT_OK;
}
 
//DAP src/pk/asn1/der/printable_string/der_length_printable_string.c
static const struct {
   int code, value;
} printable_table[] = {
{ ' ', 32 }, 
{ '\'', 39 }, 
{ '(', 40 }, 
{ ')', 41 }, 
{ '+', 43 }, 
{ ',', 44 }, 
{ '-', 45 }, 
{ '.', 46 }, 
{ '/', 47 }, 
{ '0', 48 }, 
{ '1', 49 }, 
{ '2', 50 }, 
{ '3', 51 }, 
{ '4', 52 }, 
{ '5', 53 }, 
{ '6', 54 }, 
{ '7', 55 }, 
{ '8', 56 }, 
{ '9', 57 }, 
{ ':', 58 }, 
{ '=', 61 }, 
{ '?', 63 }, 
{ 'A', 65 }, 
{ 'B', 66 }, 
{ 'C', 67 }, 
{ 'D', 68 }, 
{ 'E', 69 }, 
{ 'F', 70 }, 
{ 'G', 71 }, 
{ 'H', 72 }, 
{ 'I', 73 }, 
{ 'J', 74 }, 
{ 'K', 75 }, 
{ 'L', 76 }, 
{ 'M', 77 }, 
{ 'N', 78 }, 
{ 'O', 79 }, 
{ 'P', 80 }, 
{ 'Q', 81 }, 
{ 'R', 82 }, 
{ 'S', 83 }, 
{ 'T', 84 }, 
{ 'U', 85 }, 
{ 'V', 86 }, 
{ 'W', 87 }, 
{ 'X', 88 }, 
{ 'Y', 89 }, 
{ 'Z', 90 }, 
{ 'a', 97 }, 
{ 'b', 98 }, 
{ 'c', 99 }, 
{ 'd', 100 }, 
{ 'e', 101 }, 
{ 'f', 102 }, 
{ 'g', 103 }, 
{ 'h', 104 }, 
{ 'i', 105 }, 
{ 'j', 106 }, 
{ 'k', 107 }, 
{ 'l', 108 }, 
{ 'm', 109 }, 
{ 'n', 110 }, 
{ 'o', 111 }, 
{ 'p', 112 }, 
{ 'q', 113 }, 
{ 'r', 114 }, 
{ 's', 115 }, 
{ 't', 116 }, 
{ 'u', 117 }, 
{ 'v', 118 }, 
{ 'w', 119 }, 
{ 'x', 120 }, 
{ 'y', 121 }, 
{ 'z', 122 }, 
};

int der_printable_char_encode(int c)
{
   int x;
   for (x = 0; x < (int)(sizeof(printable_table)/sizeof(printable_table[0])); x++) {
       if (printable_table[x].code == c) {
          return printable_table[x].value;
       }
   }
   return -1;
}

int der_printable_value_decode(int v)
{
   int x;
   for (x = 0; x < (int)(sizeof(printable_table)/sizeof(printable_table[0])); x++) {
       if (printable_table[x].value == v) {
          return printable_table[x].code;
       }
   }
   return -1;
}
   
/**
  Gets length of DER encoding of Printable STRING 
  @param octets   The values you want to encode 
  @param noctets  The number of octets in the string to encode
  @param outlen   [out] The length of the DER encoding for the given string
  @return CRYPT_OK if successful
*/
int der_length_printable_string(const unsigned char *octets, unsigned long noctets, unsigned long *outlen)
{
   unsigned long x;

   LTC_ARGCHK(outlen != NULL);
   LTC_ARGCHK(octets != NULL);

   /* scan string for validity */
   for (x = 0; x < noctets; x++) {
       if (der_printable_char_encode(octets[x]) == -1) {
          return CRYPT_INVALID_ARG;
       }
   }

   if (noctets < 128) {
      /* 16 LL DD DD DD ... */
      *outlen = 2 + noctets;
   } else if (noctets < 256) {
      /* 16 81 LL DD DD DD ... */
      *outlen = 3 + noctets;
   } else if (noctets < 65536UL) {
      /* 16 82 LL LL DD DD DD ... */
      *outlen = 4 + noctets;
   } else if (noctets < 16777216UL) {
      /* 16 83 LL LL LL DD DD DD ... */
      *outlen = 5 + noctets;
   } else {
      return CRYPT_INVALID_ARG;
   }

   return CRYPT_OK;
}

//DAP src/pk/asn1/der/sequence/der_decode_sequence_ex.c
#include <stdarg.h>
/**
   Decode a SEQUENCE
   @param in       The DER encoded input
   @param inlen    The size of the input
   @param list     The list of items to decode
   @param outlen   The number of items in the list
   @param ordered  Search an unordeded or ordered list
   @return CRYPT_OK on success
*/
int der_decode_sequence_ex(const unsigned char *in, unsigned long  inlen,
                           ltc_asn1_list *list,     unsigned long  outlen, int ordered)
{
   int           err, type;
   unsigned long size, x, y, z, i, blksize=0;
   void          *data;

   LTC_ARGCHK(in   != NULL);
   LTC_ARGCHK(list != NULL);
   
   /* get blk size */
   if (inlen < 2) {
      return CRYPT_INVALID_PACKET;
   }

   /* sequence type? We allow 0x30 SEQUENCE and 0x31 SET since fundamentally they're the same structure */
   x = 0;
   if (in[x] != 0x30 && in[x] != 0x31) {
      return CRYPT_INVALID_PACKET;
   }
   ++x;

   if (in[x] < 128) {
      blksize = in[x++];
   } else if (in[x] & 0x80) {
      if (in[x] < 0x81 || in[x] > 0x83) {
         return CRYPT_INVALID_PACKET;
      }
      y = in[x++] & 0x7F;

      /* would reading the len bytes overrun? */
      if (x + y > inlen) {
         return CRYPT_INVALID_PACKET;
      }

      /* read len */
      blksize = 0;
      while (y--) {
          blksize = (blksize << 8) | (unsigned long)in[x++];
      }
  }

  /* would this blksize overflow? */
  if (x + blksize > inlen) {
     return CRYPT_INVALID_PACKET;
  }

   /* mark all as unused */
   for (i = 0; i < outlen; i++) {
       list[i].used = 0;
   }     

  /* ok read data */
   inlen = blksize;
   for (i = 0; i < outlen; i++) {
       z    = 0;
       type = list[i].type;
       size = list[i].size;
       data = list[i].data;
       if (!ordered && list[i].used == 1) { continue; }

       if (type == LTC_ASN1_EOL) { 
          break;
       }

       switch (type) {
           case LTC_ASN1_BOOLEAN:
               z = inlen;
               if ((err = der_decode_boolean(in + x, z, ((int *)data))) != CRYPT_OK) {
                   goto LBL_ERR;
               }
               if ((err = der_length_boolean(&z)) != CRYPT_OK) {
                   goto LBL_ERR;
                }
                break;
          
           case LTC_ASN1_INTEGER:
               z = inlen;
               if ((err = der_decode_integer(in + x, z, data)) != CRYPT_OK) {
                  if (!ordered) {  continue; }
                  goto LBL_ERR;
               }
               if ((err = der_length_integer(data, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;

           case LTC_ASN1_SHORT_INTEGER:
               z = inlen;
               if ((err = der_decode_short_integer(in + x, z, static_cast<long unsigned int*>(data))) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               if ((err = der_length_short_integer((static_cast<unsigned long*>(data))[0], &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               
               break;

           case LTC_ASN1_BIT_STRING:
               z = inlen;
               if ((err = der_decode_bit_string(in + x, z, static_cast<unsigned char*>(data), &size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               list[i].size = size;
               if ((err = der_length_bit_string(size, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;

           case LTC_ASN1_OCTET_STRING:
               z = inlen;
               if ((err = der_decode_octet_string(in + x, z, static_cast<unsigned char*>(data), &size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               list[i].size = size;
               if ((err = der_length_octet_string(size, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;

           case LTC_ASN1_NULL:
               if (inlen < 2 || in[x] != 0x05 || in[x+1] != 0x00) {
                  if (!ordered) { continue; }
                  err = CRYPT_INVALID_PACKET;
                  goto LBL_ERR;
               }
               z = 2;
               break;
                  
           case LTC_ASN1_OBJECT_IDENTIFIER:
               z = inlen;
               if ((err = der_decode_object_identifier(in + x, z, static_cast<long unsigned int*>(data), &size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               list[i].size = size;
               if ((err = der_length_object_identifier(static_cast<long unsigned int*>(data), size, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;

           case LTC_ASN1_IA5_STRING:
               z = inlen;
               if ((err = der_decode_ia5_string(in + x, z, static_cast<unsigned char*>(data), &size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               list[i].size = size;
               if ((err = der_length_ia5_string(static_cast<const unsigned char*>(data), size, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;


           case LTC_ASN1_PRINTABLE_STRING:
               z = inlen;
               if ((err = der_decode_printable_string(in + x, z, static_cast<unsigned char*>(data), &size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               list[i].size = size;
               if ((err = der_length_printable_string(static_cast<const unsigned char*>(data), size, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;

           case LTC_ASN1_UTF8_STRING:
               z = inlen;
               if ((err = der_decode_utf8_string(in + x, z, static_cast<wchar_t*>(data), &size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               list[i].size = size;
               if ((err = der_length_utf8_string(static_cast<const wchar_t*>(data), size, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;

           case LTC_ASN1_UTCTIME:
               z = inlen;
               if ((err = der_decode_utctime(in + x, &z, static_cast<ltc_utctime*>(data))) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               break;

           case LTC_ASN1_SET:
               z = inlen;
               if ((err = der_decode_set(in + x, z, static_cast<ltc_asn1_list*>(data), size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               if ((err = der_length_sequence(static_cast<ltc_asn1_list*>(data), size, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;
           
           case LTC_ASN1_SETOF:
           case LTC_ASN1_SEQUENCE:
               /* detect if we have the right type */
               if ((type == LTC_ASN1_SETOF && (in[x] & 0x3F) != 0x31) || (type == LTC_ASN1_SEQUENCE && (in[x] & 0x3F) != 0x30)) {
                  err = CRYPT_INVALID_PACKET;
                  goto LBL_ERR;
               }

               z = inlen;
               if ((err = der_decode_sequence(in + x, z, static_cast<ltc_asn1_list*>(data), size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               if ((err = der_length_sequence(static_cast<ltc_asn1_list*>(data), size, &z)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               break;


           case LTC_ASN1_CHOICE:
               z = inlen;
               if ((err = der_decode_choice(in + x, &z, static_cast<ltc_asn1_list*>(data), size)) != CRYPT_OK) {
                  if (!ordered) { continue; }
                  goto LBL_ERR;
               }
               break;

           default:
               err = CRYPT_INVALID_ARG;
               goto LBL_ERR;
       }
       x           += z;
       inlen       -= z;
       list[i].used = 1;
       if (!ordered) { 
	 /* restart the decoder */
	 i = 0; // was -1; BH?
       }          
   }
     
   for (i = 0; i < outlen; i++) {
      if (list[i].used == 0) {
          err = CRYPT_INVALID_PACKET;
          goto LBL_ERR;
      }
   }                
   err = CRYPT_OK;   

LBL_ERR:
   return err;
}  
 
//DAP src/pk/asn1/der/sequence/der_decode_sequence_multi.c
#include <stdarg.h>
/**
  Decode a SEQUENCE type using a VA list
  @param in    Input buffer
  @param inlen Length of input in octets
  @remark <...> is of the form <type, size, data> (int, unsigned long, void*)
  @return CRYPT_OK on success
*/  
int der_decode_sequence_multi(const unsigned char *in, unsigned long inlen, ...)
{
   int           err, type;
   unsigned long size, x;
   void          *data;
   va_list       args;
   ltc_asn1_list *list;

   LTC_ARGCHK(in    != NULL);

   /* get size of output that will be required */
   va_start(args, inlen);
   x = 0;
   for (;;) {
       type = va_arg(args, int);
       size = va_arg(args, unsigned long);
       data = va_arg(args, void*);

       if (type == LTC_ASN1_EOL) { 
          break;
       }

       switch (type) {
           case LTC_ASN1_BOOLEAN:
           case LTC_ASN1_INTEGER:
           case LTC_ASN1_SHORT_INTEGER:
           case LTC_ASN1_BIT_STRING:
           case LTC_ASN1_OCTET_STRING:
           case LTC_ASN1_NULL:
           case LTC_ASN1_OBJECT_IDENTIFIER:
           case LTC_ASN1_IA5_STRING:
           case LTC_ASN1_PRINTABLE_STRING:
           case LTC_ASN1_UTF8_STRING:
           case LTC_ASN1_UTCTIME:
           case LTC_ASN1_SET:
           case LTC_ASN1_SETOF:
           case LTC_ASN1_SEQUENCE:
           case LTC_ASN1_CHOICE:
                ++x; 
                break;
          
           default:
               va_end(args);
               return CRYPT_INVALID_ARG;
       }
   }
   va_end(args);

   /* allocate structure for x elements */
   if (x == 0) {
      return CRYPT_NOP;
   }

   list = static_cast<ltc_asn1_list*>(calloc(sizeof(*list), x));
   if (list == NULL) {
      return CRYPT_MEM;
   }

   /* fill in the structure */
   va_start(args, inlen);
   x = 0;
   for (;;) {
       type = va_arg(args, int);
       size = va_arg(args, unsigned long);
       data = va_arg(args, void*);

       if (type == LTC_ASN1_EOL) { 
          break;
       }

       switch (type) {
           case LTC_ASN1_BOOLEAN:
           case LTC_ASN1_INTEGER:
           case LTC_ASN1_SHORT_INTEGER:
           case LTC_ASN1_BIT_STRING:
           case LTC_ASN1_OCTET_STRING:
           case LTC_ASN1_NULL:
           case LTC_ASN1_OBJECT_IDENTIFIER:
           case LTC_ASN1_IA5_STRING:
           case LTC_ASN1_PRINTABLE_STRING:
           case LTC_ASN1_UTF8_STRING:
           case LTC_ASN1_UTCTIME:
           case LTC_ASN1_SEQUENCE:
           case LTC_ASN1_SET:
           case LTC_ASN1_SETOF:          
           case LTC_ASN1_CHOICE:
                list[x].type   = type;
                list[x].size   = size;
                list[x++].data = data;
                break;
         
           default:
               va_end(args);
               err = CRYPT_INVALID_ARG;
               goto LBL_ERR;
       }
   }
   va_end(args);

   err = der_decode_sequence(in, inlen, list, x);
LBL_ERR:
   XFREE(list);
   return err;
}

//DAP src/pk/asn1/der/sequence/der_length_sequence.c
/**
   Get the length of a DER sequence 
   @param list   The sequences of items in the SEQUENCE
   @param inlen  The number of items
   @param outlen [out] The length required in octets to store it 
   @return CRYPT_OK on success
*/
int der_length_sequence(ltc_asn1_list *list, unsigned long inlen,
                        unsigned long *outlen) 
{
   int           err, type;
   unsigned long size, x, y, z, i;
   void          *data;

   LTC_ARGCHK(list    != NULL);
   LTC_ARGCHK(outlen  != NULL);

   /* get size of output that will be required */
   y = 0;
   for (i = 0; i < inlen; i++) {
       type = list[i].type;
       size = list[i].size;
       data = list[i].data;

       if (type == LTC_ASN1_EOL) { 
          break;
       }

       switch (type) {
           case LTC_ASN1_BOOLEAN:
              if ((err = der_length_boolean(&x)) != CRYPT_OK) {
                 goto LBL_ERR;
              }
              y += x;
              break;
          
           case LTC_ASN1_INTEGER:
               if ((err = der_length_integer(data, &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_SHORT_INTEGER:
               if ((err = der_length_short_integer(*((unsigned long *)data), &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_BIT_STRING:
               if ((err = der_length_bit_string(size, &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_OCTET_STRING:
               if ((err = der_length_octet_string(size, &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_NULL:
               y += 2;
               break;

           case LTC_ASN1_OBJECT_IDENTIFIER:
	     if ((err = der_length_object_identifier(static_cast<long unsigned int*>(data), size, &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_IA5_STRING:
	     if ((err = der_length_ia5_string(static_cast<const unsigned char*>(data), size, &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_PRINTABLE_STRING:
	     if ((err = der_length_printable_string(static_cast<const unsigned char*>(data), size, &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_UTCTIME:
	     if ((err = der_length_utctime(static_cast<ltc_utctime*>(data), &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_UTF8_STRING:
	     if ((err = der_length_utf8_string(static_cast<const wchar_t*>(data), size, &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

           case LTC_ASN1_SET:
           case LTC_ASN1_SETOF:
           case LTC_ASN1_SEQUENCE:
	     if ((err = der_length_sequence(static_cast<ltc_asn1_list*>(data), size, &x)) != CRYPT_OK) {
                  goto LBL_ERR;
               }
               y += x;
               break;

          
           default:
               err = CRYPT_INVALID_ARG;
               goto LBL_ERR;
       }
   }

   /* calc header size */
   z = y;
   if (y < 128) {
      y += 2;
   } else if (y < 256) {
      /* 0x30 0x81 LL */
      y += 3;
   } else if (y < 65536UL) {
      /* 0x30 0x82 LL LL */
      y += 4;
   } else if (y < 16777216UL) {
      /* 0x30 0x83 LL LL LL */
      y += 5;
   } else {
      err = CRYPT_INVALID_ARG;
      goto LBL_ERR;
   }

   /* store size */
   *outlen = y;
   err     = CRYPT_OK;

LBL_ERR:
   return err;
}

//DAP src/pk/asn1/der/short_integer/der_decode_short_integer.c
/**
  Read a short integer
  @param in       The DER encoded data
  @param inlen    Size of data
  @param num      [out] The integer to decode
  @return CRYPT_OK if successful
*/
int der_decode_short_integer(const unsigned char *in, unsigned long inlen, unsigned long *num)
{
   unsigned long len, x, y;

   LTC_ARGCHK(num    != NULL);
   LTC_ARGCHK(in     != NULL);

   /* check length */
   if (inlen < 2) {
      return CRYPT_INVALID_PACKET;
   }

   /* check header */
   x = 0;
   if ((in[x++] & 0x1F) != 0x02) {
      return CRYPT_INVALID_PACKET;
   }

   /* get the packet len */
   len = in[x++];

   if (x + len > inlen) {
      return CRYPT_INVALID_PACKET;
   }

   /* read number */
   y = 0;
   while (len--) {
      y = (y<<8) | (unsigned long)in[x++];
   }
   *num = y;

   return CRYPT_OK;

}

//DAP src/pk/asn1/der/short_integer/der_length_short_integer.c
/**
  Gets length of DER encoding of num 
  @param num    The integer to get the size of 
  @param outlen [out] The length of the DER encoding for the given integer
  @return CRYPT_OK if successful
*/
int der_length_short_integer(unsigned long num, unsigned long *outlen)
{
   unsigned long z, y, len;

   LTC_ARGCHK(outlen  != NULL);

   /* force to 32 bits */
   num &= 0xFFFFFFFFUL;

   /* get the number of bytes */
   z = 0;
   y = num;
   while (y) {
     ++z;
     y >>= 8;
   }
   
   /* handle zero */
   if (z == 0) {
      z = 1;
   }

   /* we need a 0x02 to indicate it's INTEGER */
   len = 1;

   /* length byte */
   ++len;

   /* bytes in value */
   len += z;

   /* see if msb is set */
   len += (num&(1UL<<((z<<3) - 1))) ? 1 : 0;

   /* return length */
   *outlen = len; 
   
   return CRYPT_OK;
}

//DAP src/pk/asn1/der/utctime/der_decode_utctime.c
static int char_to_int(unsigned char x)
{
   switch (x)  {
      case '0': return 0;
      case '1': return 1;
      case '2': return 2;
      case '3': return 3;
      case '4': return 4;
      case '5': return 5;
      case '6': return 6;
      case '7': return 7;
      case '8': return 8;
      case '9': return 9;
   }
   return 100;
}

#define DECODE_V(y, max) \
   y  = char_to_int(buf[x])*10 + char_to_int(buf[x+1]); \
   if (y >= max) return CRYPT_INVALID_PACKET;           \
   x += 2;

/**
  Decodes a UTC time structure in DER format (reads all 6 valid encoding formats)
  @param in     Input buffer
  @param inlen  Length of input buffer in octets
  @param out    [out] Destination of UTC time structure
  @return CRYPT_OK   if successful
*/
int der_decode_utctime(const unsigned char *in, unsigned long *inlen,
                             ltc_utctime   *out)
{
   unsigned char buf[32];
   unsigned long x;
   int           y;

   LTC_ARGCHK(in    != NULL);
   LTC_ARGCHK(inlen != NULL);
   LTC_ARGCHK(out   != NULL);

   /* check header */
   if (*inlen < 2UL || (in[1] >= sizeof(buf)) || ((in[1] + 2UL) > *inlen)) {
      return CRYPT_INVALID_PACKET;
   }

   /* decode the string */
   for (x = 0; x < in[1]; x++) {
       y = der_ia5_value_decode(in[x+2]);
       if (y == -1) {
          return CRYPT_INVALID_PACKET;
       }
       buf[x] = y;
   }
   *inlen = 2 + x;


   /* possible encodings are 
YYMMDDhhmmZ
YYMMDDhhmm+hh'mm'
YYMMDDhhmm-hh'mm'
YYMMDDhhmmssZ
YYMMDDhhmmss+hh'mm'
YYMMDDhhmmss-hh'mm'

    So let's do a trivial decode upto [including] mm 
   */

    x = 0;
    DECODE_V(out->YY, 100);
    DECODE_V(out->MM, 13);
    DECODE_V(out->DD, 32);
    DECODE_V(out->hh, 24);
    DECODE_V(out->mm, 60);

    /* clear timezone and seconds info */
    out->off_dir = out->off_hh = out->off_mm = out->ss = 0;

    /* now is it Z, +, - or 0-9 */
    if (buf[x] == 'Z') {
       return CRYPT_OK;
    } else if (buf[x] == '+' || buf[x] == '-') {
       out->off_dir = (buf[x++] == '+') ? 0 : 1;
       DECODE_V(out->off_hh, 24);
       DECODE_V(out->off_mm, 60);
       return CRYPT_OK;
    }

    /* decode seconds */
    DECODE_V(out->ss, 60);

    /* now is it Z, +, - */
    if (buf[x] == 'Z') {
       return CRYPT_OK;
    } else if (buf[x] == '+' || buf[x] == '-') {
       out->off_dir = (buf[x++] == '+') ? 0 : 1;
       DECODE_V(out->off_hh, 24);
       DECODE_V(out->off_mm, 60);
       return CRYPT_OK;
    } else {
       return CRYPT_INVALID_PACKET;
    }
}

//DAP src/pk/asn1/der/utctime/der_length_utctime.c
/**
  Gets length of DER encoding of UTCTIME
  @param utctime      The UTC time structure to get the size of
  @param outlen [out] The length of the DER encoding
  @return CRYPT_OK if successful
*/
int der_length_utctime(ltc_utctime *utctime, unsigned long *outlen)
{
   LTC_ARGCHK(outlen  != NULL);
   LTC_ARGCHK(utctime != NULL);

   if (utctime->off_hh == 0 && utctime->off_mm == 0) {
      /* we encode as YYMMDDhhmmssZ */
      *outlen = 2 + 13;
   } else {
      /* we encode as YYMMDDhhmmss{+|-}hh'mm' */
      *outlen = 2 + 17;
   }

   return CRYPT_OK;
}

//DAP src/pk/asn1/der/utf8/der_decode_utf8_string.c
/**
  Store a UTF8 STRING
  @param in      The DER encoded UTF8 STRING
  @param inlen   The size of the DER UTF8 STRING
  @param out     [out] The array of utf8s stored (one per char)
  @param outlen  [in/out] The number of utf8s stored
  @return CRYPT_OK if successful
*/
int der_decode_utf8_string(const unsigned char *in,  unsigned long inlen,
                                       wchar_t *out, unsigned long *outlen)
{
   wchar_t       tmp;
   unsigned long x, y, z, len;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);

   /* must have header at least */
   if (inlen < 2) {
      return CRYPT_INVALID_PACKET;
   }

   /* check for 0x0C */
   if ((in[0] & 0x1F) != 0x0C) {
      return CRYPT_INVALID_PACKET;
   }
   x = 1;

   /* decode the length */
   if (in[x] & 0x80) {
      /* valid # of bytes in length are 1,2,3 */
      y = in[x] & 0x7F;
      if ((y == 0) || (y > 3) || ((x + y) > inlen)) {
         return CRYPT_INVALID_PACKET;
      }

      /* read the length in */
      len = 0;
      ++x;
      while (y--) {
         len = (len << 8) | in[x++];
      }
   } else {
      len = in[x++] & 0x7F;
   }

   if (len + x > inlen) {
      return CRYPT_INVALID_PACKET;
   }

   /* proceed to decode */
   for (y = 0; x < inlen; ) {
      /* get first byte */
      tmp = in[x++];
 
      /* count number of bytes */
      for (z = 0; (tmp & 0x80) && (z <= 4); z++, tmp = (tmp << 1) & 0xFF);
      
      if (z > 4 || (x + (z - 1) > inlen)) {
         return CRYPT_INVALID_PACKET;
      }

      /* decode, grab upper bits */
      tmp >>= z;

      /* grab remaining bytes */
      if (z > 1) { --z; }
      while (z-- != 0) {
         if ((in[x] & 0xC0) != 0x80) {
            return CRYPT_INVALID_PACKET;
         }
         tmp = (tmp << 6) | ((wchar_t)in[x++] & 0x3F);
      }

      if (y > *outlen) {
         *outlen = y;
         return CRYPT_BUFFER_OVERFLOW;
      }
      out[y++] = tmp;
   }
   *outlen = y;

   return CRYPT_OK;
}
 
//DAP src/pk/asn1/der/utf8/der_length_utf8_string.c
/** Return the size in bytes of a UTF-8 character
  @param c   The UTF-8 character to measure
  @return    The size in bytes
*/
unsigned long der_utf8_charsize(const wchar_t c)
{
   if (c <= 0x7F) {
      return 1;
   } else if (c <= 0x7FF) {
      return 2;
   } else if (c <= 0xFFFF) {
      return 3;
   } else {
      return 4;
   }
}

/**
  Gets length of DER encoding of UTF8 STRING 
  @param in       The characters to measure the length of
  @param noctets  The number of octets in the string to encode
  @param outlen   [out] The length of the DER encoding for the given string
  @return CRYPT_OK if successful
*/
int der_length_utf8_string(const wchar_t *in, unsigned long noctets, unsigned long *outlen)
{
   unsigned long x, len;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(outlen != NULL);

   len = 0;
   for (x = 0; x < noctets; x++) {
      if (in[x] < 0 || in[x] > 0x10FFFF) {
         return CRYPT_INVALID_ARG;
      }
      len += der_utf8_charsize(in[x]);
   }

   if (len < 128) {
      /* 0C LL DD DD DD ... */
      *outlen = 2 + len;
   } else if (len < 256) {
      /* 0C 81 LL DD DD DD ... */
      *outlen = 3 + len;
   } else if (len < 65536UL) {
      /* 0C 82 LL LL DD DD DD ... */
      *outlen = 4 + len;
   } else if (len < 16777216UL) {
      /* 0C 83 LL LL LL DD DD DD ... */
      *outlen = 5 + len;
   } else {
      return CRYPT_INVALID_ARG;
   }

   return CRYPT_OK;
}

//DAP src/pk/pkcs1/pkcs_1_mgf1.c
/**
   Perform LTC_PKCS #1 MGF1 (internal)
   @param seed        The seed for MGF1
   @param seedlen     The length of the seed
   @param hash_idx    The index of the hash desired
   @param mask        [out] The destination
   @param masklen     The length of the mask desired
   @return CRYPT_OK if successful
*/
int pkcs_1_mgf1(int                  hash_idx,
                const unsigned char *seed, unsigned long seedlen,
                      unsigned char *mask, unsigned long masklen)
{
   unsigned long hLen, x;
   ulong32       counter;
   int           err;
   hash_state    *md;
   unsigned char *buf;
 
   LTC_ARGCHK(seed != NULL);
   LTC_ARGCHK(mask != NULL);

   /* Hard code hash, DAP
     // ensure valid hash 
   if ((err = hash_is_valid(hash_idx)) != CRYPT_OK) { 
      return err;
   }
   */

   /* get hash output size */
   hLen = hash_descriptor[hash_idx].hashsize;

   /* allocate memory */
   md  = static_cast<hash_state*>(malloc(sizeof(hash_state)));
   buf = static_cast<unsigned char*>(malloc(hLen));
   if (md == NULL || buf == NULL) {
      if (md != NULL) {
         XFREE(md);
      }
      if (buf != NULL) {
         XFREE(buf);
      }
      return CRYPT_MEM;
   }

   /* start counter */
   counter = 0;

   while (masklen > 0) {
       /* handle counter */
       STORE32H(counter, buf);
       ++counter;

       /* get hash of seed || counter */
       if ((err = hash_descriptor[hash_idx].init(md)) != CRYPT_OK) {
          goto LBL_ERR;
       }
       if ((err = hash_descriptor[hash_idx].process(md, seed, seedlen)) != CRYPT_OK) {
          goto LBL_ERR;
       }
       if ((err = hash_descriptor[hash_idx].process(md, buf, 4)) != CRYPT_OK) {
          goto LBL_ERR;
       }
       if ((err = hash_descriptor[hash_idx].done(md, buf)) != CRYPT_OK) {
          goto LBL_ERR;
       }

       /* store it */
       for (x = 0; x < hLen && masklen > 0; x++, masklen--) {
          *mask++ = buf[x];
       }
   }

   err = CRYPT_OK;
LBL_ERR:

   XFREE(buf);
   XFREE(md);

   return err;
}

//DAP src/pk/pkcs1/pkcs_1_oaep_encode.c
/**
  LTC_PKCS #1 v2.00 OAEP encode
  @param msg             The data to encode
  @param msglen          The length of the data to encode (octets)
  @param lparam          A session or system parameter (can be NULL)
  @param lparamlen       The length of the lparam data
  @param modulus_bitlen  The bit length of the RSA modulus
  @param prng            An active PRNG state
  @param prng_idx        The index of the PRNG desired
  @param hash_idx        The index of the hash desired
  @param out             [out] The destination for the encoded data
  @param outlen          [in/out] The max size and resulting size of the encoded data
  @return CRYPT_OK if successful
*/
int pkcs_1_oaep_encode(const unsigned char *msg,    unsigned long msglen,
                       const unsigned char *lparam, unsigned long lparamlen,
                             unsigned long modulus_bitlen, prng_state *prng,
                             int           prng_idx,         int  hash_idx,
                             unsigned char *out,    unsigned long *outlen)
{
   unsigned char *DB, *seed, *mask;
   unsigned long hLen, x, y, modulus_len;
   int           err;

   LTC_ARGCHK(msg    != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);

   /* DAP Hardcoded now
   // test valid hash
   if ((err = hash_is_valid(hash_idx)) != CRYPT_OK) { 
      return err;
   }

   // valid prng
   if ((err = prng_is_valid(prng_idx)) != CRYPT_OK) {
      return err;
   }
   */

   hLen        = hash_descriptor[hash_idx].hashsize;
   modulus_len = (modulus_bitlen >> 3) + (modulus_bitlen & 7 ? 1 : 0);

   /* test message size */
   if ((2*hLen >= (modulus_len - 2)) || (msglen > (modulus_len - 2*hLen - 2))) {
      return CRYPT_PK_INVALID_SIZE;
   }

   /* allocate ram for DB/mask/salt of size modulus_len */
   DB   = static_cast<unsigned char*>(malloc(modulus_len));
   mask = static_cast<unsigned char*>(malloc(modulus_len));
   seed = static_cast<unsigned char*>(malloc(hLen));
   if (DB == NULL || mask == NULL || seed == NULL) {
      if (DB != NULL) {
         XFREE(DB);
      }
      if (mask != NULL) {
         XFREE(mask);
      }
      if (seed != NULL) {
         XFREE(seed);
      }
      return CRYPT_MEM;
   }

   /* get lhash */
   /* DB == lhash || PS || 0x01 || M, PS == k - mlen - 2hlen - 2 zeroes */
   x = modulus_len;
   if (lparam != NULL) {
      if ((err = hash_memory(hash_idx, lparam, lparamlen, DB, &x)) != CRYPT_OK) {
         goto LBL_ERR;
      }
   } else {
      /* can't pass hash_memory a NULL so use DB with zero length */
      if ((err = hash_memory(hash_idx, DB, 0, DB, &x)) != CRYPT_OK) {
         goto LBL_ERR;
      }
   }

   /* append PS then 0x01 (to lhash)  */
   x = hLen;
   y = modulus_len - msglen - 2*hLen - 2;
   memset(DB+x, 0, y);
   x += y;

   /* 0x01 byte */
   DB[x++] = 0x01;

   /* message (length = msglen) */
   memcpy(DB+x, msg, msglen);
   x += msglen;

   /* now choose a random seed */
   if (prng_descriptor[prng_idx].read(seed, hLen, prng) != hLen) {
      err = CRYPT_ERROR_READPRNG;
      goto LBL_ERR;
   }

   /* compute MGF1 of seed (k - hlen - 1) */
   if ((err = pkcs_1_mgf1(hash_idx, seed, hLen, mask, modulus_len - hLen - 1)) != CRYPT_OK) {
      goto LBL_ERR;
   }

   /* xor against DB */
   for (y = 0; y < (modulus_len - hLen - 1); y++) {
       DB[y] ^= mask[y]; 
   }

   /* compute MGF1 of maskedDB (hLen) */ 
   if ((err = pkcs_1_mgf1(hash_idx, DB, modulus_len - hLen - 1, mask, hLen)) != CRYPT_OK) {
      goto LBL_ERR;
   }

   /* XOR against seed */
   for (y = 0; y < hLen; y++) {
      seed[y] ^= mask[y];
   }

   /* create string of length modulus_len */
   if (*outlen < modulus_len) {
      *outlen = modulus_len;
      err = CRYPT_BUFFER_OVERFLOW;
      goto LBL_ERR;
   }

   /* start output which is 0x00 || maskedSeed || maskedDB */
   x = 0;
   out[x++] = 0x00;
   memcpy(out+x, seed, hLen);
   x += hLen;
   memcpy(out+x, DB, modulus_len - hLen - 1);
   x += modulus_len - hLen - 1;

   *outlen = x;
    
   err = CRYPT_OK;
LBL_ERR:

   XFREE(seed);
   XFREE(mask);
   XFREE(DB);

   return err;
}

//DAP src/pk/pkcs1/pkcs_1_pss_decode.c
/**
   LTC_PKCS #1 v2.00 PSS decode
   @param  msghash         The hash to verify
   @param  msghashlen      The length of the hash (octets)
   @param  sig             The signature data (encoded data)
   @param  siglen          The length of the signature data (octets)
   @param  saltlen         The length of the salt used (octets)
   @param  hash_idx        The index of the hash desired
   @param  modulus_bitlen  The bit length of the RSA modulus
   @param  res             [out] The result of the comparison, 1==valid, 0==invalid
   @return CRYPT_OK if successful (even if the comparison failed)
*/
int pkcs_1_pss_decode(const unsigned char *msghash, unsigned long msghashlen,
                      const unsigned char *sig,     unsigned long siglen,
                            unsigned long saltlen,  int           hash_idx,
                            unsigned long modulus_bitlen, int    *res)
{
   unsigned char *DB, *mask, *salt, *hash;
   unsigned long x, y, hLen, modulus_len;
   int           err;
   hash_state    md;

   LTC_ARGCHK(msghash != NULL);
   LTC_ARGCHK(res     != NULL);

   /* default to invalid */
   *res = 0;

   /* DAP hard coded 
   // ensure hash is valid 
   if ((err = hash_is_valid(hash_idx)) != CRYPT_OK) {
      return err;
   }
   */

   hLen        = hash_descriptor[hash_idx].hashsize;
   modulus_len = (modulus_bitlen>>3) + (modulus_bitlen & 7 ? 1 : 0);

   /* check sizes */
   if ((saltlen > modulus_len) || 
       (modulus_len < hLen + saltlen + 2) || (siglen != modulus_len)) {
      return CRYPT_PK_INVALID_SIZE;
   }

   /* allocate ram for DB/mask/salt/hash of size modulus_len */
   DB   = static_cast<unsigned char*>(malloc(modulus_len));
   mask = static_cast<unsigned char*>(malloc(modulus_len));
   salt = static_cast<unsigned char*>(malloc(modulus_len));
   hash = static_cast<unsigned char*>(malloc(modulus_len));
   if (DB == NULL || mask == NULL || salt == NULL || hash == NULL) {
      if (DB != NULL) {
         XFREE(DB);
      }
      if (mask != NULL) {
         XFREE(mask);
      }
      if (salt != NULL) {
         XFREE(salt);
      }
      if (hash != NULL) {
         XFREE(hash);
      }
      return CRYPT_MEM;
   }

   /* ensure the 0xBC byte */
   if (sig[siglen-1] != 0xBC) {
      err = CRYPT_INVALID_PACKET;
      goto LBL_ERR;
   }

   /* copy out the DB */
   x = 0;
   memcpy(DB, sig + x, modulus_len - hLen - 1);
   x += modulus_len - hLen - 1;

   /* copy out the hash */
   memcpy(hash, sig + x, hLen);
   x += hLen;

   /* check the MSB */
   if ((sig[0] & ~(0xFF >> ((modulus_len<<3) - (modulus_bitlen-1)))) != 0) {
      err = CRYPT_INVALID_PACKET;
      goto LBL_ERR;
   }

   /* generate mask of length modulus_len - hLen - 1 from hash */
   if ((err = pkcs_1_mgf1(hash_idx, hash, hLen, mask, modulus_len - hLen - 1)) != CRYPT_OK) {
      goto LBL_ERR;
   }

   /* xor against DB */
   for (y = 0; y < (modulus_len - hLen - 1); y++) {
      DB[y] ^= mask[y];
   }
   
   /* now clear the first byte [make sure smaller than modulus] */
   DB[0] &= 0xFF >> ((modulus_len<<3) - (modulus_bitlen-1));

   /* DB = PS || 0x01 || salt, PS == modulus_len - saltlen - hLen - 2 zero bytes */

   /* check for zeroes and 0x01 */
   for (x = 0; x < modulus_len - saltlen - hLen - 2; x++) {
       if (DB[x] != 0x00) {
          err = CRYPT_INVALID_PACKET;
          goto LBL_ERR;
       }
   }

   /* check for the 0x01 */
   if (DB[x++] != 0x01) {
      err = CRYPT_INVALID_PACKET;
      goto LBL_ERR;
   }

   /* M = (eight) 0x00 || msghash || salt, mask = H(M) */
   if ((err = hash_descriptor[hash_idx].init(&md)) != CRYPT_OK) {
      goto LBL_ERR;
   }
   zeromem(mask, 8);
   if ((err = hash_descriptor[hash_idx].process(&md, mask, 8)) != CRYPT_OK) {
      goto LBL_ERR;
   }
   if ((err = hash_descriptor[hash_idx].process(&md, msghash, msghashlen)) != CRYPT_OK) {
      goto LBL_ERR;
   }
   if ((err = hash_descriptor[hash_idx].process(&md, DB+x, saltlen)) != CRYPT_OK) {
      goto LBL_ERR;
   }
   if ((err = hash_descriptor[hash_idx].done(&md, mask)) != CRYPT_OK) {
      goto LBL_ERR;
   }

   /* mask == hash means valid signature */
   if (memcmp(mask, hash, hLen) == 0) {
      *res = 1;
   }

   err = CRYPT_OK;
LBL_ERR:

   XFREE(hash);
   XFREE(salt);
   XFREE(mask);
   XFREE(DB);

   return err;
}

//DAP src/pk/pkcs1/pkcs_1_v1_5_decode.c
/** @brief LTC_PKCS #1 v1.5 decode.
 *
 *  @param msg              The encoded data to decode
 *  @param msglen           The length of the encoded data (octets)
 *  @param block_type       Block type to use in padding (\sa ltc_pkcs_1_v1_5_blocks)
 *  @param modulus_bitlen   The bit length of the RSA modulus
 *  @param out              [out] Destination of decoding
 *  @param outlen           [in/out] The max size and resulting size of the decoding
 *  @param is_valid         [out] Boolean whether the padding was valid
 *
 *  @return CRYPT_OK if successful (even if invalid)
 */
int pkcs_1_v1_5_decode(const unsigned char *msg, 
                             unsigned long  msglen,
                                       int  block_type,
                             unsigned long  modulus_bitlen,
                             unsigned char *out, 
                             unsigned long *outlen,
                                       int *is_valid)
{
  unsigned long modulus_len, ps_len, i;
  int result;

  /* default to invalid packet */
  *is_valid = 0;

  modulus_len = (modulus_bitlen >> 3) + (modulus_bitlen & 7 ? 1 : 0);

  /* test message size */

  if ((msglen > modulus_len) || (modulus_len < 11)) {
    return CRYPT_PK_INVALID_SIZE;
  }

  /* separate encoded message */

  if ((msg[0] != 0x00) || (msg[1] != (unsigned char)block_type)) {
    result = CRYPT_INVALID_PACKET;
    goto bail;
  }

  if (block_type == LTC_LTC_PKCS_1_EME) {
    for (i = 2; i < modulus_len; i++) {
      /* separator */
      if (msg[i] == 0x00) { break; }
    }
    ps_len = i++ - 2;

    if ((i >= modulus_len) || (ps_len < 8)) {
      /* There was no octet with hexadecimal value 0x00 to separate ps from m,
       * or the length of ps is less than 8 octets.
       */
      result = CRYPT_INVALID_PACKET;
      goto bail;
    }
  } else {
    for (i = 2; i < modulus_len - 1; i++) {
       if (msg[i] != 0xFF) { break; }
    }

    /* separator check */
    if (msg[i] != 0) {
      /* There was no octet with hexadecimal value 0x00 to separate ps from m. */
      result = CRYPT_INVALID_PACKET;
      goto bail;
    }

    ps_len = i - 2;
  }

  if (*outlen < (msglen - (2 + ps_len + 1))) {
    *outlen = msglen - (2 + ps_len + 1);
    result = CRYPT_BUFFER_OVERFLOW;
    goto bail;
  }

  *outlen = (msglen - (2 + ps_len + 1));
  memcpy(out, &msg[2 + ps_len + 1], *outlen);

  /* valid packet */
  *is_valid = 1;
  result    = CRYPT_OK;
bail:
  return result;
} /* pkcs_1_v1_5_decode */

//DAP src/pk/pkcs1/pkcs_1_v1_5_encode.c
/*! \brief LTC_PKCS #1 v1.5 encode.
 *
 *  \param msg              The data to encode
 *  \param msglen           The length of the data to encode (octets)
 *  \param block_type       Block type to use in padding (\sa ltc_pkcs_1_v1_5_blocks)
 *  \param modulus_bitlen   The bit length of the RSA modulus
 *  \param prng             An active PRNG state (only for LTC_LTC_PKCS_1_EME)
 *  \param prng_idx         The index of the PRNG desired (only for LTC_LTC_PKCS_1_EME)
 *  \param out              [out] The destination for the encoded data
 *  \param outlen           [in/out] The max size and resulting size of the encoded data
 *
 *  \return CRYPT_OK if successful
 */
int pkcs_1_v1_5_encode(const unsigned char *msg, 
                             unsigned long  msglen,
                                       int  block_type,
                             unsigned long  modulus_bitlen,
                                prng_state *prng, 
                                       int  prng_idx,
                             unsigned char *out, 
                             unsigned long *outlen)
{
  unsigned long modulus_len, ps_len, i;
  unsigned char *ps;
  int result;

  /* valid block_type? */
  if ((block_type != LTC_LTC_PKCS_1_EMSA) &&
      (block_type != LTC_LTC_PKCS_1_EME)) {
     return CRYPT_PK_INVALID_PADDING;
  }

/* DAP, prng hard coded 
  if (block_type == LTC_LTC_PKCS_1_EME) {    // encryption padding, we need a valid PRNG 
    if ((result = prng_is_valid(prng_idx)) != CRYPT_OK) {
       return result;
    }
  }
*/
  modulus_len = (modulus_bitlen >> 3) + (modulus_bitlen & 7 ? 1 : 0);

  /* test message size */
  if ((msglen + 11) > modulus_len) {
    return CRYPT_PK_INVALID_SIZE;
  }

  if (*outlen < modulus_len) {
    *outlen = modulus_len;
    result = CRYPT_BUFFER_OVERFLOW;
    goto bail;
  }

  /* generate an octets string PS */
  ps = &out[2];
  ps_len = modulus_len - msglen - 3;

  if (block_type == LTC_LTC_PKCS_1_EME) {
    /* now choose a random ps */
    if (prng_descriptor[prng_idx].read(ps, ps_len, prng) != ps_len) {
      result = CRYPT_ERROR_READPRNG;
      goto bail;
    }

    /* transform zero bytes (if any) to non-zero random bytes */
    for (i = 0; i < ps_len; i++) {
      while (ps[i] == 0) {
        if (prng_descriptor[prng_idx].read(&ps[i], 1, prng) != 1) {
          result = CRYPT_ERROR_READPRNG;
          goto bail;
        }
      }
    }
  } else {
    memset(ps, 0xFF, ps_len);
  }

  /* create string of length modulus_len */
  out[0]          = 0x00;
  out[1]          = (unsigned char)block_type;  /* block_type 1 or 2 */
  out[2 + ps_len] = 0x00;
  memcpy(&out[2 + ps_len + 1], msg, msglen);
  *outlen = modulus_len;

  result  = CRYPT_OK;
bail:
  return result;
} /* pkcs_1_v1_5_encode */

//DAP src/pk/rsa/rsa_encrypt_key.c
/**
    (LTC_PKCS #1 v2.0) OAEP pad then encrypt
    @param in          The plaintext
    @param inlen       The length of the plaintext (octets)
    @param out         [out] The ciphertext
    @param outlen      [in/out] The max size and resulting size of the ciphertext
    @param lparam      The system "lparam" for the encryption
    @param lparamlen   The length of lparam (octets)
    @param prng        An active PRNG
    @param prng_idx    The index of the desired prng
    @param hash_idx    The index of the desired hash
    @param padding     Type of padding (LTC_LTC_PKCS_1_OAEP or LTC_LTC_PKCS_1_V1_5)
    @param key         The RSA key to encrypt to
    @return CRYPT_OK if successful
*/
int rsa_encrypt_key_ex(const unsigned char *in,     unsigned long inlen,
                             unsigned char *out,    unsigned long *outlen,
                       const unsigned char *lparam, unsigned long lparamlen,
                       prng_state *prng, int prng_idx, int hash_idx, int padding, rsa_key *key)
{
  unsigned long modulus_bitlen, modulus_bytelen, x;
  int           err;

  LTC_ARGCHK(in     != NULL);
  LTC_ARGCHK(out    != NULL);
  LTC_ARGCHK(outlen != NULL);
  LTC_ARGCHK(key    != NULL);

  /* valid padding? */
  if ((padding != LTC_LTC_PKCS_1_V1_5) &&
      (padding != LTC_LTC_PKCS_1_OAEP)) {
    return CRYPT_PK_INVALID_PADDING;
  }

  /* prng hard coded now, DAP
  // valid prng? 
  if ((err = prng_is_valid(prng_idx)) != CRYPT_OK) {
     return err;
  }
  */

  /* Hash Hard coded now, DAP
  if (padding == LTC_LTC_PKCS_1_OAEP) {
    // valid hash?
    if ((err = hash_is_valid(hash_idx)) != CRYPT_OK) {
       return err;
    }
  }
 */

  /* get modulus len in bits */
  modulus_bitlen = mp_count_bits( (key->N));

  /* outlen must be at least the size of the modulus */
  modulus_bytelen = mp_unsigned_bin_size( (key->N));
  if (modulus_bytelen > *outlen) {
     *outlen = modulus_bytelen;
     return CRYPT_BUFFER_OVERFLOW;
  }

  if (padding == LTC_LTC_PKCS_1_OAEP) {
    /* OAEP pad the key */
    x = *outlen;
    if ((err = pkcs_1_oaep_encode(in, inlen, lparam,
                                  lparamlen, modulus_bitlen, prng, prng_idx, hash_idx,
                                  out, &x)) != CRYPT_OK) {
       return err;
    }
  } else {
    /* LTC_PKCS #1 v1.5 pad the key */
    x = *outlen;
    if ((err = pkcs_1_v1_5_encode(in, inlen, LTC_LTC_PKCS_1_EME,
                                  modulus_bitlen, prng, prng_idx,
                                  out, &x)) != CRYPT_OK) {
      return err;
    }
  }

  /* rsa exptmod the OAEP or LTC_PKCS #1 v1.5 pad */
  return ltc_mp.rsa_me(out, x, out, outlen, PK_PUBLIC, key);
}

//DAP src/pk/rsa/rsa_exptmod.c
/** 
   Compute an RSA modular exponentiation 
   @param in         The input data to send into RSA
   @param inlen      The length of the input (octets)
   @param out        [out] The destination 
   @param outlen     [in/out] The max size and resulting size of the output
   @param which      Which exponent to use, e.g. PK_PRIVATE or PK_PUBLIC
   @param key        The RSA key to use 
   @return CRYPT_OK if successful
*/   
int rsa_exptmod(const unsigned char *in,   unsigned long inlen,
                      unsigned char *out,  unsigned long *outlen, int which,
                      rsa_key *key)
{
   void         *tmp, *tmpa, *tmpb;
   unsigned long x;
   int           err;

   LTC_ARGCHK(in     != NULL);
   LTC_ARGCHK(out    != NULL);
   LTC_ARGCHK(outlen != NULL);
   LTC_ARGCHK(key    != NULL);
  
   /* is the key of the right type for the operation? */
   if (which == PK_PRIVATE && (key->type != PK_PRIVATE)) {
      return CRYPT_PK_NOT_PRIVATE;
   }

   /* must be a private or public operation */
   if (which != PK_PRIVATE && which != PK_PUBLIC) {
      return CRYPT_PK_INVALID_TYPE;
   }

   /* init and copy into tmp */
   if ((err = mp_init_multi(&tmp, &tmpa, &tmpb, NULL)) != CRYPT_OK)                                    { return err; }
   if ((err = mp_read_unsigned_bin(tmp, (unsigned char *)in, (int)inlen)) != CRYPT_OK)                 { goto error; }

   /* sanity check on the input */
   if (mp_cmp(key->N, tmp) == LTC_MP_LT) {
      err = CRYPT_PK_INVALID_SIZE;
      goto error;
   }

   /* are we using the private exponent and is the key optimized? */
   if (which == PK_PRIVATE) {
      /* tmpa = tmp^dP mod p */
      if ((err = mp_exptmod(tmp, key->dP, key->p, tmpa)) != CRYPT_OK)                               { goto error; }

      /* tmpb = tmp^dQ mod q */
      if ((err = mp_exptmod(tmp, key->dQ, key->q, tmpb)) != CRYPT_OK)                               { goto error; }

      /* tmp = (tmpa - tmpb) * qInv (mod p) */
      if ((err = mp_sub(tmpa, tmpb, tmp)) != CRYPT_OK)                                              { goto error; }
      if ((err = mp_mulmod(tmp, key->qP, key->p, tmp)) != CRYPT_OK)                                { goto error; }

      /* tmp = tmpb + q * tmp */
      if ((err = mp_mul(tmp, key->q, tmp)) != CRYPT_OK)                                             { goto error; }
      if ((err = mp_add(tmp, tmpb, tmp)) != CRYPT_OK)                                               { goto error; }
   } else {
      /* exptmod it */
      if ((err = mp_exptmod(tmp, key->e, key->N, tmp)) != CRYPT_OK)                                { goto error; }
   }

   /* read it back */
   x = (unsigned long)mp_unsigned_bin_size(key->N);
   if (x > *outlen) {
      *outlen = x;
      err = CRYPT_BUFFER_OVERFLOW;
      goto error;
   }

   /* this should never happen ... */
   if (mp_unsigned_bin_size(tmp) > mp_unsigned_bin_size(key->N)) {
      err = CRYPT_ERROR;
      goto error;
   }
   *outlen = x;

   /* convert it */
   zeromem(out, x);
   if ((err = mp_to_unsigned_bin(tmp, out+(x-mp_unsigned_bin_size(tmp)))) != CRYPT_OK)               { goto error; }

   /* clean up and return */
   err = CRYPT_OK;
error:
   mp_clear_multi(tmp, tmpa, tmpb, NULL);
   return err;
}

//DAP src/pk/rsa/rsa_free.c
/**
  Free an RSA key from memory
  @param key   The RSA key to free
*/
void rsa_free(rsa_key *key)
{
   LTC_ARGCHKVD(key != NULL);
   mp_clear_multi(key->e, key->d, key->N, key->dQ, key->dP, key->qP, key->p, key->q, NULL);
}

//DAP src/pk/rsa/rsa_import.c
/**
  Import an RSAPublicKey or RSAPrivateKey [two-prime only, only support >= 1024-bit keys, defined in LTC_PKCS #1 v2.1]
  @param in      The packet to import from
  @param inlen   It's length (octets)
  @param key     [out] Destination for newly imported key
  @return CRYPT_OK if successful, upon error allocated memory is freed
*/
int rsa_import(const unsigned char *in, unsigned long inlen, rsa_key *key)
{
   int           err;
   void         *zero;
   unsigned char *tmpbuf;
   unsigned long  t, x, y, z, tmpoid[16];
   ltc_asn1_list ssl_pubkey_hashoid[2];
   ltc_asn1_list ssl_pubkey[2];

   LTC_ARGCHK(in          != NULL);
   LTC_ARGCHK(key         != NULL);
   LTC_ARGCHK(ltc_mp.name != NULL);

   /* init key */
   if ((err = mp_init_multi(&key->e, &key->d, &key->N, &key->dQ, 
                            &key->dP, &key->qP, &key->p, &key->q, NULL)) != CRYPT_OK) {
      return err;
   }

   /* see if the OpenSSL DER format RSA public key will work */
   tmpbuf = static_cast<unsigned char*>(calloc(1, MAX_RSA_SIZE*8));
   if (tmpbuf == NULL) {
       err = CRYPT_MEM;
       goto LBL_ERR;
   }

   /* this includes the internal hash ID and optional params (NULL in this case) */
   LTC_SET_ASN1(ssl_pubkey_hashoid, 0, LTC_ASN1_OBJECT_IDENTIFIER, tmpoid,                sizeof(tmpoid)/sizeof(tmpoid[0]));   
   LTC_SET_ASN1(ssl_pubkey_hashoid, 1, LTC_ASN1_NULL,              NULL,                  0);

   /* the actual format of the SSL DER key is odd, it stores a RSAPublicKey in a **BIT** string ... so we have to extract it
      then proceed to convert bit to octet 
    */
   LTC_SET_ASN1(ssl_pubkey, 0,         LTC_ASN1_SEQUENCE,          &ssl_pubkey_hashoid,   2);
   LTC_SET_ASN1(ssl_pubkey, 1,         LTC_ASN1_BIT_STRING,        tmpbuf,                MAX_RSA_SIZE*8);

   if (der_decode_sequence(in, inlen,
                           ssl_pubkey, 2UL) == CRYPT_OK) {

      /* ok now we have to reassemble the BIT STRING to an OCTET STRING.  Thanks OpenSSL... */
      for (t = y = z = x = 0; x < ssl_pubkey[1].size; x++) {
          y = (y << 1) | tmpbuf[x];
          if (++z == 8) {
             tmpbuf[t++] = (unsigned char)y;
             y           = 0;
             z           = 0;
          }
      }

      /* now it should be SEQUENCE { INTEGER, INTEGER } */
      if ((err = der_decode_sequence_multi(tmpbuf, t,
                                           LTC_ASN1_INTEGER, 1UL, key->N, 
                                           LTC_ASN1_INTEGER, 1UL, key->e, 
                                           LTC_ASN1_EOL,     0UL, NULL)) != CRYPT_OK) {
         XFREE(tmpbuf);
         goto LBL_ERR;
      }
      XFREE(tmpbuf);
      key->type = PK_PUBLIC;
      return CRYPT_OK;
   }
   XFREE(tmpbuf);

   /* not SSL public key, try to match against LTC_PKCS #1 standards */
   if ((err = der_decode_sequence_multi(in, inlen, 
                                  LTC_ASN1_INTEGER, 1UL, key->N, 
                                  LTC_ASN1_EOL,     0UL, NULL)) != CRYPT_OK) {
      goto LBL_ERR;
   }

   if (mp_cmp_d(key->N, 0) == LTC_MP_EQ) {
      if ((err = mp_init(&zero)) != CRYPT_OK) { 
         goto LBL_ERR;
      }
      /* it's a private key */
      if ((err = der_decode_sequence_multi(in, inlen, 
                          LTC_ASN1_INTEGER, 1UL, zero, 
                          LTC_ASN1_INTEGER, 1UL, key->N, 
                          LTC_ASN1_INTEGER, 1UL, key->e,
                          LTC_ASN1_INTEGER, 1UL, key->d, 
                          LTC_ASN1_INTEGER, 1UL, key->p, 
                          LTC_ASN1_INTEGER, 1UL, key->q, 
                          LTC_ASN1_INTEGER, 1UL, key->dP,
                          LTC_ASN1_INTEGER, 1UL, key->dQ, 
                          LTC_ASN1_INTEGER, 1UL, key->qP, 
                          LTC_ASN1_EOL,     0UL, NULL)) != CRYPT_OK) {
         mp_clear(zero);
         goto LBL_ERR;
      }
      mp_clear(zero);
      key->type = PK_PRIVATE;
   } else if (mp_cmp_d(key->N, 1) == LTC_MP_EQ) {
      /* we don't support multi-prime RSA */
      err = CRYPT_PK_INVALID_TYPE;
      goto LBL_ERR;
   } else {
      /* it's a public key and we lack e */
      if ((err = der_decode_sequence_multi(in, inlen, 
                                     LTC_ASN1_INTEGER, 1UL, key->N, 
                                     LTC_ASN1_INTEGER, 1UL, key->e, 
                                     LTC_ASN1_EOL,     0UL, NULL)) != CRYPT_OK) {
         goto LBL_ERR;
      }
      key->type = PK_PUBLIC;
   }
   return CRYPT_OK;
LBL_ERR:
   mp_clear_multi(key->d,  key->e, key->N, key->dQ, key->dP, key->qP, key->p, key->q, NULL);
   return err;
}

//DAP src/pk/rsa/rsa_verify_hash.c
/**
  LTC_PKCS #1 de-sign then v1.5 or PSS depad
  @param sig              The signature data
  @param siglen           The length of the signature data (octets)
  @param hash             The hash of the message that was signed
  @param hashlen          The length of the hash of the message that was signed (octets)
  @param padding          Type of padding (LTC_LTC_PKCS_1_PSS or LTC_LTC_PKCS_1_V1_5)
  @param hash_idx         The index of the desired hash
  @param saltlen          The length of the salt used during signature
  @param stat             [out] The result of the signature comparison, 1==valid, 0==invalid
  @param key              The public RSA key corresponding to the key that performed the signature
  @return CRYPT_OK on success (even if the signature is invalid)
*/
int rsa_verify_hash_ex(const unsigned char *sig,      unsigned long siglen,
                       const unsigned char *hash,     unsigned long hashlen,
                             int            padding,
                             int            hash_idx, unsigned long saltlen,
                             int           *stat,     rsa_key      *key)
{
  unsigned long modulus_bitlen, modulus_bytelen, x;
  int           err;
  unsigned char *tmpbuf;

  LTC_ARGCHK(hash  != NULL);
  LTC_ARGCHK(sig   != NULL);
  LTC_ARGCHK(stat  != NULL);
  LTC_ARGCHK(key   != NULL);

  /* default to invalid */
  *stat = 0;

  /* valid padding? */

  if ((padding != LTC_LTC_PKCS_1_V1_5) &&
      (padding != LTC_LTC_PKCS_1_PSS)) {
    return CRYPT_PK_INVALID_PADDING;
  }

  if (padding == LTC_LTC_PKCS_1_PSS) {
    /* DAP Hard coded
    // valid hash ?
    if ((err = hash_is_valid(hash_idx)) != CRYPT_OK) {
       return err;
    }
    */
  }

  /* get modulus len in bits */
  modulus_bitlen = mp_count_bits( (key->N));

  /* outlen must be at least the size of the modulus */
  modulus_bytelen = mp_unsigned_bin_size( (key->N));
  if (modulus_bytelen != siglen) {
     return CRYPT_INVALID_PACKET;
  }

  /* allocate temp buffer for decoded sig */
  tmpbuf = static_cast<unsigned char*>(malloc(siglen));
  if (tmpbuf == NULL) {
     return CRYPT_MEM;
  }

  /* RSA decode it  */
  x = siglen;
  if ((err = ltc_mp.rsa_me(sig, siglen, tmpbuf, &x, PK_PUBLIC, key)) != CRYPT_OK) {
     XFREE(tmpbuf);
     return err;
  }

  /* make sure the output is the right size */
  if (x != siglen) {
     XFREE(tmpbuf);
     return CRYPT_INVALID_PACKET;
  }

  if (padding == LTC_LTC_PKCS_1_PSS) {
    /* PSS decode and verify it */
    err = pkcs_1_pss_decode(hash, hashlen, tmpbuf, x, saltlen, hash_idx, modulus_bitlen, stat);
  } else {
    /* LTC_PKCS #1 v1.5 decode it */
    unsigned char *out;
    unsigned long outlen, loid[16];
    int           decoded;
    ltc_asn1_list digestinfo[2], siginfo[2];

    /* not all hashes have OIDs... so sad */
    if (hash_descriptor[hash_idx].OIDlen == 0) {
       err = CRYPT_INVALID_ARG;
       goto bail_2;
    }

    /* allocate temp buffer for decoded hash */
    outlen = ((modulus_bitlen >> 3) + (modulus_bitlen & 7 ? 1 : 0)) - 3;
    out    = static_cast<unsigned char*>(malloc(outlen));
    if (out == NULL) {
      err = CRYPT_MEM;
      goto bail_2;
    }

    if ((err = pkcs_1_v1_5_decode(tmpbuf, x, LTC_LTC_PKCS_1_EMSA, modulus_bitlen, out, &outlen, &decoded)) != CRYPT_OK) {
      XFREE(out);       
      goto bail_2;
    }

    /* now we must decode out[0...outlen-1] using ASN.1, test the OID and then test the hash */
    /* construct the SEQUENCE 
      SEQUENCE {
         SEQUENCE {hashoid OID
                   blah    NULL
         }
         hash    OCTET STRING 
      }
   */
    LTC_SET_ASN1(digestinfo, 0, LTC_ASN1_OBJECT_IDENTIFIER, loid, sizeof(loid)/sizeof(loid[0]));
    LTC_SET_ASN1(digestinfo, 1, LTC_ASN1_NULL,              NULL,                          0);
    LTC_SET_ASN1(siginfo,    0, LTC_ASN1_SEQUENCE,          digestinfo,                    2);
    LTC_SET_ASN1(siginfo,    1, LTC_ASN1_OCTET_STRING,      tmpbuf,                        siglen);
   
    if ((err = der_decode_sequence(out, outlen, siginfo, 2)) != CRYPT_OK) {
       XFREE(out);
       goto bail_2;
    }

    /* test OID */
    if ((digestinfo[0].size == hash_descriptor[hash_idx].OIDlen) &&
        (memcmp(digestinfo[0].data, hash_descriptor[hash_idx].OID, sizeof(unsigned long) * hash_descriptor[hash_idx].OIDlen) == 0) &&
        (siginfo[1].size == hashlen) &&
        (memcmp(siginfo[1].data, hash, hashlen) == 0)) {
       *stat = 1;
    }

    XFREE(out);
  }

bail_2:
  XFREE(tmpbuf);
  return err;
}

//DAP src/prngs/rng_get_bytes.c
/* on *NIX read /dev/random */
static unsigned long rng_nix(unsigned char *buf, unsigned long len, 
                             void (*callback)(void))
{
    FILE *f;
    unsigned long x;
    f = fopen("/dev/urandom", "rb");
    if (f == NULL)
       f = fopen("/dev/random", "rb");

    if (f == NULL) {
       return 0;
    }
    
    /* disable buffering */
    if (setvbuf(f, NULL, _IONBF, 0) != 0) {
       fclose(f);
       return 0;
    }   
 
    x = (unsigned long)fread(buf, 1, (size_t)len, f);
    fclose(f);
    return x;
}

/**
  Read the system RNG
  @param out       Destination
  @param outlen    Length desired (octets)
  @param callback  Pointer to void function to act as "callback" when RNG is slow.  This can be NULL
  @return Number of octets read
*/     
unsigned long rng_get_bytes(unsigned char *out, unsigned long outlen, 
                            void (*callback)(void))
{
   unsigned long x;

   LTC_ARGCHK(out != NULL);

   x = rng_nix(out, outlen, callback);   if (x != 0) { return x; }
   return 0;
}

//DAP src/prngs/sprng.c
/* A secure PRNG using the RNG functions.  Basically this is a
 * wrapper that allows you to use a secure RNG as a PRNG
 * in the various other functions.
 */

struct ltc_prng_descriptor prng_descriptor[TAB_SIZE] = {
  { NULL, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL }
};

const struct ltc_prng_descriptor sprng_desc =
{
    "sprng", 0,
    &sprng_start,
    &sprng_add_entropy,
    &sprng_ready,
    &sprng_read,
    &sprng_done,
    &sprng_export,
    &sprng_import,
    &sprng_test
};

/**
  Start the PRNG
  @param prng     [out] The PRNG state to initialize
  @return CRYPT_OK if successful
*/  
int sprng_start(prng_state *prng)
{
   return CRYPT_OK;  
}

/**
  Add entropy to the PRNG state
  @param in       The data to add
  @param inlen    Length of the data to add
  @param prng     PRNG state to update
  @return CRYPT_OK if successful
*/  
int sprng_add_entropy(const unsigned char *in, unsigned long inlen, prng_state *prng)
{
   return CRYPT_OK;
}

/**
  Make the PRNG ready to read from
  @param prng   The PRNG to make active
  @return CRYPT_OK if successful
*/  
int sprng_ready(prng_state *prng)
{
   return CRYPT_OK;
}

/**
  Read from the PRNG
  @param out      Destination
  @param outlen   Length of output
  @param prng     The active PRNG to read from
  @return Number of octets read
*/  
unsigned long sprng_read(unsigned char *out, unsigned long outlen, prng_state *prng)
{
   LTC_ARGCHK(out != NULL);
   return rng_get_bytes(out, outlen, NULL);
}

/**
  Terminate the PRNG
  @param prng   The PRNG to terminate
  @return CRYPT_OK if successful
*/  
int sprng_done(prng_state *prng)
{
   return CRYPT_OK;
}

/**
  Export the PRNG state
  @param out       [out] Destination
  @param outlen    [in/out] Max size and resulting size of the state
  @param prng      The PRNG to export
  @return CRYPT_OK if successful
*/  
int sprng_export(unsigned char *out, unsigned long *outlen, prng_state *prng)
{
   LTC_ARGCHK(outlen != NULL);

   *outlen = 0;
   return CRYPT_OK;
}
 
/**
  Import a PRNG state
  @param in       The PRNG state
  @param inlen    Size of the state
  @param prng     The PRNG to import
  @return CRYPT_OK if successful
*/  
int sprng_import(const unsigned char *in, unsigned long inlen, prng_state *prng)
{
   return CRYPT_OK;
}

/**
  PRNG self-test
  @return CRYPT_OK if successful, CRYPT_NOP if self-testing has been disabled
*/  
int sprng_test(void)
{
   return CRYPT_OK;
}

//DAP 
