/* The MIT License

   Copyright (c) 20082-2012 by Heng Li <lh3@me.com>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#include "public_parameters.h"
#include "HardClipReads.h"
#include <fstream>
#include <map>
using namespace std;

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <zlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <math.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

typedef struct {
	int n, m;
	uint64_t *a;
} reglist_t;

#include "khash.h"
KHASH_MAP_INIT_STR(reg, reglist_t)

typedef kh_reg_t reghash_t;

reghash_t *stk_reg_read(const char *fn)
{
	reghash_t *h = kh_init(reg);
	gzFile fp;
	kstream_t *ks;
	int dret;
	kstring_t *str;
	// read the list 
	str = (kstring_t *)calloc(1, sizeof(kstring_t));
	fp = strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	ks = ks_init(fp);
	while (ks_getuntil(ks, 0, str, &dret) >= 0) {
		int beg = -1, end = -1;
		reglist_t *p;
		khint_t k = kh_get(reg, h, str->s);
		if (k == kh_end(h)) {
			int ret;
			char *s = strdup(str->s);
			k = kh_put(reg, h, s, &ret);
			memset(&kh_val(h, k), 0, sizeof(reglist_t));
		}
		p = &kh_val(h, k);
		if (dret != '\n') {
			if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
				beg = atoi(str->s);
				if (dret != '\n') {
					if (ks_getuntil(ks, 0, str, &dret) > 0 && isdigit(str->s[0])) {
						end = atoi(str->s);
						if (end < 0) end = -1;
					}
				}
			}
		}
		// skip the rest of the line
		if (dret != '\n') while ((dret = ks_getc(ks)) > 0 && dret != '\n');
		if (end < 0 && beg > 0) end = beg, beg = beg - 1; // if there is only one column
		if (beg < 0) beg = 0, end = INT_MAX;
		if (p->n == p->m) {
			p->m = p->m? p->m<<1 : 4;
			p->a = (uint64_t*)realloc(p->a, p->m * 8);
		}
		p->a[p->n++] = (uint64_t)beg<<32 | end;
	}
	ks_destroy(ks);
	gzclose(fp);
	free(str->s); free(str);
	return h;
}

void stk_reg_destroy(reghash_t *h)
{
	khint_t k;
	if (h == 0) return;
	for (k = 0; k < kh_end(h); ++k) {
		if (kh_exist(h, k)) {
			free(kh_val(h, k).a);
			free((char*)kh_key(h, k));
		}
	}
	kh_destroy(reg, h);
}

/* constant table */

unsigned char seq_nt16_table[256] = {
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15 /*'-'*/,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9,  0,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

char *seq_nt16_rev_table = "XACMGRSVTWYHKDBN";
unsigned char seq_nt16to4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
unsigned char seq_nt16comp_table[] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
int bitcnt_table[] = { 4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
char comp_tab[] = {
	  0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
	 16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
	 32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
	 48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
	 64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
	 64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

static void stk_printstr(const kstring_t *s, unsigned line_len)
{
	if (line_len != UINT_MAX && line_len != 0) {
		int i, rest = s->l;
		for (i = 0; i < s->l; i += line_len, rest -= line_len) {
			putchar('\n');
			if (rest > line_len) fwrite(s->s + i, 1, line_len, stdout);
			else fwrite(s->s + i, 1, rest, stdout);
		}
		putchar('\n');
	} else {
		putchar('\n');
		puts(s->s);
	}
}

void stk_printseq(const kseq_t *s, int line_len)
{
	putchar(s->qual.l? '@' : '>');
	fputs(s->name.s, stdout);
	if (s->comment.l) {
		putchar(' '); fputs(s->comment.s, stdout);
	}
	stk_printstr(&s->seq, line_len);
	if (s->qual.l) {
		putchar('+');
		stk_printstr(&s->qual, line_len);
	}
}

/* 
   64-bit Mersenne Twister pseudorandom number generator. Adapted from:

     http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/C-LANG/mt19937-64.c

   which was written by Takuji Nishimura and Makoto Matsumoto and released
   under the 3-clause BSD license. 
*/

typedef uint64_t krint64_t;

struct _krand_t;
typedef struct _krand_t krand_t;

#define KR_NN 312
#define KR_MM 156
#define KR_UM 0xFFFFFFFF80000000ULL /* Most significant 33 bits */
#define KR_LM 0x7FFFFFFFULL /* Least significant 31 bits */

struct _krand_t {
	int mti;
	krint64_t mt[KR_NN];
};

static void kr_srand0(krint64_t seed, krand_t *kr)
{
	kr->mt[0] = seed;
	for (kr->mti = 1; kr->mti < KR_NN; ++kr->mti) 
		kr->mt[kr->mti] = 6364136223846793005ULL * (kr->mt[kr->mti - 1] ^ (kr->mt[kr->mti - 1] >> 62)) + kr->mti;
}

krand_t *kr_srand(krint64_t seed)
{
	krand_t *kr;
	kr = (krand_t *)malloc(sizeof(krand_t));
	kr_srand0(seed, kr);
	return kr;
}

krint64_t kr_rand(krand_t *kr)
{
	krint64_t x;
	static const krint64_t mag01[2] = { 0, 0xB5026F5AA96619E9ULL };
    if (kr->mti >= KR_NN) {
		int i;
		if (kr->mti == KR_NN + 1) kr_srand0(5489ULL, kr);
        for (i = 0; i < KR_NN - KR_MM; ++i) {
            x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
            kr->mt[i] = kr->mt[i + KR_MM] ^ (x>>1) ^ mag01[(int)(x&1)];
        }
        for (; i < KR_NN - 1; ++i) {
            x = (kr->mt[i] & KR_UM) | (kr->mt[i+1] & KR_LM);
            kr->mt[i] = kr->mt[i + (KR_MM - KR_NN)] ^ (x>>1) ^ mag01[(int)(x&1)];
        }
        x = (kr->mt[KR_NN - 1] & KR_UM) | (kr->mt[0] & KR_LM);
        kr->mt[KR_NN - 1] = kr->mt[KR_MM - 1] ^ (x>>1) ^ mag01[(int)(x&1)];
        kr->mti = 0;
    }
    x = kr->mt[kr->mti++];
    x ^= (x >> 29) & 0x5555555555555555ULL;
    x ^= (x << 17) & 0x71D67FFFEDA60000ULL;
    x ^= (x << 37) & 0xFFF7EEE000000000ULL;
    x ^= (x >> 43);
    return x;
}

#define kr_drand(_kr) ((kr_rand(_kr) >> 11) * (1.0/9007199254740992.0))


/* subseq */
int stk_subseq(const char* freads, const char* fid, const char* fout_raw)
{
	khash_t(reg) *h = kh_init(reg);
	gzFile fp;
	kseq_t *seq;
	int l, i, j, is_tab = 0, line = 0;
	khint_t k;

	h = stk_reg_read(fid);
	FILE * pFile;
	pFile = fopen (fout_raw,"w");
	// subseq 
	fp=gzopen(freads, "r");
	seq = kseq_init(fp);
	while ((l = kseq_read(seq)) >= 0) {
		reglist_t *p;
		k = kh_get(reg, h, seq->name.s);
		if (k == kh_end(h)) continue;
		p = &kh_val(h, k);
		for (i = 0; i < p->n; ++i) {
			int beg = p->a[i]>>32, end = p->a[i];
			if (beg >= seq->seq.l) {
				fprintf(stderr, "[subseq] %s: %d >= %ld\n", seq->name.s, beg, seq->seq.l);
				continue;
			}
			if (end > seq->seq.l) end = seq->seq.l;
			
			if (is_tab == 0) {
				//printf("%c%s", seq->qual.l == seq->seq.l? '@' : '>', seq->name.s);
				fprintf(pFile, "%s", seq->name.s);
				if (beg > 0 || (int)p->a[i] != INT_MAX) {
					if (end == INT_MAX) {
						if (beg) fprintf(pFile,":%d", beg+1);
					} else fprintf(pFile,":%d-%d", beg+1, end);
				}
			} else fprintf(pFile,"%s\t%d\t", seq->name.s, beg + 1);

			if (end > seq->seq.l) end = seq->seq.l;
			for (j = 0; j < end - beg; ++j) {
				if (is_tab == 0 && (j == 0 || (line > 0 && j % line == 0))) fprintf(pFile,"%c",'\t');//putchar('\t'); /*putchar('\n');*/
				//putchar(seq->seq.s[j + beg]);
				fprintf(pFile,"%c",seq->seq.s[j + beg]);
			}
			//putchar('\n');
			if (seq->qual.l != seq->seq.l || is_tab) continue;
			//printf("+");
			for (j = 0; j < end - beg; ++j) {
				if (j == 0 || (line > 0 && j % line == 0)) {/*putchar('\n');*/}
				//putchar(seq->qual.s[j + beg]);
			}
			fprintf(pFile,"%c",'\n');//putchar('\n');
		}
	}

	fclose(pFile);
	// free
	kseq_destroy(seq);
	gzclose(fp);
	stk_reg_destroy(h);
	return 0;
}

void HardClipReads::traceRawReadById(string fleft_fastq, string fright_fastq)
{
	//trace left hard-clip reads
	stk_subseq(fleft_fastq.c_str(), fname_left_hclip.c_str(), fname_lhclip_raw.c_str());
	
	//trace right hard-clip reads
	stk_subseq(fright_fastq.c_str(), fname_right_hclip.c_str(), fname_rhclip_raw.c_str());
}
