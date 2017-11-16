/**
 * PForDelta compression, from
 * @inproceedings{Zukowski:2006:SRC:1129754.1129919,
     author = {Zukowski, Marcin and Heman, Sandor and Nes, Niels and Boncz, Peter},
     title = {Super-Scalar RAM-CPU Cache Compression},
     booktitle = {Proceedings of the 22Nd International Conference on Data Engineering},
     series = {ICDE '06},
     year = {2006},
     isbn = {0-7695-2570-9},
     pages = {59--},
     url = {http://dx.doi.org/10.1109/ICDE.2006.150},
     doi = {10.1109/ICDE.2006.150},
     acmid = {1129919},
     publisher = {IEEE Computer Society},
     address = {Washington, DC, USA},
    }
 * as implemented by Shuai Ding and Tosten Suel in the Block-Max WAND method described in
 * @inproceedings{Ding:2011:FTD:2009916.2010048,
     author = {Ding, Shuai and Suel, Torsten},
     title = {Faster Top-k Document Retrieval Using Block-max Indexes},
     booktitle = {Proceedings of the 34th International ACM SIGIR Conference on Research and Development in Information Retrieval},
     series = {SIGIR '11},
     year = {2011},
     isbn = {978-1-4503-0757-4},
     location = {Beijing, China},
     pages = {993--1002},
     numpages = {10},
     url = {http://doi.acm.org/10.1145/2009916.2010048},
     doi = {10.1145/2009916.2010048},
     acmid = {2010048},
     publisher = {ACM},
     address = {New York, NY, USA},
     keywords = {block-max index, early termination, inverted index, ir query processing, top-k query processing},
    }
 * http://www.jinruhe.com/.
 */

#ifndef BLOCK_SIZE
    #define BLOCK_SIZE 64
#endif
#define BS BLOCK_SIZE
#define FRAC 0.10
#define S 16

void pack(unsigned int *v, unsigned int b, unsigned int n, unsigned int *w);
int pack_encode(unsigned int **w, unsigned int *p, int num);

void unpack2(unsigned int *p, unsigned int *w);
void unpack3(unsigned int *p, unsigned int *w);
void unpack4(unsigned int *p, unsigned int *w);
void unpack5(unsigned int *p, unsigned int *w);
void unpack6(unsigned int *p, unsigned int *w);
void unpack7(unsigned int *p, unsigned int *w);
void unpack8(unsigned int *p, unsigned int *w);
void unpack9(unsigned int *p, unsigned int *w);
void unpack10(unsigned int *p, unsigned int *w);
void unpack12(unsigned int *p, unsigned int *w);
void unpack16(unsigned int *p, unsigned int *w);
void unpack20(unsigned int *p, unsigned int *w);
void unpack32(unsigned int *p, unsigned int *w);

typedef void (*pf)(unsigned int *p, unsigned int *w);
unsigned int *pack_decode(unsigned int *_p, unsigned int *_w, int flag,
                          int pref, int base);

void prefix(unsigned int *_p, unsigned int s);
void prefix2(unsigned int *p);
void prefix3(unsigned int *p);
