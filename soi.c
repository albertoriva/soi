#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define NDATA 1000
#define NITERS 100000
#define MAXCNT 16
#define PNGFILE "venn.png"

#define FMTDEF 0
#define FMTTAB 1

struct Params {
  int *buf;
  int ndata;
  int niters;
  int nsets;
  int nA, nB, nC, nD;
  int *counts, *targets, *successes; /* Vectors of MAXCNT elements holding counts of intersection occurrences, target intersection, successes */
  float *expected;
  int outfmt;
  char *Rfile;
  char *imgfile;
} params;

struct Params *P;

/*                     0    1    2    3     4    5     6     7      8    9     10    11     12    13     14     15 */
char *KEYS[MAXCNT] = {"N", "A", "B", "AB", "C", "AC", "BC", "ABC", "D", "AD", "BD", "ABD", "CD", "ACD", "BCD", "ABCD"};

char *idxToKey(int idx) {
  if ((idx >= 0) && (idx < MAXCNT)) {
    return KEYS[idx];
  } else {
    fprintf(stderr, "Cannot decode index %d!\n", idx);
    exit(1);
  }
}

int keyToIdx(char *key) {
  int i;

  for (i = 0; i < MAXCNT; i++) {
    if (!strcmp(key, KEYS[i])) {
      return i;
    }
  }
  return -1;
}

void verifyIdxKey() {
  /* For debugging only */
  int i;
  char *key;
  for (i = 0; i < MAXCNT; i++) {
    key = idxToKey(i);
    printf("%d -> %s -> %d\n", i, key, keyToIdx(key));
  }
}

void parsePair(char *s, int *data) {
  /* parse the string S which should be of the form X=n, where X is one of A, B, C, D, or a pair, or a triple, or ABCD. Returns the 
     bitwise value corresponding to X in data[0] and n in data[1]. */
  int i, len;
  char code[10], numb[10];

  len = strlen(s);
  for (i = 0; i < len; i++) {
    if (s[i] == '=') break;
  }
  if (i == len) {
    fprintf(stderr, "Argument %s should be of the form S=N.\n", s);
    exit(1);
  }
  strcpy(code, s);
  code[i] = '\0';
  data[0] = keyToIdx(code);
  data[1] = atoi(s+i+1);
}

void zerofill(int *vec, int n) {
  int i;
  for (i = 0; i < n; i++) {
    vec[i] = 0;
  }
}

void dump() {
  int i;
  printf("N: %d\n", P->ndata);
  printf("Iterations: %d\n", P->niters);
  printf("nA: %d\n", P->nA);
  printf("nB: %d\n", P->nB);
  printf("nC: %d\n", P->nC);
  printf("nD: %d\n", P->nD);

  printf("Targets: ");
  for (i = 0; i < MAXCNT; i++) {
    printf("%d, ", P->targets[i]);
  }
  printf("\n");
}

void usage() {
  printf("Usage: soi [-h] [-t] [-i I] [-n N] specs...\n\n");
  printf("This program evaluates whether the intersection between different sets is statistically significant.\n\n");
  printf("`specs' consists of entries of the form S=n, where S describes a set or an intersection, and n is the\n");
  printf("number of elements it contains. S can be one of: N, A, B, C, D, AB, AC, AD, BC, BD, CD, ABC, ABD, ACD,\n");
  printf("BCD, ABCD. N indicates the total number of items in the domain (can also be set with the -n option).\n");
  printf("\nExample: soi N=1000 A=100 B=200 AB=30\n");
  printf("\nOptions:\n");
  printf(" -h         | Print this help message.\n");
  printf(" -t         | Output results in tab-delimited format.\n");
  printf(" -n N       | Set the total number of items to N (default: %d).\n", NDATA);
  printf(" -i I       | Set the number of iterations to I (default: %d).\n", NITERS);
  printf(" -r rfile   | Write an R script to plot the specified sets to `rfile'.\n");
  printf(" -i imgfile | The R script will save image to png file `imgfile' (default: %s).\n", PNGFILE);
}

void setExpected() {
  float fA, fB, fC, fD, *e;

  e = malloc(MAXCNT*sizeof(float));
  fA = 1.0*P->nA/P->ndata;
  fB = 1.0*P->nB/P->ndata;
  fC = 1.0*P->nC/P->ndata;
  fD = 1.0*P->nD/P->ndata;
  e[3] = fA*fB;
  e[5] = fA*fC;
  e[6] = fB*fC;
  e[7] = fA*fB*fC;
  e[9] = fA*fD;
  e[10] = fB*fD;
  e[11] = fA*fB*fD;
  e[12] = fC*fD;
  e[13] = fA*fC*fD;
  e[14] = fB*fC*fD;
  e[15] = fA*fB*fC*fD;
  P->expected = e;
}

void initialize(int argc, char *argv[]) {
  int i;
  int data[2];
  char *a, *next;

  next = "";

  P = malloc(sizeof(params));
  P->outfmt = FMTDEF;
  P->counts = malloc(MAXCNT*sizeof(int)); /* Counters for occurrence of each intersection in a single iteration */
  P->targets = malloc(MAXCNT*sizeof(int)); /* Observed size of each intersection */
  P->successes = malloc(MAXCNT*sizeof(int)); /* How many times each intersection count is above target */
  P->ndata = NDATA;
  P->niters = NITERS;
  P->nsets = 0;
  P->nA = 0;
  P->nB = 0;
  P->nC = 0;
  P->nD = 0;

  P->Rfile = "";
  P->imgfile = PNGFILE;

  zerofill(P->targets, MAXCNT);

  for (i = 1; i < argc; i++) {
    a = argv[i];

    if (!strcmp(a, "-h")) {
      usage();
      exit(0);
    } else if (!strcmp(next, "-i")) {
      P->niters = atoi(a);
      next = "";
    } else if (!strcmp(next, "-n")) {
      P->ndata = atoi(a);
      next = "";
    } else if (!strcmp(next, "-r")) {
      P->Rfile = a;
      next = "";
    } else if (!strcmp(a, "-i")) {
      next = a;
    } else if (!strcmp(a, "-n")) {
      next = a;
    } else if (!strcmp(a, "-r")) {
      next = a;
    } else if (!strcmp(a, "-t")) {
      P->outfmt = FMTTAB;
    } else {
      parsePair(a, data);
      if (data[0] == 0) {
	P->ndata = data[1];
      } else if (data[0] == 1) {
	P->nA = data[1];
      } else if (data[0] == 2) {
	P->nB = data[1];
      } else if (data[0] == 4) {
	P->nC = data[1];
      } else if (data[0] == 8) {
	P->nD = data[1];
      } else {
	P->targets[data[0]] = data[1];
      }
    }
  }
  P->buf = malloc(P->ndata*sizeof(int));
  setExpected();
}

/* Returns number of non-empty sets. */
int num_sets() {
  int n = 0;
  if (P->nA > 0) n++;
  if (P->nB > 0) n++;
  if (P->nC > 0) n++;
  if (P->nD > 0) n++;
  return n;
}
      
void randomFill(int na, int v) {
  int *buf;
  int N, r, found=0;
  
  buf = P->buf;
  N = P->ndata;

  while (found < na) {
    r = rand() % N;
    if ((buf[r] & v) == 0) {
      buf[r] += v;
      found++;
      /* printf("%d, %d\n", buf[r], found); */
    }
  }
  /* printf("Filled with %dx%d\n", na, v); */
}

int countIntersecting() {
  int i;

  zerofill(P->buf, P->ndata);

  if (P->nA > 0) {
    randomFill(P->nA, 1);
  }
  if (P->nB > 0) {
    randomFill(P->nB, 2);
  }
  if (P->nC > 0) {
    randomFill(P->nC, 4);
  }
  if (P->nD > 0) {
    randomFill(P->nD, 8);
  }

  /* zero counts */
  zerofill(P->counts, MAXCNT);

  /* count how many times we see each intersection */
  for (i = 0; i < P->ndata; i++) {
    P->counts[P->buf[i]]++;
  }

  /* print it 
  for (i = 0; i < MAXCNT; i++) {
    printf("%d: %d\n", i, P->counts[i]);
  }
  */
}

/* Functions to write R script */

void writeRscript2(FILE *f) {
  fprintf(f, "venn.plot <- draw.pairwise.venn(\n");
  fprintf(f, "area1 = %d,\n", P->nA);
  fprintf(f, "area2 = %d,\n", P->nB);
  fprintf(f, "cross.area = %d,\n", P->targets[3]);
  fprintf(f, "category = c(\"A\", \"B\"),\n");
  fprintf(f, "fill = c(\"orange\", \"blue\"),\n");
  fprintf(f, "lty = \"solid\",\n");
  fprintf(f, "euler.d = TRUE,\n");
  fprintf(f, "scaled = TRUE\n");
  fprintf(f, ");\n\n");
}

void writeRscript3(FILE *f) {
  fprintf(f, "venn.plot <- draw.triple.venn(\n");
  fprintf(f, "area1 = %d,\n", P->nA);
  fprintf(f, "area2 = %d,\n", P->nB);
  fprintf(f, "area3 = %d,\n", P->nC);
  fprintf(f, "n12 = %d,\n", P->targets[3]);
  fprintf(f, "n13 = %d,\n", P->targets[5]);
  fprintf(f, "n23 = %d,\n", P->targets[6]);
  fprintf(f, "n123 = %d,\n", P->targets[7]);
  fprintf(f, "category = c(\"A\", \"B\", \"C\"),\n");
  fprintf(f, "fill = c(\"orange\", \"red\", \"green\"),\n");
  fprintf(f, "lty = \"solid\",\n");
  fprintf(f, "cex = 2,\n");
  fprintf(f, "cat.cex = 2,\n");
  fprintf(f, "cat.col = c(\"orange\", \"red\", \"green\")\n");
  fprintf(f, ");\n\n");
}

void writeRscript4(FILE *f) {
  fprintf(f, "venn.plot <- draw.quad.venn(\n");
  fprintf(f, "area1 = %d,\n", P->nA);
  fprintf(f, "area2 = %d,\n", P->nB);
  fprintf(f, "area3 = %d,\n", P->nC);
  fprintf(f, "area4 = %d,\n", P->nD);
  fprintf(f, "n12 = %d,\n", P->targets[3]);
  fprintf(f, "n13 = %d,\n", P->targets[5]);
  fprintf(f, "n14 = %d,\n", P->targets[9]);
  fprintf(f, "n23 = %d,\n", P->targets[6]);
  fprintf(f, "n24 = %d,\n", P->targets[10]);
  fprintf(f, "n34 = %d,\n", P->targets[12]);
  fprintf(f, "n123 = %d,\n", P->targets[7]);
  fprintf(f, "n124 = %d,\n", P->targets[11]);
  fprintf(f, "n134 = %d,\n", P->targets[13]);
  fprintf(f, "n234 = %d,\n", P->targets[14]);
  fprintf(f, "n1234 = %d,\n", P->targets[15]);
  fprintf(f, "category = c(\"A\", \"B\", \"C\", \"D\"),\n");
  fprintf(f, "fill = c(\"orange\", \"red\", \"green\", \"blue\"),\n");
  fprintf(f, "lty = \"solid\",\n");
  fprintf(f, "cex = 2,\n");
  fprintf(f, "cat.cex = 2,\n");
  fprintf(f, "cat.col = c(\"orange\", \"red\", \"green\", \"blue\")\n");
  fprintf(f, ");\n\n");
}

void writeRscript() {
  int n = 0;
  FILE *f;
  n = num_sets();

  f = fopen(P->Rfile, "w");
  fprintf(f, "# Generated by soi - (c) 2016, A. Riva\n");
  fprintf(f, "library(VennDiagram)\n\n");
  
  switch (n) {
  case 2:
    writeRscript2(f);
    break;
  case 3:
    writeRscript3(f);
    break;
  case 4:
    writeRscript4(f);
    break;
  }

  fprintf(f, "# Write image to file\n");
  fprintf(f, "png(filename=\"%s\");\n", P->imgfile);
  fprintf(f, "grid.draw(venn.plot);\n");
  fprintf(f, "dev.off()\n");
}

/* Main */

int main(int argc, char *argv[]) {
  int it, i;
  char *key;
  float p;

  srand(time(0));
  initialize(argc, argv);
  /* dump(); */

  for (it = 0; it < P->niters; it++) {
    countIntersecting();
    for (i = 0; i < MAXCNT; i++) {
      if (P->counts[i] >= P->targets[i]) {
	P->successes[i]++;
      }
    }
  }
  /* dump(); */

  printf("# N = %d\n", P->ndata);
  for (i = 0; i < MAXCNT; i++) {
    if (P->targets[i] > 0) {
      key = idxToKey(i);
      if (strlen(key) > 1) {
	p = 1.0*(1 + P->successes[i])/(1 + P->niters);
	/* printf("%s: Expected=%d, Observed=%d, Successes=%d/%d, P=%f\n", key, (int)(P->ndata*P->expected[i]), P->targets[i], P->successes[i], P->niters, (1.0*P->successes[i])/P->niters); */
	switch (P->outfmt) {
	case FMTDEF:
	  printf("%s: Expected=%d, Observed=%d, P=%f\n", key, (int)(P->ndata*P->expected[i]), P->targets[i], p);
	  break;
	case FMTTAB:
	  printf("%s\t%d\t%d\t%f\n", key, (int)(P->ndata*P->expected[i]), P->targets[i], p);
	  break;
	}
      }
    }
  }

  if (strlen(P->Rfile) > 0) {
    writeRscript();
  }
}

