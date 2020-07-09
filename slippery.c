/* slippery.c by Kaeding 2018-2019 */
/* slippery hill-climbing method for ciphertext-only attack on
   periodic polyalphabetic substution ciphers */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "monograms.h"  /* single-letter frequencies */
#include "tetragrams.h" /* tetragram frequencies */

#define MAXTEXTLEN    10000
#define MAXKEYLEN     100
#define IOCTHRESHOLD  1.65
#define IOCMULTTHRESH 1.2

char alphabet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

/* the index of coincidence (IoC) for a text */
double index_of_coincidence(char* text) {
  int counts[26],total=0,i,length,numer=0;

  for (i=0;i<26;i++)
    counts[i] = 0;
  length = strlen(text);
  for (i=0;i<length;i++)
    counts[text[i]-'A']++;
  for (i=0;i<26;i++) {
    numer += counts[i]*(counts[i]-1);
    total += counts[i];
    }
  return (26.*numer)/(total*(total-1));
  }

/* the fitness of a text, based on tetragram frequencies */
double fitness(char *text) {
  int length,i,count=0;
  double result=0.;

  length = strlen(text);
  for (i=0;i<length-3;i++) {
    result += tetragrams[(text[i+0]-'A')*26*26*26
                        +(text[i+1]-'A')*26*26
                        +(text[i+2]-'A')*26
                        +(text[i+3]-'A')];
    count++;
    }
  return result/count;
  }

/* the position of a character in a string */
int position(char c, char* s) {
  int i,length;

  length = strlen(s);
  for (i=0;i<length;i++)
    if (c == s[i])
      return i;
  return -1;
  }

/* the decryption function for the polyalphabetic cipher */
void decrypt(char* c, char* p, char s[MAXKEYLEN][26], int keylen) {
  int length,i;

  length = strlen(c);
  p[length] = '\0';
  for (i=0;i<length;i++)
    p[i] = alphabet[position(c[i],s[i%keylen])];
  return;
  }

/* find the single-letter frequencies of a text */
void monogram_frequencies(char* text, double* freqs) {
  int i,total=0,length;

  for (i=0;i<26;i++)
    freqs[i] = 0;
  length = strlen(text);
  for (i=0;i<length;i++) {
    freqs[text[i]-'A']++;
    total++;
    }
  for (i=0;i<26;i++)
    freqs[i] *= (1./total);
  return;
  }

/* swap two random characters in an alphabet */
void random_swap(char s[26]) {
  int i,j;
  char temp;

  i = j = random()%26;
  while (i == j)
    j = random()%26;
  temp = s[i];
  s[i] = s[j];
  s[j] = temp;
  return;
  }

/* randomize an alphabet */
void randomize(char s[26]) {
  int i,j;

  for (i=0;i<26;i++)
    s[i] = -1;
  for (i=0;i<26;i++) {
    j = random() % 26;
    while (s[j] != -1)
      j = random() % 26;
    s[j] = alphabet[i];
    }
  return;
  }

/* copy the set of key alphabets */
void copy_keys(char source[MAXKEYLEN][26],
               char target[MAXKEYLEN][26],int keylen) {
  int i,j;

  for (i=0;i<keylen;i++)
    for (j=0;j<26;j++)
      target[i][j] = source[i][j];
  return;
  }

int main(int argc, char** argv) {
  char c[MAXTEXTLEN];          /* the ciphertext                 */
  char p[MAXTEXTLEN];          /* the plaintext                  */
  char slice[MAXTEXTLEN];      /* one slice of ciphertext        */
  char bestp[MAXTEXTLEN];      /* best plaintext so far          */
  char pk[MAXKEYLEN][26];      /* the parent key alphabets       */
  char ck[MAXKEYLEN][26];      /* the child key alphabets        */
  char bestk[MAXKEYLEN][26];   /* best keys so far               */
  double fitp;                 /* fitness of parent              */
  double fitc;                 /* fitness of child               */
  double bestf;                /* best fitness so far            */
  double ioc;                  /* index of coincidence           */
  double oldioc=99.;           /* previous IoC                   */
  long int count,bigcount=0;   /* counters for loops             */
  int period=0;                /* period of the cipher           */
  int i,j,k;
  int length;                  /* length of the ciphertext       */
  int found=0;                 /* have we found the period?      */
  int indexf,indexr;           /* used in setting initial keys   */
  double maxf,maxr;            /* used in setting initial keys   */
  double reference[26];        /* holds reference monogram freqs */
  double freqs[MAXKEYLEN][26]; /* monogram frequencies of slice  */

  strcpy(c,argv[1]);
  length = strlen(c);

  srandom(time(0));

  /* find the period and cut the ciphertext into slices */
  while (!found) {
    period++;
    ioc = 0.;
    for (i=0;i<period;i++) {
      for (j=0;j<length/period;j++) {
        slice[j] = c[period*j+i];
        }
      slice[j] = '\0';
      ioc += index_of_coincidence(slice);
      }
    ioc /= period;
    if ((ioc > IOCTHRESHOLD) && (ioc > IOCMULTTHRESH*oldioc))
      found = 1;
    oldioc = ioc;
    }

  /* uncomment the next line if you want to tell the program
     what the period is as a parameter on the command line */
  /* period = atoi(argv[2]); */

  /* set the initial key alphabets */
  for (i=0;i<period;i++) {
    monogram_frequencies(slice,freqs[i]);
    for (j=0;j<26;j++) {
      pk[i][j] = -1.;
      reference[j] = monograms[j];
      }
    for (j=0;j<26;j++) {
      maxf = maxr = -1.;
      indexf = indexr = 0;
      for (k=0;k<26;k++) {
        if (freqs[i][k] > maxf) {
          indexf = k;
          maxf = freqs[i][k];
          }
        if (reference[k] > maxr) {
          indexr = k;
          maxr = reference[k];
          }
        }
      pk[i][indexr] = alphabet[indexf];
      freqs[i][indexf] = -1.;
      reference[indexr] = -1;
      }
    }

  decrypt (c,p,pk,period);
  bestf = fitness(p);

  while (bigcount < 5000000*period*period/length) /* main loop */
    for (j=0;j<period;j++) {
      randomize(pk[j]);
      decrypt(c,p,pk,period);
      fitp = fitness(p);
      count = 0;
      while (count < 1000) { /* inner loop */
        copy_keys(pk,ck,period);
        random_swap(ck[j]);
        decrypt(c,p,ck,period);
        fitc = fitness(p);
        if (fitc > fitp) {
          copy_keys(ck,pk,period);
          fitp = fitc;
          count = 0;
          }
        else
          count++;
        if (fitc > bestf) {
          copy_keys(ck,bestk,period);
          bestf = fitc;
          bigcount = 0;
          strcpy(bestp,p);
          }
        else
          bigcount++;
        }
      }

  /* print the results */
  printf("%s\n",bestp);
  printf("key alphabets:\n");
  for (i=0;i<period;i++) {
    printf("    [");
    for (j=0;j<26;j++)
      printf("%c",bestk[i][j]);
    printf("]\n");
    }
  printf("fitness: %8.4f\n",bestf);

  return 0;
  }
