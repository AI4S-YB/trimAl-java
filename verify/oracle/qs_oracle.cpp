#include <cstdio>

// Exact copy of trimAl v1.5.1 utils::quicksort(int**) (utils.cpp:338-374).
void quicksort(int **vect, int ini, int fin) {
  float elem_div; int i, j;
  if ((ini >= fin) || (fin < 0)) return;
  elem_div = vect[fin][0]; i = ini - 1; j = fin;
  while (true) {
    while (vect[++i][0] < elem_div) if (i == fin) break;
    while (vect[--j][0] > elem_div) if (j == 0) break;
    if (i < j) { int *t = vect[i]; vect[i] = vect[j]; vect[j] = t; } else break;
  }
  { int *t = vect[i]; vect[i] = vect[fin]; vect[fin] = t; }
  quicksort(vect, ini, i - 1); quicksort(vect, i + 1, fin);
}

int main() {
  // adversarial case 1: many equal keys (order of ties is exchange-history dependent)
  int lens[24] = {5,3,5,5,1,5,3,3,5,1,1,5,3,5,1,3,5,5,3,1,5,3,1,5};
  int* seqs[24];
  for (int i = 0; i < 24; i++) { seqs[i] = new int[2]; seqs[i][0] = lens[i]; seqs[i][1] = i; }
  quicksort(seqs, 0, 23);
  for (int i = 0; i < 24; i++) printf("%d ", seqs[i][1]);
  printf("\n");
  // case 2: reverse sorted
  for (int i = 0; i < 24; i++) { seqs[i][0] = 24 - i; seqs[i][1] = i; }
  quicksort(seqs, 0, 23);
  for (int i = 0; i < 24; i++) printf("%d ", seqs[i][1]);
  printf("\n");
  // case 3: all equal
  for (int i = 0; i < 24; i++) { seqs[i][0] = 7; seqs[i][1] = i; }
  quicksort(seqs, 0, 23);
  for (int i = 0; i < 24; i++) printf("%d ", seqs[i][1]);
  printf("\n");
  return 0;
}
