/* ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
 ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****

    trimAl v1.5.1 (Java 11 port): a tool for automated trimming of multiple
    sequence alignments, together with its readAl and statal companions.

    This file is part of the trimAl Java port, a derivative work of trimAl
    (https://github.com/inab/trimal).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** *****
***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** ***** */

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
