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

package trimal.stat;

import java.io.IOException;

import trimal.util.CppOut;
import trimal.util.Parallel;

/**
 * Faithful transcription of trimAl v1.5.1 sequencesMatrix.cpp — per-sequence residue
 * numbering (nth non-gap residue, gaps = 0) used by comparison and backtranslation.
 */
public final class SequencesMatrix {

  int resNumber;
  int seqsNumber;
  int[][] matrix;
  String[] seqsName;

  public SequencesMatrix(String[] sequences, String[] names, int sequenNumber, int residNumber) throws IOException {
    int i, j, k;

    resNumber = residNumber;
    seqsNumber = sequenNumber;

    seqsName = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++) seqsName[i] = names[i];

    matrix = new int[sequenNumber][];
    // Each unit fills only its own row (per-residue numbering is independent), so
    // running rows on daemon workers is bit-identical to the serial loop.
    if (Parallel.threads() > 1 && sequenNumber >= 2 * Parallel.threads()
        && (long) sequenNumber * residNumber > 20_000_000L) {
      Parallel.forIndex(0, sequenNumber, row -> {
        int[] values = new int[residNumber];
        int counter = 1;
        for (int jj = 0; jj < residNumber; jj++) {
          if (sequences[row].charAt(jj) != '-') {
            values[jj] = counter;
            counter++;
          } else {
            values[jj] = 0;
          }
        }
        matrix[row] = values;
      });
    } else {
      for (i = 0; i < sequenNumber; i++) {
        matrix[i] = new int[residNumber];
        k = 1;
        for (j = 0; j < residNumber; j++) {
          if (sequences[i].charAt(j) != '-') {
            matrix[i][j] = k;
            k++;
          } else {
            matrix[i][j] = 0;
          }
        }
      }
    }
  }

  public int getResidNumber() {
    return resNumber;
  }

  public void printMatrix(CppOut cout) throws java.io.IOException {
    int i, j, k;
    for (i = 0; i < resNumber; i += 20) {
      for (j = 0; j < seqsNumber; j++) {
        for (k = i; k < (20 + i) && k < resNumber; k++) {
          cout.setw(4).p(matrix[j][k]).p(" ");
        }
        cout.endl();
      }
      cout.endl();
    }
  }

  /** getColumn(column, vector): whole column or 0 padding for out-of-range columns. */
  public void getColumn(int column, int[] vector) {
    int i;
    if (column < resNumber)
      for (i = 0; i < seqsNumber; i++)
        vector[i] = matrix[i][column];
    else
      for (i = 0; i < seqsNumber; i++)
        vector[i] = 0;
  }

  /** getColumn(value, row, vector): find value in row, return that column; -1 when absent. */
  public void getColumn(int value, int row, int[] vector) {
    int i, j;
    for (i = 0; i < resNumber; i++)
      if (matrix[row][i] == value) break;

    if (i < resNumber)
      for (j = 0; j < seqsNumber; j++)
        vector[j] = matrix[j][i];
    else
      for (j = 0; j < seqsNumber; j++)
        vector[j] = -1;
  }

  public void setOrder(int[] order) {
    // sequencesMatrix.cpp:162-178 copies residue rows only. Names deliberately
    // retain their original order, even though getSequence() then uses that order.
    int[][] matrixAux = new int[seqsNumber][resNumber];
    for (int i = 0; i < seqsNumber; i++)
      System.arraycopy(matrix[order[i]], 0, matrixAux[i], 0, resNumber);
    for (int i = 0; i < seqsNumber; i++)
      System.arraycopy(matrixAux[i], 0, matrix[i], 0, resNumber);
  }

  /** getSequence(name, vector): first row matching the name; false when absent. */
  public boolean getSequence(String sequence, int[] vector) {
    int i, j, seqPos = -1;
    for (i = 0; i < seqsNumber; i++) {
      if (seqsName[i].equals(sequence)) {
        seqPos = i;
        break;
      }
    }
    if (seqPos == -1) return false;
    for (j = 0; j < resNumber; j++) vector[j] = matrix[seqPos][j];
    return true;
  }
}
