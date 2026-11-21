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
import trimal.util.Utils;

/**
 * Faithful transcription of trimAl v1.5.1 statisticsGaps.cpp. The alignment matrix is
 * byte[sequence][column] with ISO-8859-1 bytes (C++ used std::string rows).
 */
public final class StatisticsGaps {

  int columnLength; // species
  int columns;      // aminos
  int maxGaps;
  int halfWindow;
  int dataType;
  int[] gapsInColumn;
  int[] aminosXInColumn;
  int[] gapsWindow;
  int[] numColumnsWithGaps;

  public StatisticsGaps(String[] alignmentMatrix, int species, int aminos, int dataType_) throws IOException {
    columnLength = species;
    columns = aminos;
    maxGaps = 0;
    halfWindow = 0;
    dataType = dataType_;

    char indet = (dataType == Utils.AAType) ? 'X' : 'N';

    gapsInColumn = new int[columns];
    aminosXInColumn = new int[aminos];
    gapsWindow = new int[columns];
    numColumnsWithGaps = new int[species + 1];

    // Column counts are independent integer sums over each column, so running
    // columns on daemon workers is bit-identical to the serial loop; the shared
    // histogram and maxGaps are then rebuilt in one cheap serial pass.
    if (Parallel.threads() > 1 && columns >= 2 * Parallel.threads()
        && (long) columns * columnLength > 20_000_000L) {
      Parallel.forIndex(0, columns, i -> {
        int gaps = 0, indetCount = 0;
        for (int j = 0; j < columnLength; j++) {
          char c = Utils.charAt(alignmentMatrix[j], i);
          if (c == '-')
            gaps++;
          else if (c == indet)
            indetCount++;
        }
        gapsInColumn[i] = gaps;
        aminosXInColumn[i] = indetCount;
      });
    } else {
      for (int i = 0; i < columns; i++) {
        for (int j = 0; j < columnLength; j++) {
          char c = Utils.charAt(alignmentMatrix[j], i);
          if (c == '-')
            gapsInColumn[i]++;
          else if (c == indet)
            aminosXInColumn[i]++;
        }
      }
    }

    for (int i = 0; i < columns; i++) {
      numColumnsWithGaps[gapsInColumn[i]]++;
      gapsWindow[i] = gapsInColumn[i];
      if (gapsWindow[i] > maxGaps) maxGaps = gapsWindow[i];
    }
  }

  public boolean applyWindow(int _halfWindow) {
    if ((halfWindow == _halfWindow) || (_halfWindow > columns / 4))
      return false;

    Utils.initlVect(numColumnsWithGaps, columnLength + 1, 0);
    Utils.initlVect(gapsWindow, columns, 0);

    maxGaps = 0;
    halfWindow = _halfWindow;
    int window = (2 * halfWindow + 1);

    for (int i = 0; i < columns; i++) {
      gapsWindow[i] = 0;
      for (int j = i - halfWindow; j <= i + halfWindow; j++) {
        if (j < 0)
          gapsWindow[i] += gapsInColumn[-j];
        else if (j >= columns)
          gapsWindow[i] += gapsInColumn[((2 * columns - j) - 2)];
        else
          gapsWindow[i] += gapsInColumn[j];
      }

      gapsWindow[i] = Utils.roundInt((double) gapsWindow[i] / window);
      numColumnsWithGaps[gapsWindow[i]]++;

      if (gapsWindow[i] > maxGaps)
        maxGaps = gapsWindow[i];
    }
    return true;
  }

  public int[] getGapsWindow() {
    return gapsWindow;
  }

  public double calcCutPoint(float minInputAlignment, float gapThreshold) {
    double cuttingPoint_MinimumConserv, cuttingPoint_gapThreshold;
    int i, acum;

    cuttingPoint_gapThreshold = (double) columnLength * gapThreshold;

    cuttingPoint_MinimumConserv = Utils.roundInt((double) (columns * minInputAlignment) / 100.0);
    if (cuttingPoint_MinimumConserv > columns)
      cuttingPoint_MinimumConserv = columns;

    for (i = 0, acum = 0; i < columnLength; i++) {
      acum += numColumnsWithGaps[i];
      if (acum >= cuttingPoint_MinimumConserv)
        break;
    }

    if (numColumnsWithGaps[i] != 0)
      cuttingPoint_MinimumConserv =
          (double) (i - (float) (acum - cuttingPoint_MinimumConserv) / numColumnsWithGaps[i]);
    else
      cuttingPoint_MinimumConserv = 0;

    return Utils.max(cuttingPoint_MinimumConserv, cuttingPoint_gapThreshold);
  }

  public int calcCutPoint2ndSlope() {
    // Upstream scans beyond an empty histogram when columns == 0 (C++ UB,
    // sometimes returning arbitrary output rather than crashing). There is no
    // meaningful slope without columns; avoid exposing a Java OOB stacktrace.
    if (columns == 0) return 0;
    float maxSlope = -1;
    int prev, pprev, maxIter, act = 0, max = 0;

    float[] secondSlopeVector = new float[maxGaps + 1];
    Utils.initlVect(secondSlopeVector, maxGaps, -1.0f);
    maxIter = maxGaps + 1;

    for (act = 0, max = 0; numColumnsWithGaps[act] == 0; act++)
      max = act + 1;

    act = 0;
    while (act < maxIter) {

      while ((numColumnsWithGaps[act]) == 0)
        act++;
      pprev = act;
      if ((act + 1) >= maxIter)
        break;

      do {
        act++;
      } while ((numColumnsWithGaps[act]) == 0);
      prev = act;
      if ((act + 1) >= maxIter)
        break;

      do {
        act++;
      } while ((numColumnsWithGaps[act]) == 0);
      if (act >= maxIter)
        break;

      secondSlopeVector[act] = ((float) (act - pprev) / columnLength);
      secondSlopeVector[act] /= ((float) (numColumnsWithGaps[act] + numColumnsWithGaps[prev]) / columns);

      if (secondSlopeVector[pprev] != -1.0) {
        if ((secondSlopeVector[act] / secondSlopeVector[pprev]) > maxSlope) {
          maxSlope = (secondSlopeVector[act] / secondSlopeVector[pprev]);
          max = pprev;
        }
      } else if (secondSlopeVector[prev] != -1.0) {
        if ((secondSlopeVector[act] / secondSlopeVector[prev]) > maxSlope) {
          maxSlope = (secondSlopeVector[act] / secondSlopeVector[prev]);
          max = pprev;
        }
      }
      act = prev;
    }

    return max;
  }

  public void printGapsColumns(CppOut cout) throws IOException {
    int[] vectAux = new int[columns];

    if (halfWindow == 0)
      Utils.copyVect(gapsInColumn, vectAux, columns);
    else
      Utils.copyVect(gapsWindow, vectAux, columns);

    cout.p("| Residue\t  % Gaps \t   Gap Score   |").endl();
    cout.p("| Number \t         \t               |").endl();
    cout.p("+----------------------------------------------+").endl();
    cout.precision(10);

    for (int i = 0; i < columns; i++)
      cout.p("  ").setw(5).p(i).p("\t\t").setw(10).p((vectAux[i] * 100.0) / columnLength)
          .p("\t").setw(7).p(1 - ((vectAux[i] * 1.0) / columnLength)).endl();
  }

  public void printGapsAcl(CppOut cout) throws IOException {
    int acm, i;

    cout.p("| Number of\t        \t|\t Cumulative \t% Cumulative\t|\tNumber of Gaps\t  % Gaps  \tGap Score  |").endl();
    cout.p("| Residues \t% Length\t|\tNumberResid.\t   Length   \t|\t  per Column  \tper Column\tper Column |").endl();
    cout.p("+-------------------------------+-----------------------------")
        .p("----------+--------------------------------------------------+").endl();
    cout.precision(10);

    for (i = 0, acm = 0; i <= maxGaps; i++) {
      if (numColumnsWithGaps[i] != 0) {
        acm += numColumnsWithGaps[i];
        cout.p("  ").left().p(numColumnsWithGaps[i]).p("\t\t").setw(10).p((numColumnsWithGaps[i] * 100.0) / columns)
            .p("\t\t").p(acm).p("\t\t").setw(10).p((acm * 100.0) / columns)
            .p("\t\t").p(i).p("\t\t").setw(10).p((i * 1.0) / columnLength).p("\t").setw(10).p(1 - ((i * 1.0) / columnLength)).endl();
      }
    }
  }
}
