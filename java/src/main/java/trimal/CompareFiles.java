package trimal;

import java.io.IOException;

import trimal.util.CppOut;
import trimal.util.Utils;

/**
 * Faithful transcription of trimAl v1.5.1 compareFiles.cpp.
 */
public final class CompareFiles {

  private CompareFiles() {}

  public static int algorithm(Alignment[] vectAlignments, String[] fileNames,
      float[] columnsValue, int numAlignments, boolean verbosity) throws IOException {

    int[] numResiduesAlig, correspNames, columnSeqMatrix, columnSeqMatrixAux;
    int i, j, k, l, m, numSeqs, pairRes, hits, alig = 0;
    float max = 0, value = 0;
    float[][] vectHits = new float[numAlignments][];
    boolean appearErrors = false;
    String[] names;

    numSeqs = vectAlignments[0].getNumSpecies();

    names = new String[numSeqs];
    correspNames = new int[numSeqs];
    numResiduesAlig = new int[numAlignments];
    columnSeqMatrix = new int[numSeqs];
    columnSeqMatrixAux = new int[numSeqs];

    for (i = 1; i < numAlignments; i++) {
      if (numSeqs != vectAlignments[i].getNumSpecies()) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The files to compare do not have ")
              .p("the same number of sequences").endl().endl();
        }
        appearErrors = true;
        break;
      }

      vectAlignments[i].getSequencesNames(names);
      if (!vectAlignments[0].getSeqNameOrder(names, correspNames)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The files to compare do not")
              .p(" have the sequence names").endl().endl();
        }
        appearErrors = true;
        break;
      }
    }

    for (i = 1; ((i < numAlignments) && (!appearErrors)); i++) {
      vectAlignments[i].getSequencesNames(names);
      vectAlignments[0].getSeqNameOrder(names, correspNames);
      vectAlignments[i].setSeqMatrixOrder(correspNames);
    }

    for (i = 0; ((i < numAlignments) && (!appearErrors)); i++)
      numResiduesAlig[i] = vectAlignments[i].getNumAminos();

    for (i = 0; ((i < numAlignments) && (!appearErrors)); i++, value = 0) {

      if (verbosity) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.endl();
        }
      }

      vectHits[i] = new float[numResiduesAlig[i]];
      Utils.initlVect(vectHits[i], numResiduesAlig[i], 0);

      for (j = 0, pairRes = 0, hits = 0; j < numResiduesAlig[i]; j++, pairRes = 0, hits = 0) {

        vectAlignments[i].getColumnSeqMatrix(j, columnSeqMatrix);

        for (k = 0; k < numSeqs; k++) {
          if (columnSeqMatrix[k] != 0) {
            for (l = 0; l < i; l++) {
              vectAlignments[l].getColumnSeqMatrix(columnSeqMatrix[k], k, columnSeqMatrixAux);
              for (m = k + 1; m < numSeqs; m++)
                if (columnSeqMatrix[m] != 0) {
                  if (columnSeqMatrix[m] == columnSeqMatrixAux[m])
                    hits++;
                  pairRes++;
                }
            }

            for (l = i + 1; l < numAlignments; l++) {
              vectAlignments[l].getColumnSeqMatrix(columnSeqMatrix[k], k, columnSeqMatrixAux);
              for (m = k + 1; m < numSeqs; m++)
                if (columnSeqMatrix[m] != 0) {
                  if (columnSeqMatrix[m] == columnSeqMatrixAux[m])
                    hits++;
                  pairRes++;
                }
            }
          }
        }

        if (pairRes != 0) {
          vectHits[i][j] += (float) ((1.0 * hits) / pairRes);
          value += vectHits[i][j];
        }
      }

      if (verbosity) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.p("File:\t\t").p(fileNames[i]).endl().p("Values:\t\tSequences: ").p(numSeqs)
              .p("\tResidues: ").p(numResiduesAlig[i]).p("\tPond. Hits: ").setw(8)
              .p(value).p("\t%Consistency: ").p(value / numResiduesAlig[i]).endl();
        }
      }

      if ((value / numResiduesAlig[i]) > max) {
        alig = i;
        max = value / numResiduesAlig[i];
      }
    }

    if ((verbosity) && (!appearErrors)) {
      try (CppOut cout = CppOut.toStdout()) {
        cout.p("\t\t\t\t\t--------------").endl();
        cout.endl().p("File Selected:\t").p(fileNames[alig]).endl().p("Value:\t\t").p(max).endl().endl();
      }
    }

    if ((columnsValue != null) && (!appearErrors)) {
      Utils.initlVect(columnsValue, numResiduesAlig[alig], -1);
      for (i = 0; i < numResiduesAlig[alig]; i++)
        columnsValue[i] = vectHits[alig][i];
    }

    if (appearErrors)
      return -1;
    return alig;
  }

  public static boolean forceComparison(Alignment[] vectAlignments, int numAlignments,
      Alignment selected, float[] columnsValue) throws IOException {

    int[] correspNames, columnSeqMatrix, columnSeqMatrixAux;
    int i, j, k, ll, numResidues, numSeqs, pairRes, hit;
    boolean appearErrors = false;
    String[] names;

    numResidues = selected.getNumAminos();
    numSeqs = selected.getNumSpecies();

    Utils.initlVect(columnsValue, numResidues, 0);

    names = new String[numSeqs];
    correspNames = new int[numSeqs];
    columnSeqMatrix = new int[numSeqs];
    columnSeqMatrixAux = new int[numSeqs];

    for (i = 0; i < numAlignments; i++) {
      if (numSeqs != vectAlignments[i].getNumSpecies()) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The files to compare do not have ")
              .p("the same number of sequences").endl().endl();
        }
        appearErrors = true;
        break;
      }

      vectAlignments[i].getSequencesNames(names);
      if (!selected.getSeqNameOrder(names, correspNames)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The files to compare do not")
              .p(" have the sequence names").endl().endl();
        }
        appearErrors = true;
        break;
      }
    }

    for (i = 0; i < numAlignments; i++) {
      vectAlignments[i].getSequencesNames(names);
      selected.getSeqNameOrder(names, correspNames);
      vectAlignments[i].setSeqMatrixOrder(correspNames);
    }

    for (i = 0, pairRes = 0, hit = 0; ((i < numResidues) && (!appearErrors)); i++, pairRes = 0, hit = 0) {

      Utils.initlVect(columnSeqMatrix, numSeqs, 0);
      selected.getColumnSeqMatrix(i, columnSeqMatrix);

      for (j = 0; j < numSeqs; j++) {
        if (columnSeqMatrix[j] != 0) {
          for (k = 0; k < numAlignments; k++) {
            Utils.initlVect(columnSeqMatrixAux, numSeqs, 0);
            vectAlignments[k].getColumnSeqMatrix(columnSeqMatrix[j], j, columnSeqMatrixAux);

            for (ll = j + 1; ll < numSeqs; ll++)
              if (columnSeqMatrix[ll] != 0) {
                if (columnSeqMatrix[ll] == columnSeqMatrixAux[ll])
                  hit++;
                pairRes++;
              }
          }
        }
      }
      if (pairRes != 0)
        columnsValue[i] += (float) ((1.0 * hit) / pairRes);
    }

    if (appearErrors)
      return false;
    return true;
  }

  public static boolean applyWindow(int columns, int halfWindow, float[] columnsValue) {
    int i, j, window;

    if (halfWindow > columns / 4)
      return false;
    else
      window = (2 * halfWindow + 1);

    float[] vectAux = new float[columns];
    Utils.copyVect(columnsValue, vectAux, columns);

    for (i = 0; i < columns; i++) {
      columnsValue[i] = 0;
      for (j = i - halfWindow; j <= i + halfWindow; j++) {
        if (j < 0)
          columnsValue[i] += vectAux[-j];
        else if (j >= columns)
          columnsValue[i] += vectAux[((2 * columns - j) - 2)];
        else
          columnsValue[i] += vectAux[j];
      }
      columnsValue[i] /= window;
    }

    return true;
  }

  public static void printStatisticsFileColumns(int numAminos, float[] compareVect) throws IOException {
    try (CppOut cout = CppOut.toStdout()) {
      cout.p("| Residue\tConsistency |").endl();
      cout.p("| Number \t   Value    |").endl();
      cout.p("+---------------------------+").endl();
      cout.precision(10);

      for (int i = 0; i < numAminos; i++)
        cout.p("  ").setw(5).p(i + 1).p("\t")
            .p("\t").p(compareVect[i]).endl();
    }
  }

  public static void printStatisticsFileAcl(int numAminos, float[] compareVect) throws IOException {
    float refer;
    float[] vectAux = new float[numAminos];
    int i, num;

    Utils.copyVect(compareVect, vectAux, numAminos);
    Utils.quicksort(vectAux, 0, numAminos - 1);

    try (CppOut cout = CppOut.toStdout()) {
      cout.p("| Number of\t        \t|\t Cumulative \t% ")
          .p("Cumulative\t|   Consistency   |").endl();
      cout.p("| Residues \t% Length\t|\tNumberResid.\t   ")
          .p("Length   \t|      Value      |").endl();
      cout.p("+-------------------------------+------------")
          .p("---------------------------+-----------------+")
          .endl();
      cout.precision(10);

      refer = vectAux[0];
      num = 1;

      for (i = 1; i < numAminos; i++) {
        if (refer != vectAux[i]) {
          cout.p("  ").p(num).p("\t\t").setw(10).p((float) num / numAminos * 100.0)
              .p("\t\t").p(i).p("\t\t").setw(10).p((float) i / numAminos * 100.0)
              .p("\t").setw(15).p(refer).endl();
          refer = vectAux[i];
          num = 1;
        } else
          num++;
      }

      cout.p("  ").p(num).p("\t\t").setw(10).p((float) num / numAminos * 100.0)
          .p("\t\t").p(i).p("\t\t").setw(10).p((float) i / numAminos * 100.0)
          .p("\t").setw(15).p(refer).endl();
    }
  }
}
