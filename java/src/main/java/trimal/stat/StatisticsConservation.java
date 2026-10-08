package trimal.stat;

import java.io.IOException;

import trimal.util.CppOut;
import trimal.util.Parallel;
import trimal.util.Utils;

/**
 * Transcription of trimAl v1.5.1 statisticsConservation.cpp. Float/double
 * promotion and accumulation order are checked against the real C++ modules.
 * C++ exp(float) is emulated explicitly rather than through the double overload.
 *
 * The two O(sequences^2 x columns) loops (matrixIdentity construction and the
 * per-column Q/MDK vectors) are parallelized across rows/columns. Units are
 * independent and write only distinct slots while keeping their serial float
 * accumulation order, so results are bit-identical to the serial loop; small
 * workloads and -Dtrimal.threads=1 run the plain serial code.
 */
public final class StatisticsConservation {

  int columns;
  int sequences;
  int dataType;
  int halfWindow;
  float[] Q;
  float[] MDK;
  float[] MDK_Window;
  float[] matrixIdentity; // row-major: identityAt(j, k) == matrixIdentity[j * sequences + k]
  SimilarityMatrix simMatrix;

  public StatisticsConservation(String[] alignmentMatrix, int species, int aminos, int dataType_) throws IOException {
    columns = aminos;
    sequences = species;
    dataType = dataType_;
    halfWindow = -1;

    Q = new float[columns];
    MDK = new float[columns];
    MDK_Window = new float[columns];
    matrixIdentity = new float[sequences * sequences];

    simMatrix = null;

    calculateMatrixIdentity(alignmentMatrix);
  }

  private void calculateMatrixIdentity(String[] alignmentMatrix) throws IOException {
    char indet = (dataType == Utils.AAType) ? 'X' : 'N';

    long pairs = (long) sequences * (sequences - 1) / 2;
    if (Parallel.threads() > 1 && sequences >= 2 * Parallel.threads()
        && pairs * columns > 50_000_000L) {
      Parallel.forIndex(0, sequences, i -> computeIdentityRow(alignmentMatrix, indet, i));
    } else {
      for (int i = 0; i < sequences; i++)
        computeIdentityRow(alignmentMatrix, indet, i);
    }
  }

  /** Owns exactly the (i, j>i) pairs: no other unit reads or writes these slots. */
  private void computeIdentityRow(String[] alignmentMatrix, char indet, int i) {
    for (int j = i + 1; j < sequences; j++) {
      int sum = 0, length = 0;
      for (int k = 0; k < columns; k++) {
        char ci = alignmentMatrix[i].charAt(k);
        char cj = alignmentMatrix[j].charAt(k);
        if ((ci != '-') && (ci != indet)) {
          if ((cj != '-') && (cj != indet))
            if (cj == ci)
              sum++;
          length++;
        } else if ((cj != '-') && (cj != indet))
          length++;
      }
      float value = (float) (100.0 - ((float) sum / length) * 100.0);
      matrixIdentity[j * sequences + i] = value;
      matrixIdentity[i * sequences + j] = value;
    }
  }

  public boolean calculateVectors(String[] alignmentMatrix, int[] gaps) throws IOException {
    char indet = (dataType == Utils.AAType) ? 'X' : 'N';

    if (simMatrix == null)
      return false;

    // getDistance reports undefined symbols on stderr in scan order, so keep the
    // serial loop whenever any residue could trigger that message.
    boolean parallel = Parallel.threads() > 1
        && columns >= 2 * Parallel.threads()
        && (long) columns * sequences * (sequences - 1) / 2 > 20_000_000L
        && !hasUndefinedSymbol(alignmentMatrix, indet);

    if (parallel)
      Parallel.forIndex(0, columns, i -> computeColumn(alignmentMatrix, gaps, indet, i));
    else
      for (int i = 0; i < columns; i++)
        computeColumn(alignmentMatrix, gaps, indet, i);

    return true;
  }

  /** Owns exactly the output slots Q[i]/MDK[i] for its column. */
  private void computeColumn(String[] alignmentMatrix, int[] gaps, char indet, int i) throws IOException {
    float num = 0, den = 0;
    for (int j = 0; j < sequences; j++) {
      char cj = alignmentMatrix[j].charAt(i);
      if ((cj != '-') && (cj != indet))
        for (int k = j + 1; k < sequences; k++) {
          char ck = alignmentMatrix[k].charAt(i);
          if ((ck != '-') && (ck != indet)) {
            num += matrixIdentity[j * sequences + k] * simMatrix.getDistance(cj, ck);
            den += matrixIdentity[j * sequences + k];
          }
        }
    }
    Q[i] = (den == 0) ? 0 : num / den;
    MDK[i] = Utils.expf(-Q[i]);

    if (gaps != null)
      if (((float) gaps[i] / sequences) >= 0.8) MDK[i] = 0;

    if (MDK[i] > 1) MDK[i] = 1;
  }

  /** Over-approximates whether any residue could make getDistance print (conservative). */
  private boolean hasUndefinedSymbol(String[] alignmentMatrix, char indet) {
    for (int j = 0; j < sequences; j++) {
      String row = alignmentMatrix[j];
      for (int k = 0; k < columns; k++) {
        char c = row.charAt(k);
        if ((c != '-') && (c != indet) && !simMatrix.isDefined(c))
          return true;
      }
    }
    return false;
  }

  public boolean applyWindow(int _halfWindow) {
    if ((halfWindow == _halfWindow) || (_halfWindow > columns / 4))
      return false;

    halfWindow = _halfWindow;
    int window = 2 * halfWindow + 1;

    for (int i = 0; i < columns; i++) {
      for (int j = i - halfWindow; j <= i + halfWindow; j++) {
        if (j < 0) MDK_Window[i] += MDK[-j];
        else if (j >= columns) MDK_Window[i] += MDK[((2 * columns - j) - 2)];
        else MDK_Window[i] += MDK[j];
      }
      MDK_Window[i] = MDK_Window[i] / (float) window;
    }
    return true;
  }

  public boolean isDefinedWindow() {
    return (halfWindow != -1);
  }

  public float[] getMdkwVector() {
    return MDK_Window;
  }

  public float[] getMdkVector() {
    return MDK;
  }

  public boolean setSimilarityMatrix(SimilarityMatrix sm) {
    if (sm == null)
      return false;
    simMatrix = sm;
    return true;
  }

  public boolean isSimMatrixDef() {
    return (simMatrix != null);
  }

  public double calcCutPoint(float minInputAlignment, float similThreshold) {
    double cuttingPoint_MinimumConserv, cuttingPoint_SimilThreshold;
    int i, highestPos;

    float[] vectAux = new float[columns];

    Utils.copyVect(MDK_Window, vectAux, columns);
    Utils.quicksort(vectAux, 0, columns - 1);

    for (i = columns - 1; i >= 0; i--)
      if (vectAux[i] < similThreshold)
        break;
    // With no value below the threshold C++ reads vectAux[-1], which is UB.
    // Zero is a compatibility fallback, not a portable or guaranteed C++ value.
    // verify/numeric records ordinary-build output separately from sanitizer UB.
    cuttingPoint_SimilThreshold = (i >= 0) ? vectAux[i] : 0.0f;

    highestPos = (int) ((double) (columns - 1) * (100.0 - minInputAlignment) / 100.0);
    highestPos = highestPos < (columns - 1) ? highestPos : columns - 1;
    cuttingPoint_MinimumConserv = vectAux[highestPos];

    return (cuttingPoint_MinimumConserv < cuttingPoint_SimilThreshold
        ? cuttingPoint_MinimumConserv
        : cuttingPoint_SimilThreshold);
  }

  public void printConservationColumns(CppOut cout) throws IOException {
    cout.p("| Residue\t Similarity  |").endl();
    cout.p("| Number \t    Value    |").endl();
    cout.p("+----------------------------+").endl();
    cout.precision(10);

    if (MDK_Window != null)
      for (int i = 0; i < columns; i++)
        cout.p("  ").setw(5).p(i).p("\t\t").setw(7).p(MDK_Window[i]).endl();
    else
      for (int i = 0; i < columns; i++)
        cout.p("  ").setw(5).p(i).p("\t\t").setw(7).p(MDK[i]).endl();
  }

  public void printConservationAcl(CppOut cout) throws IOException {
    float[] vectAux = new float[columns];

    if (MDK_Window != null) Utils.copyVect(MDK_Window, vectAux, columns);
    else Utils.copyVect(MDK, vectAux, columns);

    Utils.quicksort(vectAux, 0, columns - 1);

    cout.p("| Number of\t        \t|\t Cumulative \t% Cumulative\t|   Similarity   |").endl();
    cout.p("| Residues \t% Length\t|\tNumberResid.\t   Length   \t|     Value      |").endl();
    cout.p("+-------------------------------+---------------------------------------+----------------+").endl();
    cout.precision(10);

    float refer = vectAux[columns - 1];
    int acm = 0, num = 1;

    for (int i = columns - 2; i >= 0; i--) {
      acm++;
      if (refer != vectAux[i]) {
        cout.p("  ").p(num).p("\t\t").setw(10).p((float) num / columns * 100.0)
            .p("\t\t").p(acm).p("\t\t").setw(10).p((float) acm / columns * 100.0).p("\t")
            .setw(15).p(refer).endl();
        refer = vectAux[i];
        num = 1;
      } else num++;
    }
    acm++;
    cout.p("  ").p(num).p("\t\t").setw(10).p((float) num / columns * 100.0)
        .p("\t\t").p(acm).p("\t\t").setw(10).p((float) acm / columns * 100.0).p("\t")
        .setw(15).p(refer).endl();
  }
}
