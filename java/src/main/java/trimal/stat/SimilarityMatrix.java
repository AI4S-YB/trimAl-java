package trimal.stat;

import java.io.IOException;
import trimal.util.CppOut;
import trimal.util.Utils;
import trimal.util.Values;

/**
 * Faithful transcription of trimAl v1.5.1 similarityMatrix.cpp. All arithmetic is
 * float; the Euclidean distance loop keeps the original j/i/k loop order because the
 * k-sum accumulation order is observable in the float result.
 */
public final class SimilarityMatrix {

  private static final String WHITESPACE = " \t\n\r\f\u000B";

  private int numPositions = 0;
  private int[] vhash = null;      // TAMABC entries
  private float[][] simMat = null;
  private float[][] distMat = null;

  public int getNumPositions() { return numPositions; }

  public void memoryAllocation(int nPos) {
    if (numPositions != 0) memoryDeletion();
    numPositions = nPos;
    vhash = new int[Utils.TAMABC];
    simMat = new float[nPos][nPos];
    distMat = new float[nPos][nPos];
  }

  private void memoryDeletion() {
    numPositions = 0;
    vhash = null;
    simMat = null;
    distMat = null;
  }

  /**
   * loadSimMatrix: first line stripped of spaces/tabs defines the symbol table;
   * rows may or may not carry a leading symbol; asymmetric input is symmetrized;
   * then Euclidean distances are derived.
   */
  public boolean loadSimMatrix(String fn) {
    try {
      byte[] raw = java.nio.file.Files.readAllBytes(java.nio.file.Path.of(fn));
      // CRLF and any byte go through ISO-8859-1 like the C++ ifstream
      StringBuilder content = new StringBuilder(raw.length);
      for (byte b : raw) content.append((char) (b & 0xFF));

      // first line = up to '\n' (getline), max LINE_LENGTH-1 chars
      int nl = content.indexOf("\n");
      if (nl < 0) nl = content.length();
      String firstLine = content.substring(0, Math.min(nl, Utils.LINE_LENGTH - 1));
      StringBuilder symbols = new StringBuilder(Utils.removeSpaces(firstLine));
      for (int i = 0; i < symbols.length(); i++) symbols.setCharAt(i, asciiUpper(symbols.charAt(i)));
      String listSym = symbols.toString();
      int nul = listSym.indexOf('\0');
      if (nul >= 0) listSym = listSym.substring(0, nul);
      memoryAllocation(listSym.length());

      for (int i = 0; i < Utils.TAMABC; i++) vhash[i] = -1;

      for (int i = 0; i < numPositions; i++) {
        char c = listSym.charAt(i);
        if (c >= 'a' && c <= 'z') c = Character.toUpperCase(c);
        if (c >= 'A' && c <= 'Z') {
          if (vhash[c - 'A'] != -1) {
            memoryDeletion();
            return false;
          }
          vhash[c - 'A'] = i;
        } else {
          memoryDeletion();
          return false;
        }
      }

      // Numeric operator>> does not tokenize like atof: it can consume only a
      // prefix, set failbit, or leave a previous destination unchanged at EOF.
      MatrixInput input = new MatrixInput(content.toString(), Math.min(nl + 1, content.length()));
      if (nl >= Utils.LINE_LENGTH) input.failed = true; // getline buffer overflow
      boolean firstColumn = true;
      String first = ""; // C++ reuses the same char buffer after a failed extraction
      for (int i = 0; i < numPositions; i++) {
        int j = 0;
        first = input.word(first);
        char rowSym = first.isEmpty() ? '\0' : first.charAt(0);
        if (firstColumn) {
          rowSym = asciiUpper(rowSym);
          if (!first.isEmpty()) first = rowSym + first.substring(1);
          boolean numeric = rowSym >= '0' && rowSym <= '9'
              || rowSym == '-' && first.length() > 1 && first.charAt(1) >= '0' && first.charAt(1) <= '9';
          if (numeric && i > 0) { memoryDeletion(); return false; }
          if (rowSym >= 'A' && rowSym <= 'Z') {
            if (vhash[rowSym - 'A'] == -1) { memoryDeletion(); return false; }
          } else if (numeric) {
            firstColumn = false;
            j = 1;
            simMat[i][0] = (float) Utils.atof(first);
            rowSym = listSym.charAt(i);
            first = rowSym + first.substring(1);
          }
        } else {
          j = 1;
          // C++ does not uppercase this branch: a lower-case token is passed to
          // atof and becomes zero rather than triggering the uppercase check.
          if (rowSym >= 'A' && rowSym <= 'Z' && i > 0) { memoryDeletion(); return false; }
          simMat[i][0] = (float) Utils.atof(first);
          rowSym = listSym.charAt(i);
          first = rowSym + (first.isEmpty() ? "" : first.substring(1));
        }
        int row = rowSym >= 'A' && rowSym <= 'Z' ? vhash[rowSym - 'A'] : -1;
        // Invalid row symbols index arbitrary C++ memory. Reject rather than
        // inventing a Java exception or claiming deterministic UB compatibility.
        if (row < 0) { memoryDeletion(); return false; }
        for (; j < numPositions; j++) simMat[row][j] = input.number(simMat[row][j]);
      }

      symmetrize();
      computeDistances();
      return true;
    } catch (IOException e) {
      return false;
    }
  }

  private static char asciiUpper(char c) { return c >= 'a' && c <= 'z' ? (char) (c - 32) : c; }

  /** C-locale libstdc++ formatted extraction with sticky EOF/fail state. */
  private static final class MatrixInput {
    final String data;
    int pos;
    boolean failed;
    MatrixInput(String data, int pos) { this.data = data; this.pos = pos; }
    private boolean ready() {
      if (failed) return false;
      while (pos < data.length() && WHITESPACE.indexOf(data.charAt(pos)) >= 0) pos++;
      if (pos == data.length()) { failed = true; return false; }
      return true;
    }
    String word(String previous) {
      if (!ready()) return previous;
      int begin = pos;
      while (pos < data.length() && WHITESPACE.indexOf(data.charAt(pos)) < 0) pos++;
      return data.substring(begin, pos);
    }
    float number(float previous) {
      if (!ready()) return previous;
      int begin = pos;
      if (data.charAt(pos) == '+' || data.charAt(pos) == '-') pos++;
      int digits = 0;
      while (pos < data.length() && digit(data.charAt(pos))) { pos++; digits++; }
      if (pos < data.length() && data.charAt(pos) == '.') {
        pos++;
        while (pos < data.length() && digit(data.charAt(pos))) { pos++; digits++; }
      }
      if (digits == 0) { failed = true; return 0; }
      if (pos < data.length() && (data.charAt(pos) == 'e' || data.charAt(pos) == 'E')) {
        pos++;
        if (pos < data.length() && (data.charAt(pos) == '+' || data.charAt(pos) == '-')) pos++;
        int exponentStart = pos;
        while (pos < data.length() && digit(data.charAt(pos))) pos++;
        if (exponentStart == pos) { failed = true; return 0; }
      }
      String token = data.substring(begin, pos);
      float value;
      try { value = Float.parseFloat(token); }
      catch (NumberFormatException error) { failed = true; return 0; }
      // num_get<float> clamps overflow and sets failbit; atof's HUGE_VAL differs.
      if (Float.isInfinite(value)) {
        failed = true;
        return value < 0 ? -Float.MAX_VALUE : Float.MAX_VALUE;
      }
      return value;
    }
    private static boolean digit(char c) { return c >= '0' && c <= '9'; }
  }

  private void symmetrize() {
    for (int i = 0; i < numPositions; i++) {
      for (int j = i + 1; j < numPositions; j++) {
        if (simMat[i][j] != simMat[j][i]) {
          float value = (float) ((simMat[i][j] + simMat[j][i]) / 2.0);
          simMat[i][j] = value;
          simMat[j][i] = value;
        }
      }
    }
  }

  /** j outer, i inner, k ascending — the float sum depends on this order. */
  private void computeDistances() {
    for (int j = 0; j < numPositions; j++) {
      for (int i = 0; i < numPositions; i++) {
        if ((i != j) && (distMat[i][j] == 0.0)) {
          float sum = 0;
          for (int k = 0; k < numPositions; k++) {
            float d = simMat[k][j] - simMat[k][i];
            sum += d * d;
          }
          sum = (float) Math.sqrt(sum);
          distMat[i][j] = sum;
          distMat[j][i] = sum;
        }
      }
    }
  }

  public void defaultAASimMatrix() {
    memoryAllocation(20);
    for (int i = 0; i < Utils.TAMABC; i++) vhash[i] = -1;
    for (int i = 0; i < numPositions; i++)
      vhash[Values.listAASym.charAt(i) - 'A'] = i;
    for (int i = 0; i < numPositions; i++)
      for (int j = 0; j < numPositions; j++)
        simMat[i][j] = Values.defaultAAMatrix[i][j];
    computeDistances();
  }

  public void defaultNTSimMatrix() {
    memoryAllocation(5);
    for (int i = 0; i < Utils.TAMABC; i++) vhash[i] = -1;
    for (int i = 0; i < numPositions; i++)
      vhash[Values.listNTSym.charAt(i) - 'A'] = i;
    for (int i = 0; i < numPositions; i++)
      for (int j = 0; j < numPositions; j++)
        simMat[i][j] = Values.defaultNTMatrix[i][j];
    computeDistances();
  }

  public void defaultNTDegeneratedSimMatrix() {
    memoryAllocation(15);
    for (int i = 0; i < Utils.TAMABC; i++) vhash[i] = -1;
    for (int i = 0; i < numPositions; i++)
      vhash[Values.listNTDegenerateSym.charAt(i) - 'A'] = i;
    for (int i = 0; i < numPositions; i++)
      for (int j = 0; j < numPositions; j++)
        simMat[i][j] = Values.defaultNTDegeneratedMatrix[i][j];
    computeDistances();
  }

  public void alternativeSimilarityMatrices(int matrixCode, int datatype) {
    switch (datatype) {
      case Utils.AAType:
        memoryAllocation(20);
        break;
      case Utils.DNAType:
      case Utils.RNAType:
        memoryAllocation(5);
        break;
      case Utils.DNADeg:
      case Utils.RNADeg:
        memoryAllocation(15);
        break;
      default:
        break;
    }

    for (int i = 0; i < Utils.TAMABC; i++) vhash[i] = -1;

    for (int i = 0; i < numPositions; i++) {
      switch (datatype) {
        case Utils.AAType:
          vhash[Values.listAASym.charAt(i) - 'A'] = i;
          break;
        case Utils.DNAType:
        case Utils.RNAType:
          vhash[Values.listNTSym.charAt(i) - 'A'] = i;
          break;
        case Utils.DNADeg:
        case Utils.RNADeg:
          vhash[Values.listNTDegenerateSym.charAt(i) - 'A'] = i;
          break;
        default:
          break;
      }
    }

    for (int i = 0; i < numPositions; i++) {
      for (int j = 0; j < numPositions; j++) {
        switch (matrixCode) {
          case 1:
            simMat[i][j] = alternativeIdentityValue(i, j);
            break;
          default:
            break;
        }
      }
    }

    computeDistances();
  }

  private static float alternativeIdentityValue(int i, int j) {
    // The AA alternative option allocates 20x20 but the upstream constant is
    // only 15x15. That C++ indexing is UB. In the shipped g++ -O2 ELF the array
    // is contiguous, followed by 7 zero padding floats and defaultAAMatrix.
    // Preserve that measured reference layout rather than emitting a Java OOB
    // stacktrace. This compatibility rule is not portable C++ semantics.
    int offset = i * 15 + j;
    if (offset < 225) return Values.alternative1NTDegeneratedMatrix[offset / 15][offset % 15];
    if (offset < 232) return 0;
    offset -= 232;
    return Values.defaultAAMatrix[offset / 20][offset % 20];
  }

  /** cerr << setw(8) << setprecision(4) << right << simMat[i][j]. */
  public void printMatrix() throws IOException {
    try (CppOut err = CppOut.toStderr()) {
      for (int i = 0; i < numPositions; i++) {
        for (int j = 0; j < numPositions; j++)
          err.setw(8).right().precision(4).p(simMat[i][j]);
        err.endl();
      }
    }
  }

  public float getDistance(char a, char b) throws IOException {
    char chA = asciiUpper(a);
    char chB = asciiUpper(b);

    int numa, numb;
    if (chA >= 'A' && chA <= 'Z') numa = vhash[chA - 'A'];
    else {
      try (CppOut err = CppOut.toStderr()) {
        err.p("Error: the symbol '").p(a).p("' is incorrect").endl();
      }
      return -1;
    }

    if (chB >= 'A' && chB <= 'Z') numb = vhash[chB - 'A'];
    else {
      try (CppOut err = CppOut.toStderr()) {
        err.p("Error: the symbol '").p(b).p("' is incorrect").endl();
      }
      return -1;
    }

    if (numa == -1) {
      try (CppOut err = CppOut.toStderr()) {
        err.p("Error: the symbol '").p(a).p("' accesing the matrix is not defined in this object").endl();
      }
      return -1;
    }

    if (numb == -1) {
      try (CppOut err = CppOut.toStderr()) {
        err.p("Error: the symbol '").p(b).p("' accesing the matrix is not defined in this object").endl();
      }
      return -1;
    }

    return distMat[numa][numb];
  }
}
