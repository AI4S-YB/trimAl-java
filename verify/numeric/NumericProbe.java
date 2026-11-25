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

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Locale;
import java.util.StringTokenizer;
import trimal.util.CppFormat;
import trimal.util.Utils;

/** Batched bit-level probe of production Java modules, compiled with --release 11. */
public final class NumericProbe {
  private static float asFloat(String raw) { return Float.intBitsToFloat((int) Long.parseUnsignedLong(raw)); }
  private static double asDouble(String raw) { return Double.longBitsToDouble(Long.parseUnsignedLong(raw)); }
  private static String bits(float v) { return Integer.toUnsignedString(Float.floatToRawIntBits(v)); }
  private static String bits(double v) { return Long.toUnsignedString(Double.doubleToRawLongBits(v)); }
  private static void floats(StringBuilder out, String label, float[] v) {
    out.append('|').append(label).append('=');
    for (int i = 0; i < v.length; ++i) { if (i > 0) out.append(','); out.append(bits(v[i])); }
  }
  private static void ints(StringBuilder out, String label, int[] v) {
    out.append('|').append(label).append('=');
    for (int i = 0; i < v.length; ++i) { if (i > 0) out.append(','); out.append(v[i]); }
  }
  private static void matrix(SimilarityMatrix m, int mode) {
    if (mode == 1) m.defaultNTSimMatrix();
    else if (mode == 3) m.defaultAASimMatrix();
    else if (mode == 4) m.defaultNTDegeneratedSimMatrix();
    else m.alternativeSimilarityMatrices(1, mode - 10);
  }
  private static String run(String line) throws Exception {
    StringTokenizer in = new StringTokenizer(line);
    String cmd = in.nextToken(), id = in.nextToken();
    StringBuilder out = new StringBuilder(id);
    switch (cmd) {
      case "M": {
        int mode = Integer.parseInt(in.nextToken()); String alphabet = in.nextToken(), path = in.nextToken();
        SimilarityMatrix m = new SimilarityMatrix(); boolean ok = true;
        if (mode == 0) ok = m.loadSimMatrix(path); else matrix(m, mode);
        out.append("|ok=").append(ok ? 1 : 0);
        if (ok) {
          float[] d = new float[alphabet.length() * alphabet.length()]; int pos = 0;
          for (char a : alphabet.toCharArray()) for (char b : alphabet.toCharArray()) d[pos++] = m.getDistance(a, b);
          floats(out, "dist", d);
        }
        break;
      }
      case "A": {
        SimilarityMatrix m = new SimilarityMatrix(); m.alternativeSimilarityMatrices(1, Utils.AAType);
        java.lang.reflect.Field sf = SimilarityMatrix.class.getDeclaredField("simMat"), df = SimilarityMatrix.class.getDeclaredField("distMat");
        sf.setAccessible(true); df.setAccessible(true);
        float[] sim = new float[400], dist = new float[400]; int pos = 0;
        float[][] sm = (float[][]) sf.get(m), dm = (float[][]) df.get(m);
        for (int i = 0; i < 20; ++i) for (int j = 0; j < 20; ++j) { sim[pos] = sm[i][j]; dist[pos++] = dm[i][j]; }
        out.append("|size=").append(m.getNumPositions()); floats(out, "sim", sim); floats(out, "dist", dist); break;
      }
      case "U": {
        int cols = Integer.parseInt(in.nextToken()); float v = asFloat(in.nextToken()), threshold = asFloat(in.nextToken());
        String row = new String(new char[cols]).replace('\0', 'A');
        StatisticsConservation c = new StatisticsConservation(new String[]{row, row}, 2, cols, Utils.DNAType);
        java.util.Arrays.fill(c.getMdkwVector(), v);
        out.append("|cut=").append(bits(c.calcCutPoint(0, threshold))); break;
      }
      case "G": {
        int n = Integer.parseInt(in.nextToken()); String[] rows = new String[n]; java.util.Arrays.fill(rows, "");
        StatisticsGaps g = new StatisticsGaps(rows, n, 0, Utils.DNAType);
        out.append("|slope=").append(g.calcCutPoint2ndSlope()); break;
      }
      case "S": {
        int datatype = Integer.parseInt(in.nextToken()), n = Integer.parseInt(in.nextToken()), cols = Integer.parseInt(in.nextToken());
        int hw = Integer.parseInt(in.nextToken()), mode = Integer.parseInt(in.nextToken());
        float min = asFloat(in.nextToken()), sim = asFloat(in.nextToken()), gap = asFloat(in.nextToken());
        int secondWindow = Integer.parseInt(in.nextToken());
        String[] rows = new String[n]; for (int i = 0; i < n; ++i) rows[i] = in.nextToken();
        StatisticsGaps g = new StatisticsGaps(rows, n, cols, datatype); boolean gw = g.applyWindow(hw);
        ints(out, "gaps", g.getGapsWindow());
        out.append("|gw=").append(gw ? 1 : 0).append("|slope=").append(g.calcCutPoint2ndSlope());
        out.append("|gapcut=").append(bits(g.calcCutPoint(min, gap)));
        SimilarityMatrix m = new SimilarityMatrix(); matrix(m, mode);
        StatisticsConservation c = new StatisticsConservation(rows, n, cols, datatype); c.setSimilarityMatrix(m);
        out.append("|vectors=").append(c.calculateVectors(rows, g.getGapsWindow()) ? 1 : 0);
        floats(out, "q", c.Q); floats(out, "mdk", c.MDK);
        float[] identity = new float[n * n]; int pos = 0;
        for (int q = 0; q < n * n; q++) identity[pos++] = c.matrixIdentity[q];
        floats(out, "identity", identity);
        out.append("|cw=").append(c.applyWindow(hw) ? 1 : 0);
        floats(out, "window", c.getMdkwVector());
        out.append("|conscut=").append(bits(c.calcCutPoint(min, sim)));
        if (secondWindow >= 0) {
          out.append("|gw2=").append(g.applyWindow(secondWindow) ? 1 : 0); ints(out, "gaps2", g.getGapsWindow());
          out.append("|cw2=").append(c.applyWindow(secondWindow) ? 1 : 0); floats(out, "window2", c.getMdkwVector());
        }
        break;
      }
      case "F": {
        String type = in.nextToken(), raw = in.nextToken(); int p = Integer.parseInt(in.nextToken()); boolean fixed = in.nextToken().equals("1");
        String text;
        if (type.equals("f")) text = fixed ? CppFormat.formatFixed(asFloat(raw), p) : CppFormat.format(asFloat(raw), p);
        else text = fixed ? CppFormat.formatFixed(asDouble(raw), p) : CppFormat.format(asDouble(raw), p);
        out.append("|text=").append(text); break;
      }
      case "T": {
        float v = asFloat(in.nextToken());
        out.append("|exp=").append(bits(Utils.expf(v))).append("|log10=").append(bits(Utils.log10f(v)))
           .append("|sqrt=").append(bits((float) Math.sqrt(v)))
           .append("|expDouble=").append(bits((float) Math.exp(v)))
           .append("|log10Double=").append(bits((float) Math.log10(v))).append("|floatOverloads=1");
        break;
      }
      case "W": {
        // double pow() at the port's two call sites: the strict method's similarity cut
        // (Alignment.cleanCombMethods) and the standard deviation of the column identities.
        int site = Integer.parseInt(in.nextToken());
        float a = asFloat(in.nextToken()), b = asFloat(in.nextToken());
        if (site == 1) {
          double inic = Utils.log10f(a), fin = Utils.log10f(b);
          double vlr = ((inic - fin) / 10) + fin;
          out.append("|simcut=").append(bits((float) Math.pow(10, vlr)));
        } else {
          out.append("|term=").append(bits((float) Math.pow((a - b), 2)));
        }
        break;
      }
      case "I": {
        String hex = in.nextToken(); StringBuilder s = new StringBuilder();
        for (int i = 0; i < hex.length(); i += 2) s.append((char) Integer.parseInt(hex.substring(i, i + 2), 16));
        out.append("|atoi=").append(Utils.atoi(s.toString())).append("|atof=").append(bits(Utils.atof(s.toString())))
           .append("|number=").append(Utils.isNumber(s.toString()) ? 1 : 0); break;
      }
      case "R": {
        double v = asDouble(in.nextToken());
        out.append("|inf=").append(Utils.roundToInf(v)).append("|nearest=").append(Utils.roundInt(v)).append("|sup=").append(Utils.roundToSup(v)); break;
      }
      case "Q": {
        char type = in.nextToken().charAt(0); int n = Integer.parseInt(in.nextToken());
        if (type == 'f') {
          float[] v = new float[n]; for (int i = 0; i < n; ++i) v[i] = asFloat(in.nextToken());
          Utils.quicksort(v, 0, n - 1); floats(out, "sorted", v);
        } else if (type == 'p') {
          int[][] v = new int[n][2]; for (int i = 0; i < n; ++i) { v[i][0] = Integer.parseInt(in.nextToken()); v[i][1] = i; }
          Utils.quicksort(v, 0, n - 1); int[] indices = new int[n]; for (int i = 0; i < n; ++i) indices[i] = v[i][1]; ints(out, "sorted", indices);
        } else {
          int[] v = new int[n]; for (int i = 0; i < n; ++i) v[i] = Integer.parseInt(in.nextToken());
          Utils.quicksort(v, 0, n - 1); ints(out, "sorted", v);
        }
        break;
      }
      default: throw new IllegalArgumentException(cmd);
    }
    return out.toString();
  }
  public static void main(String[] args) throws Exception {
    Locale.setDefault(Locale.ROOT);
    try (BufferedReader input = new BufferedReader(new InputStreamReader(System.in))) {
      String line;
      while ((line = input.readLine()) != null) {
        try { System.out.println(run(line)); }
        catch (Exception error) { System.out.println(line.split(" ")[1] + "|error=" + error.getClass().getSimpleName()); }
      }
    }
  }
}
