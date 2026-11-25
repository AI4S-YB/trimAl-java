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

package trimal.util;

/**
 * Faithful transcription of trimAl v1.5.1 utils.cpp (and defines.h constants),
 * including its own quicksort (last-element pivot, order-significant for int[][]),
 * roundInt((int)(x+0.5)) and the character-level helpers.
 */
public final class Utils {

  // defines.h
  public static final int DNAType = 1;
  public static final int RNAType = 2;
  public static final int AAType  = 3;
  public static final int DNADeg  = 4;
  public static final int RNADeg  = 5;

  public static final int SINGLE = 1;
  public static final int MULTI  = 2;

  public static final int GAPPYOUT = 1;
  public static final int STRICT   = 2;

  public static final String DELIMITERS     = "   \t\n";
  public static final String OTHDELIMITERS  = "   \t\n,:";
  public static final String OTH2DELIMITERS = "   \n,:;";

  public static final int HTMLBLOCKS = 120;
  public static final int PHYLIPDISTANCE = 10;

  // similarityMatrix.cpp local defines
  public static final int NUMAMINOS = 20;
  public static final int TAMABC = 28;
  public static final int LINE_LENGTH = 256;
  public static final int REFER = 65;

  private Utils() {}

  public static void initlVect(int[] vector, int tam, int valor) {
    for (int i = 0; i < tam; i++) vector[i] = valor;
  }

  public static void initlVect(float[] vector, int tam, float valor) {
    for (int i = 0; i < tam; i++) vector[i] = valor;
  }

  public static void copyVect(int[] vect1, int[] vect2, int tam) {
    for (int i = 0; i < tam; i++) vect2[i] = vect1[i];
  }

  public static void copyVect(float[] vect1, float[] vect2, int tam) {
    for (int i = 0; i < tam; i++) vect2[i] = vect1[i];
  }

  /** (int) number — truncation toward zero. */
  public static int roundToInf(double number) {
    return x86IntCast(number);
  }

  /** (int)((double) number + 0.5) — half-up for non-negative values. */
  public static int roundInt(double number) {
    return x86IntCast(number + 0.5);
  }

  /** (int)((double) number + 1.0). */
  public static int roundToSup(double number) {
    return x86IntCast(number + 1.0);
  }

  /**
   * The reference g++/x86-64 runtime uses cvttsd2si: out-of-range and NaN produce
   * INT_MIN. Such C++ casts are undefined by the language; this is the observed
   * binary behavior, not a portable definition of an out-of-range C++ cast.
   */
  private static int x86IntCast(double number) {
    if (!(number >= -2147483649.0 && number < 2147483648.0))
      return Integer.MIN_VALUE;
    return (int) number;
  }

  /** The quiet NaN an invalid float operation produces on the x86-64 reference:
   *  sign bit set, printed by the reference as "-nan". */
  private static final float X86_INDEFINITE_NAN = Float.intBitsToFloat(0xFFC00000);

  /**
   * float division with the NaN of the x86-64 reference. An invalid division with no
   * NaN operand (0/0, inf/inf) yields that platform's "indefinite" quiet NaN there,
   * while the JVM leaves the sign of such a NaN to the CPU - and ARM's default NaN
   * is positive, so the port would print "nan" where the reference prints "-nan".
   * A NaN operand propagates unchanged, as it does on both. On x86 this is the
   * hardware result anyway.
   */
  public static float x86Div(float a, float b) {
    float q = a / b;
    if (q != q && a == a && b == b)
      return X86_INDEFINITE_NAN;
    return q;
  }

  /**
   * std::string operator[] semantics for column-indexed reads: index == size() is
   * defined in C++11 and yields the null character. main() only allows
   * -clusters/-maxidentity on unaligned input, whose shorter sequences are then read
   * up to the first sequence's length; beyond-size reads stay '\0' here (native
   * reads heap bytes there — noted divergence, unreachable for the shipped options).
   */
  public static char charAt(String s, int i) {
    return i < s.length() ? s.charAt(i) : '\0';
  }

  public static int max(int x, int y) { return x > y ? x : y; }
  public static float max(float x, float y) { return x > y ? x : y; }
  public static double max(double x, double y) { return x > y ? x : y; }
  public static int min(int x, int y) { return x < y ? x : y; }
  public static float min(float x, float y) { return x < y ? x : y; }
  public static double min(double x, double y) { return x < y ? x : y; }

  /** utils::isNumber — lenient C++ scan, transcribed statement by statement. */
  public static boolean isNumber(String num) {
    int nul = num.indexOf('\0');
    int tam = nul < 0 ? num.length() : nul;
    int flt = 1, expn = 1, sgn = 1;

    for (int i = 0; i < tam; i++) {
      char c = num.charAt(i);
      if (c == '.' && flt != 0)
        flt = 0;
      else if ((c == 'e' || c == 'E') && expn != 0)
        expn = 0;
      else if ((c == '+' || c == '-') && sgn != 0) {
        if (expn == 0) sgn = 0;
      } else if (c > '9' || c < '0')
        return false;
    }
    return true;
  }

  public static boolean compare(String a, String b) {
    return a.equals(b);
  }

  public static String removeSpaces(String in) {
    StringBuilder outb = new StringBuilder(in.length());
    for (int i = 0; i < in.length(); i++) {
      char c = in.charAt(i);
      if (c != ' ' && c != '\t') outb.append(c);
    }
    return outb.toString();
  }

  // ------------------------------------------------------------------ sorts

  /**
   * The three quicksorts below transcribe utils.cpp:230-374 verbatim (last-element
   * pivot). Recursion is replaced by an explicit stack: partition boundaries depend
   * only on the array contents and subranges never overlap, so the final arrangement
   * is identical to the recursive form — but a near-sorted 100k-column input no
   * longer overflows the thread stack before the reference binary does.
   */
  private static final class IntStack {
    private int[] values = new int[16];
    private int top = 0;

    void push(int v) {
      if (top == values.length) {
        int[] grown = new int[values.length * 2];
        System.arraycopy(values, 0, grown, 0, values.length);
        values = grown;
      }
      values[top++] = v;
    }

    int pop() {
      return values[--top];
    }

    boolean isEmpty() {
      return top == 0;
    }
  }

  public static void quicksort(float[] vect, int ini, int fin) {
    IntStack bounds = new IntStack();
    bounds.push(ini);
    bounds.push(fin);
    while (!bounds.isEmpty()) {
      int hi = bounds.pop();
      int lo = bounds.pop();
      if ((lo >= hi) || (hi < 0)) continue;

      float elem_div = vect[hi];
      int i = lo - 1;
      int j = hi;

      while (true) {
        while (vect[++i] < elem_div)
          if (i == hi) break;

        while (vect[--j] > elem_div)
          if (j == 0) break;

        if (i < j) swapF(vect, i, j);
        else break;
      }

      swapF(vect, i, hi);
      bounds.push(lo);
      bounds.push(i - 1);
      bounds.push(i + 1);
      bounds.push(hi);
    }
  }

  private static void swapF(float[] a, int i, int j) {
    float t = a[i]; a[i] = a[j]; a[j] = t;
  }

  public static void quicksort(int[] vect, int ini, int fin) {
    IntStack bounds = new IntStack();
    bounds.push(ini);
    bounds.push(fin);
    while (!bounds.isEmpty()) {
      int hi = bounds.pop();
      int lo = bounds.pop();
      if ((lo >= hi) || (hi < 0)) continue;

      int elem_div = vect[hi];
      int i = lo - 1;
      int j = hi;

      while (true) {
        while (vect[++i] < elem_div)
          if (i == hi) break;

        while (vect[--j] > elem_div)
          if (j == 0) break;

        if (i < j) swapI(vect, i, j);
        else break;
      }

      swapI(vect, i, hi);
      bounds.push(lo);
      bounds.push(i - 1);
      bounds.push(i + 1);
      bounds.push(hi);
    }
  }

  private static void swapI(int[] a, int i, int j) {
    int t = a[i]; a[i] = a[j]; a[j] = t;
  }

  /** int** overload: pivot on vect[k][0], comparisons promote int to float. */
  public static void quicksort(int[][] vect, int ini, int fin) {
    IntStack bounds = new IntStack();
    bounds.push(ini);
    bounds.push(fin);
    while (!bounds.isEmpty()) {
      int hi = bounds.pop();
      int lo = bounds.pop();
      if ((lo >= hi) || (hi < 0)) continue;

      float elem_div = vect[hi][0];
      int i = lo - 1;
      int j = hi;

      while (true) {
        while (vect[++i][0] < elem_div) if (i == hi) break;
        while (vect[--j][0] > elem_div) if (j == 0) break;

        if (i < j) swapII(vect, i, j);
        else break;
      }

      swapII(vect, i, hi);
      bounds.push(lo);
      bounds.push(i - 1);
      bounds.push(i + 1);
      bounds.push(hi);
    }
  }

  private static void swapII(int[][] a, int i, int j) {
    int[] t = a[i]; a[i] = a[j]; a[j] = t;
  }

  // ------------------------------------------------------------------ misc

  public static String getReverse(String toReverse) {
    StringBuilder line = new StringBuilder(toReverse.length());
    for (int i = toReverse.length() - 1; i >= 0; i--)
      line.append(toReverse.charAt(i));
    return line.toString();
  }

  public static String removeCharacter(char c, String line) {
    int pos = line.indexOf(c, 0);
    StringBuilder sb = new StringBuilder(line);
    while (pos >= 0) {
      sb.deleteCharAt(pos);
      pos = sb.indexOf(String.valueOf(c), pos);
    }
    return sb.toString();
  }

  /**
   * utils::checkTypeAlignment — k is the count of non-gap/non-? chars; an all-gap
   * sequence yields 0/0 = NaN so the first branch is false and nothing is counted
   * (replicates C++ float division semantics).
   */
  public static int checkTypeAlignment(int seqNumber, int residNumber, String[] sequences) {
    String listRNA = "AGCUNagcun";
    String listDNA = "AGCTNagctn";
    String degeneratedCodes = "MmRrWwSsYyKkVvHhDdBb";

    int gDNA = 0, gRNA = 0, extDNA = 0, extRNA = 0;
    for (int i = 0; i < seqNumber; i++) {
      int n = sequences[i].length();
      int hitDNA = 0, hitRNA = 0, degenerate = 0, k = 0;
      for (int j = 0; j < n; j++) {
        char c = sequences[i].charAt(j);
        if (c != '-' && c != '.' && c != '?') {
          k++;
          for (int l = 0; l < listDNA.length(); l++)
            if (listDNA.charAt(l) == c) hitDNA++;
          for (int l = 0; l < listRNA.length(); l++)
            if (listRNA.charAt(l) == c) hitRNA++;
          for (int l = 0; l < degeneratedCodes.length(); l++)
            if (degeneratedCodes.charAt(l) == c) degenerate++;
        }
      }

      float ratioDNA = (float) (degenerate + hitDNA) / k;
      float ratioRNA = (float) (degenerate + hitRNA) / k;

      if (ratioDNA < 0.95 && ratioRNA < 0.95)
        return AAType;
      else if (hitRNA > hitDNA && degenerate == 0)
        gRNA++;
      else if (hitRNA > hitDNA && degenerate != 0)
        extRNA++;
      else if (hitRNA < hitDNA && degenerate == 0)
        gDNA++;
      else if (hitRNA < hitDNA && degenerate != 0)
        extDNA++;
    }

    if (extDNA != 0 && extDNA > extRNA)
      return DNADeg;
    else if (extRNA != 0 && extDNA < extRNA)
      return RNADeg;
    else if (gRNA > gDNA)
      return RNAType;
    else
      return DNAType;
  }

  /** utils::readNumbers_StartEnd — returns null when the line is not exactly "a,b". */
  public static int[] readNumbersStartEnd(String line) {
    int nElems = 0;
    int comma = -1;
    while ((comma = line.indexOf(',', comma + 1)) != -1)
      nElems += 2;

    if (nElems != 2)
      return null;

    int[] numbers = new int[2];
    comma = line.indexOf(",", 0);
    numbers[0] = atoi(line.substring(0, comma));
    numbers[1] = atoi(line.substring(comma + 1));
    return numbers;
  }

  /**
   * utils::readNumbers — expands "a-b,c,..." pairs into [nElems, a1,b1, a2,b2, ...].
   * Returns null on negative or inverted ranges (C++ leaks its static buffer; a fresh
   * array is behaviorally identical for single use).
   */
  public static int[] readNumbers(String line) {
    int nElems = 0;
    int comma = -1;
    while ((comma = line.indexOf(',', comma + 1)) != -1)
      nElems += 2;
    nElems += 2;

    int[] numbers = new int[nElems + 1];
    numbers[0] = nElems;

    int init = 0;
    int i = 1;

    do {
      comma = line.indexOf(",", init);
      int separ = line.indexOf("-", init);

      if (((separ < comma) || (comma == -1)) && (separ != -1)) {
        numbers[i++] = atoi(line.substring(init, separ));
        numbers[i++] = atoi(subUntil(line, separ + 1, comma));
        init = comma + 1;
      } else if ((separ > comma) || (separ == -1)) {
        numbers[i++] = atoi(subUntil(line, init, comma));
        numbers[i++] = atoi(subUntil(line, init, comma));
        init = comma + 1;
      }

      if (numbers[i - 2] < 0)
        return null;
      if (numbers[i - 1] < numbers[i - 2])
        return null;
      if (comma == -1)
        break;
    } while (true);

    return numbers;
  }

  /** line.substr(a, b - a) with C++ substr clamping (b beyond end -> until end). */
  private static String subUntil(String line, int a, int b) {
    if (a > line.length()) a = line.length();
    if (b < 0 || b > line.length()) b = line.length();
    if (a >= b) return "";
    return line.substring(a, b);
  }

  /** Whitespace recognized by C isspace/iostream in the reference C locale. */
  public static boolean isCppWhitespace(char c) {
    return c == ' ' || (c >= '\t' && c <= '\r');
  }

  /**
   * glibc atoi is (int) strtol(..., 10): saturate at signed 64-bit limits, then
   * narrow to 32 bits (not saturation at Integer.MAX_VALUE). Integer overflow is
   * outside the portable atoi contract; these results are calibrated on glibc.
   */
  public static int atoi(String s) {
    int i = 0, n = s.length();
    while (i < n && isCppWhitespace(s.charAt(i))) i++;
    boolean negative = false;
    if (i < n && (s.charAt(i) == '+' || s.charAt(i) == '-')) {
      negative = s.charAt(i++) == '-';
    }
    // Accumulate negatively so Long.MIN_VALUE itself is representable.
    long limit = negative ? Long.MIN_VALUE : -Long.MAX_VALUE;
    long cutoff = limit / 10;
    long value = 0;
    boolean overflow = false;
    while (i < n && s.charAt(i) >= '0' && s.charAt(i) <= '9') {
      int digit = s.charAt(i++) - '0';
      if (!overflow) {
        if (value < cutoff || value * 10 < limit + digit) overflow = true;
        else value = value * 10 - digit;
      }
    }
    if (overflow) return (int) (negative ? Long.MIN_VALUE : Long.MAX_VALUE);
    return (int) (negative ? value : -value);
  }

  /** C atof/strtod in the C locale: longest valid decimal/hexadecimal prefix. */
  public static double atof(String s) {
    int n = s.length(), i = 0;
    while (i < n && isCppWhitespace(s.charAt(i))) i++;
    int start = i;
    boolean negative = false;
    if (i < n && (s.charAt(i) == '+' || s.charAt(i) == '-')) negative = s.charAt(i++) == '-';
    if (s.regionMatches(true, i, "inf", 0, 3))
      return negative ? Double.NEGATIVE_INFINITY : Double.POSITIVE_INFINITY;
    if (s.regionMatches(true, i, "nan", 0, 3)) {
      long payload = 0;
      int begin = i + 3;
      if (begin < n && s.charAt(begin) == '(') {
        int end = s.indexOf(')', begin + 1);
        if (end >= 0) {
          String p = s.substring(begin + 1, end);
          try {
            int radix = 10;
            if (p.startsWith("0x") || p.startsWith("0X")) { p = p.substring(2); radix = 16; }
            else if (p.startsWith("0") && p.length() > 1) { p = p.substring(1); radix = 8; }
            payload = Long.parseUnsignedLong(p, radix);
          } catch (NumberFormatException ignored) { }
        }
      }
      return Double.longBitsToDouble((negative ? Long.MIN_VALUE : 0L)
          | 0x7ff8000000000000L | (payload & 0x0007ffffffffffffL));
    }
    boolean hex = false;
    if (i + 2 <= n && s.charAt(i) == '0' && i + 1 < n
        && (s.charAt(i + 1) == 'x' || s.charAt(i + 1) == 'X')) {
      int next = i + 2;
      if (next < n && (hexDigit(s.charAt(next)) || (s.charAt(next) == '.'
          && next + 1 < n && hexDigit(s.charAt(next + 1))))) {
        hex = true;
        i = next;
      }
    }
    int count = 0;
    while (i < n && (hex ? hexDigit(s.charAt(i)) : decimalDigit(s.charAt(i)))) { i++; count++; }
    if (i < n && s.charAt(i) == '.') {
      i++;
      while (i < n && (hex ? hexDigit(s.charAt(i)) : decimalDigit(s.charAt(i)))) { i++; count++; }
    }
    if (count == 0) return 0;
    boolean exponent = false;
    if (i < n && (hex ? (s.charAt(i) == 'p' || s.charAt(i) == 'P') : (s.charAt(i) == 'e' || s.charAt(i) == 'E'))) {
      int next = i + 1;
      if (next < n && (s.charAt(next) == '+' || s.charAt(next) == '-')) next++;
      int begin = next;
      while (next < n && decimalDigit(s.charAt(next))) next++;
      if (next > begin) { i = next; exponent = true; }
    }
    String token = s.substring(start, i);
    if (hex && !exponent) token += "p0";
    try { return Double.parseDouble(token); }
    catch (NumberFormatException error) { return 0; }
  }

  private static boolean decimalDigit(char c) { return c >= '0' && c <= '9'; }
  private static boolean hexDigit(char c) {
    return decimalDigit(c) || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
  }

  /*
   * Single-precision exp/log kernels adapted from glibc 2.35:
   * sysdeps/ieee754/flt-32/{e_expf,e_exp2f_data,e_logf,e_logf_data}.c
   * https://sourceware.org/git/?p=glibc.git;a=tree;f=sysdeps/ieee754/flt-32;hb=glibc-2.35
   *
   * Copyright (C) 2017-2022 Free Software Foundation, Inc.
   * This file is part of the GNU C Library.
   *
   * The GNU C Library is free software; you can redistribute it and/or
   * modify it under the terms of the GNU Lesser General Public
   * License as published by the Free Software Foundation; either
   * version 2.1 of the License, or (at your option) any later version.
   *
   * The GNU C Library is distributed in the hope that it will be useful,
   * but WITHOUT ANY WARRANTY; without even the implied warranty of
   * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   * Lesser General Public License for more details.
   *
   * You should have received a copy of the GNU Lesser General Public
   * License along with the GNU C Library; if not, see
   * <https://www.gnu.org/licenses/>.
   */
  private static final long[] EXPF_TABLE = {
      0x3ff0000000000000L, 0x3fefd9b0d3158574L, 0x3fefb5586cf9890fL, 0x3fef9301d0125b51L,
      0x3fef72b83c7d517bL, 0x3fef54873168b9aaL, 0x3fef387a6e756238L, 0x3fef1e9df51fdee1L,
      0x3fef06fe0a31b715L, 0x3feef1a7373aa9cbL, 0x3feedea64c123422L, 0x3feece086061892dL,
      0x3feebfdad5362a27L, 0x3feeb42b569d4f82L, 0x3feeab07dd485429L, 0x3feea47eb03a5585L,
      0x3feea09e667f3bcdL, 0x3fee9f75e8ec5f74L, 0x3feea11473eb0187L, 0x3feea589994cce13L,
      0x3feeace5422aa0dbL, 0x3feeb737b0cdc5e5L, 0x3feec49182a3f090L, 0x3feed503b23e255dL,
      0x3feee89f995ad3adL, 0x3feeff76f2fb5e47L, 0x3fef199bdd85529cL, 0x3fef3720dcef9069L,
      0x3fef5818dcfba487L, 0x3fef7c97337b9b5fL, 0x3fefa4afa2a490daL, 0x3fefd0765b6e4540L
  };

  /** The C++ exp(float) overload, not a double exp followed by a float cast. */
  public static strictfp float expf(float x) {
    int raw = Float.floatToRawIntBits(x);
    int top = (raw >>> 20) & 0x7ff;
    if (top >= (Float.floatToRawIntBits(88.0f) >>> 20)) {
      if (raw == 0xff800000) return 0;
      if (top >= 0x7f8) return x + x;
      if (x > 0x1.62e42ep6f) return Float.POSITIVE_INFINITY;
      if (x < -0x1.9fe368p6f) return 0;
    }
    double z = (0x1.71547652b82fep0 * 32) * x;
    double kd = z + 0x1.8p52;
    long ki = Double.doubleToRawLongBits(kd);
    kd -= 0x1.8p52;
    double r = z - kd;
    double scale = Double.longBitsToDouble(EXPF_TABLE[(int) (ki & 31)] + (ki << 47));
    z = (0x1.c6af84b912394p-5 / 32 / 32 / 32) * r + (0x1.ebfce50fac4f3p-3 / 32 / 32);
    double r2 = r * r;
    double y = (0x1.62e42ff0c52d6p-1 / 32) * r + 1;
    y = z * r2 + y;
    return (float) (y * scale);
  }

  private static final double[][] LOGF_TABLE = {
      {0x1.661ec79f8f3bep0, -0x1.57bf7808caadep-2},
      {0x1.571ed4aaf883dp0, -0x1.2bef0a7c06ddbp-2},
      {0x1.49539f0f010bp0, -0x1.01eae7f513a67p-2},
      {0x1.3c995b0b80385p0, -0x1.b31d8a68224e9p-3},
      {0x1.30d190c8864a5p0, -0x1.6574f0ac07758p-3},
      {0x1.25e227b0b8eap0, -0x1.1aa2bc79c81p-3},
      {0x1.1bb4a4a1a343fp0, -0x1.a4e76ce8c0e5ep-4},
      {0x1.12358f08ae5bap0, -0x1.1973c5a611cccp-4},
      {0x1.0953f419900a7p0, -0x1.252f438e10c1ep-5},
      {0x1p0, 0x0p0},
      {0x1.e608cfd9a47acp-1, 0x1.aa5aa5df25984p-5},
      {0x1.ca4b31f026aap-1, 0x1.c5e53aa362eb4p-4},
      {0x1.b2036576afce6p-1, 0x1.526e57720db08p-3},
      {0x1.9c2d163a1aa2dp-1, 0x1.bc2860d22477p-3},
      {0x1.886e6037841edp-1, 0x1.1058bc8a07ee1p-2},
      {0x1.767dcf5534862p-1, 0x1.4043057b6ee09p-2}
  };

  private static strictfp float logf(float x) {
    int ix = Float.floatToRawIntBits(x);
    if (ix == 0x3f800000) return 0;
    if (Integer.compareUnsigned(ix - 0x00800000, 0x7f800000 - 0x00800000) >= 0) {
      if (ix * 2 == 0) return Float.NEGATIVE_INFINITY;
      if (ix == 0x7f800000) return x;
      if ((ix & 0x80000000) != 0 || Integer.compareUnsigned(ix * 2, 0xff000000) >= 0)
        return Float.intBitsToFloat(0x7fc00000);
      ix = Float.floatToRawIntBits(x * 0x1p23f) - (23 << 23);
    }
    int tmp = ix - 0x3f330000;
    int i = (tmp >>> 19) & 15;
    int k = tmp >> 23;
    int iz = ix - (tmp & (0x1ff << 23));
    double z = Float.intBitsToFloat(iz);
    double r = z * LOGF_TABLE[i][0] - 1;
    double y0 = LOGF_TABLE[i][1] + k * 0x1.62e42fefa39efp-1;
    double r2 = r * r;
    double y = 0x1.5575b0be00b6ap-2 * r - 0x1.ffffef20a4123p-2;
    y = -0x1.00ea348b88334p-2 * r2 + y;
    return (float) (y * r2 + (y0 + r));
  }

  /*
   * e_log10f.c -- float version of e_log10.c.
   *
   * ====================================================
   * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
   *
   * Developed at SunPro, a Sun Microsystems, Inc. business.
   * Permission to use, copy, modify, and distribute this
   * software is freely granted, provided that this notice
   * is preserved.
   * ====================================================
   */
  /** The C++ log10(float) overload used by the strict trimming method. */
  public static strictfp float log10f(float x) {
    int hx = Float.floatToRawIntBits(x);
    // The public glibc wrapper preserves either-sign NaN payloads before the
    // finite kernel; x+x quiets signaling NaNs without canonicalizing the sign.
    if ((hx & 0x7fffffff) > 0x7f800000) return x + x;
    int k = 0;
    if (hx < 0x00800000) {
      if ((hx & 0x7fffffff) == 0) return Float.NEGATIVE_INFINITY;
      if (hx < 0) return Float.intBitsToFloat(0x7fc00000);
      k -= 25;
      x *= 0x1p25f;
      hx = Float.floatToRawIntBits(x);
    }
    if (hx >= 0x7f800000) return x + x;
    k += (hx >>> 23) - 127;
    int i = k >>> 31;
    hx = (hx & 0x007fffff) | ((0x7f - i) << 23);
    float y = k + i;
    x = Float.intBitsToFloat(hx);
    float z = y * Float.intBitsToFloat(0x355427db)
        + Float.intBitsToFloat(0x3ede5bd9) * logf(x);
    return z + y * Float.intBitsToFloat(0x3e9a2080);
  }

  /**
   * utils::determineColor — CLUSTAL coloring scheme, returns one of o/y/b/r/g/m/p/c/w.
   */
  public static char determineColor(char res, String column) {
    if (Character.toUpperCase(res) == 'G')
      return 'o';
    else if (Character.toUpperCase(res) == 'P')
      return 'y';
    else if (res != '-') {
      switch (Character.toUpperCase(res)) {
        case 'W': case 'L': case 'V': case 'I': case 'M': case 'F':
          if (lookForPattern(column, "p", 0.5f)) return 'b';
          else if (lookForPattern(column, "wlvimafcyhp", 0.6f)) return 'b';
          else return 'w';
        case 'A':
          if (lookForPattern(column, "p", 0.5f)) return 'b';
          else if (lookForPattern(column, "wlvimafcyhp", 0.6f)) return 'b';
          else if (lookForPattern(column, "t", 0.85f)) return 'b';
          else if (lookForPattern(column, "s", 0.85f)) return 'b';
          else if (lookForPattern(column, "g", 0.85f)) return 'b';
          else return 'w';
        case 'C':
          if (lookForPattern(column, "p", 0.5f)) return 'b';
          else if (lookForPattern(column, "wlvimafcyhp", 0.6f)) return 'b';
          else if (lookForPattern(column, "s", 0.85f)) return 'b';
          else if (lookForPattern(column, "c", 0.85f)) return 'p';
          else return 'w';
        case 'K': case 'R':
          if (lookForPattern(column, "kr", 0.6f)) return 'r';
          else if (lookForPattern(column, "q", 0.85f)) return 'r';
          else return 'w';
        case 'T':
          if (lookForPattern(column, "ts", 0.5f)) return 'g';
          else if (lookForPattern(column, "wlvimafcyhp", 0.6f)) return 'g';
          else return 'w';
        case 'S':
          if (lookForPattern(column, "ts", 0.5f)) return 'g';
          else if (lookForPattern(column, "wlvimafcyhp", 0.8f)) return 'g';
          else return 'w';
        case 'N':
          if (lookForPattern(column, "n", 0.5f)) return 'g';
          else if (lookForPattern(column, "d", 0.85f)) return 'g';
          else return 'w';
        case 'Q':
          if (lookForPattern(column, "qe", 0.5f)) return 'g';
          else if (lookForPattern(column, "kr", 0.6f)) return 'g';
          else return 'w';
        case 'D':
          if (lookForPattern(column, "de", 0.5f)) return 'm';
          else if (lookForPattern(column, "n", 0.5f)) return 'm';
          else return 'w';
        case 'E':
          if (lookForPattern(column, "de", 0.5f)) return 'm';
          else if (lookForPattern(column, "qe", 0.5f)) return 'm';
          else return 'w';
        case 'H': case 'Y':
          if (lookForPattern(column, "p", 0.5f)) return 'c';
          else if (lookForPattern(column, "wlvimafcyhp", 0.5f)) return 'c';
          else return 'w';
        default:
          break;
      }
    }
    return 'w';
  }

  public static boolean lookForPattern(String column, String dataset, float level) {
    float count = 0;
    for (int i = 0; i < column.length(); i++) {
      for (int j = 0; j < dataset.length(); j++) {
        if (Character.toUpperCase(column.charAt(i)) == Character.toUpperCase(dataset.charAt(j))) {
          count++;
          break;
        }
      }
    }
    return (count / column.length()) >= level;
  }

  /**
   * utils::trimLine (utils.cpp:455-536) — removes "..." and [...] comments from a
   * sequence line. Returns null (after printing the exact C++ error) on mismatched
   * comment markers or when nothing remains.
   */
  public static String trimLine(String nlineIn) {
    StringBuilder nline = new StringBuilder(nlineIn);
    int pos, next;

    pos = -1;
    while (true) {
      pos = nline.indexOf("\"", pos + 1);
      if (pos == -1)
        break;

      next = nline.lastIndexOf("\"");

      if (nline.indexOf("\"", pos + 1) == next) {
        nline.delete(pos, next + 1);
        pos = -1;
      }

      if (pos == next) {
        System.err.println("\nERROR: Possible (\") mismatch for comments");
        return null;
      }
    }

    while (true) {
      pos = -1;
      next = -1;
      while ((pos = nline.indexOf("[", pos + 1)) != -1)
        next = pos;

      if (next == -1) {
        if (nline.indexOf("]", 0) == -1)
          break;
        System.err.println("\nERROR: Brackets (]) mismatch found");
        return null;
      }

      pos = nline.indexOf("]", next + 1);
      if (pos == -1) {
        System.err.println("\nERROR: Brackets ([) mismatch found");
        return null;
      }

      nline.delete(next, pos + 1);
    }

    if (nline.length() == 0)
      return null;

    return nline.toString();
  }
}
