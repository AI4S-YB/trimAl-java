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

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;

/**
 * Emulation of C++ iostream floating-point output (operator<< with default flags,
 * arbitrary precision, and the `fixed` manipulator), following printf %g / %f
 * semantics: round the exact binary value to P significant digits (ties-to-even on
 * the exact decimal expansion), then choose fixed vs scientific style based on the
 * exponent of the ROUNDED value, strip trailing zeros (default style only).
 */
public final class CppFormat {

  private CppFormat() {}

  /** Equivalent of `std::ostream << (double)value` with default precision 6. */
  public static String format(double v) {
    return format(v, 6, false);
  }

  /** Equivalent of `std::ostream << (float)value` (promoted to double, precision 6). */
  public static String format(float v) {
    return format(v, 6, false);
  }

  /** Equivalent of `std::ostream << value` with `out.precision(P)` in effect. */
  public static String format(double v, int P) {
    return format(v, P, false);
  }

  /** Equivalent of `std::ostream << (float)value` with `out.precision(P)` in effect. */
  public static String format(float v, int P) {
    return format(v, P, false);
  }

  /** Equivalent of `std::ostream << value` after the `std::fixed` manipulator. */
  public static String formatFixed(double v, int P) {
    return format(v, P, true);
  }

  /** Equivalent of `std::ostream << (float)value` after `std::fixed`. */
  public static String formatFixed(float v, int P) {
    return format(v, P, true);
  }

  private static String format(double v, int P, boolean fixed) {
    boolean neg = (Double.doubleToRawLongBits(v) < 0); // includes signed NaNs and -0.0
    if (Double.isNaN(v)) return neg ? "-nan" : "nan";
    if (Double.isInfinite(v)) return neg ? "-inf" : "inf";
    // libstdc++ forwards precision to printf: a negative precision means the
    // default six digits; %g treats zero precision as one significant digit.
    if (P < 0) P = 6;
    if (!fixed && P == 0) P = 1;
    BigDecimal exact = new BigDecimal(Math.abs(v));
    String s;
    if (exact.signum() == 0) {
      s = fixed ? zeroFixed(P) : "0";
    } else {
      s = fixed ? fixedPositive(exact, P) : formatPositive(exact, P);
    }
    return neg ? "-" + s : s;
  }

  private static String zeroFixed(int P) {
    if (P == 0) return "0";
    StringBuilder sb = new StringBuilder("0.");
    for (int i = 0; i < P; ++i) sb.append('0');
    return sb.toString();
  }

  private static String fixedPositive(BigDecimal exact, int P) {
    // %.*f: exactly P digits after the decimal point, HALF_EVEN rounding of the
    // exact binary expansion (glibc behavior).
    BigDecimal r = exact.setScale(P, RoundingMode.HALF_EVEN);
    return r.toPlainString();
  }

  private static String formatPositive(BigDecimal exact, int P) {
    // Round the exact binary value to P significant digits (HALF_EVEN), like glibc.
    BigDecimal rounded = exact.round(new MathContext(P, RoundingMode.HALF_EVEN));
    if (rounded.signum() == 0) return "0";
    // Scientific exponent of the leading digit of the ROUNDED value.
    int exp = rounded.precision() - rounded.scale() - 1;
    // Significant digit string, no leading zeros; pad to P digits.
    String digits = rounded.unscaledValue().toString();
    StringBuilder d = new StringBuilder(digits);
    while (d.length() < P) d.append('0');
    if (exp >= -4 && exp < P) {
      // Fixed style with P-1-exp digits after the decimal point.
      StringBuilder sb = new StringBuilder();
      if (exp >= 0) {
        sb.append(d, 0, Math.min(exp + 1, P));
        for (int i = P; i < exp + 1; ++i) sb.append('0');
        String frac = exp + 1 < P ? d.substring(exp + 1) : "";
        sb.append('.').append(frac);
      } else {
        sb.append("0.");
        for (int i = 0; i < -exp - 1; ++i) sb.append('0');
        sb.append(d);
      }
      return stripZeros(sb.toString());
    }
    // Scientific style: d.ddd with P-1 digits after the point.
    StringBuilder sb = new StringBuilder();
    sb.append(d.charAt(0)).append('.').append(d, 1, P);
    String mantissa = stripZeros(sb.toString());
    return mantissa + "e" + (exp < 0 ? "-" : "+") + exp2(Math.abs(exp));
  }

  private static String stripZeros(String s) {
    if (s.indexOf('.') < 0) return s;
    int end = s.length();
    while (end > 0 && s.charAt(end - 1) == '0') --end;
    if (end > 0 && s.charAt(end - 1) == '.') --end;
    return s.substring(0, end);
  }

  private static String exp2(int e) {
    String s = Integer.toString(e);
    if (s.length() < 2) s = "0" + s;
    return s;
  }
}
