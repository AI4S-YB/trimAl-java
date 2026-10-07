package trimal.util;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

/** Golden-value tests calibrated against the reference C++ binary / glibc. */
public class CppFormatTest {

  @Test
  public void testDefaultPrecision6() {
    assertEquals("0.5", CppFormat.format(0.5));
    assertEquals("1", CppFormat.format(1.0));
    assertEquals("0.333333", CppFormat.format(1.0 / 3.0));
    assertEquals("100", CppFormat.format(100.0));
    assertEquals("1e-05", CppFormat.format(1e-5));
    assertEquals("1.23457e+08", CppFormat.format(123456789.0));
    assertEquals("3.14159", CppFormat.format(Math.PI));
    assertEquals("0", CppFormat.format(0.0));
    assertEquals("-0", CppFormat.format(-0.0));
    assertEquals("nan", CppFormat.format(Double.NaN));
    assertEquals("inf", CppFormat.format(Double.POSITIVE_INFINITY));
    assertEquals("-inf", CppFormat.format(Double.NEGATIVE_INFINITY));
    assertEquals("30.4167", CppFormat.format(30.416666666666668));
    assertEquals("0.506944", CppFormat.format(0.5069444444444444));
  }

  @Test
  public void testFloatPromotion() {
    // operator<<(float) promotes to double, but still uses precision 6.
    assertEquals("1.41421", CppFormat.format(1.4142135f));
    // MDK values: float(exp(-x))
    float mdk = (float) Math.exp(-0.25);
    assertEquals("0.778801", CppFormat.format(mdk));
  }

  @Test
  public void testPrecision10() {
    assertEquals("30.41666667", CppFormat.format(730.0 / 24.0, 10));
    assertEquals("0.7419354839", CppFormat.format(23.0 / 31.0, 10));
    assertEquals("1", CppFormat.format(1.0, 10));
    assertEquals("0.1", CppFormat.format(0.1, 10));
    // float MDK bit pattern through precision 10
    assertEquals("0.7788007855", CppFormat.format((float) Math.exp(-0.25), 10));
  }

  @Test
  public void testZeroAndNegativePrecision() {
    assertEquals("1e+02", CppFormat.format(123.45, 0));
    assertEquals("2", CppFormat.format(1.5, 0));
    assertEquals("2", CppFormat.formatFixed(2.5, 0)); // printf uses ties-to-even
    assertEquals("1.23457", CppFormat.format(1.23456789, -1));
    assertEquals("1.234568", CppFormat.formatFixed(1.23456789, -1));
    assertEquals("-0.0000", CppFormat.formatFixed(-0.0, 4));
  }

  @Test
  public void testSignedNan() {
    assertEquals("-nan", CppFormat.format(Double.longBitsToDouble(0xfff8000000000001L)));
    assertEquals("-nan", CppFormat.format(Float.intBitsToFloat(0xffc01234)));
    assertEquals("nan", CppFormat.format(Float.intBitsToFloat(0x7fc01234), 10));
    assertEquals("-nan", CppFormat.formatFixed(Double.longBitsToDouble(0xfff8000000000001L), 4));
  }

  @Test
  public void testFixed4() {
    assertEquals("0.5000", CppFormat.formatFixed(0.5, 4));
    assertEquals("1.0000", CppFormat.formatFixed(1.0, 4));
    assertEquals("-0.3333", CppFormat.formatFixed(-1.0 / 3.0, 4));
    assertEquals("nan", CppFormat.formatFixed(Double.NaN, 4));
  }
}
