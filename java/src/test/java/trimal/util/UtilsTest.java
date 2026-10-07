package trimal.util;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;

/** Values captured from verify/oracle/numeric_oracle.cpp linked to real utils.cpp/libm. */
public class UtilsTest {
  @Test
  public void atoiUsesCWhitespaceAndGlibcNarrowing() {
    assertEquals(-23, Utils.atoi("\t\r\n\u000b\f-23tail"));
    assertEquals(0, Utils.atoi("\u001c12"));
    assertEquals(0, Utils.atoi("\u00a012"));
    assertEquals(Integer.MIN_VALUE, Utils.atoi("2147483648"));
    assertEquals(-2147483647, Utils.atoi("2147483649"));
    assertEquals(-1, Utils.atoi("4294967295"));
    assertEquals(0, Utils.atoi("4294967296"));
    assertEquals(Integer.MAX_VALUE, Utils.atoi("-2147483649"));
    assertEquals(-1, Utils.atoi("9223372036854775808"));
    assertEquals(0, Utils.atoi("-9223372036854775809"));
    assertEquals(12, Utils.atoi("12\0x"));
  }

  @Test
  public void atofParsesLongestDecimalHexAndSpecialPrefix() {
    assertEquals(-23.0, Utils.atof("\t\r\n\u000b\f-23tail"));
    assertEquals(1500.0, Utils.atof("+1.5e3tail"));
    assertEquals(2.0, Utils.atof("2e+"));
    assertEquals(4.0, Utils.atof("0x1p2"));
    assertEquals(-6.0, Utils.atof("-0x1.8p+2"));
    assertEquals(1.0, Utils.atof("0x1p"));
    assertEquals(0.03125, Utils.atof("0x.8p-4"));
    assertEquals(0, Utils.atof(".e4"));
    assertEquals(Double.POSITIVE_INFINITY, Utils.atof("infinity!"));
    assertEquals(Double.NEGATIVE_INFINITY, Utils.atof("-INFINITY"));
    assertEquals(0xfff8000000000000L, Double.doubleToRawLongBits(Utils.atof("-nan(payload)")));
    assertEquals(0x7ff800000000002aL, Double.doubleToRawLongBits(Utils.atof("nan(42)")));
    assertEquals(Long.MIN_VALUE, Double.doubleToRawLongBits(Utils.atof("-0")));
  }

  @Test
  public void lenientNumberScanAndDelimitersMatchCStringRules() {
    assertTrue(Utils.isNumber(""));
    assertTrue(Utils.isNumber("1\0bad"));
    assertTrue(Utils.isNumber("++++"));
    assertFalse(Utils.isNumber("1e--2"));
    assertEquals("   \n,:;", Utils.OTH2DELIMITERS);
    assertArrayEquals(new int[]{Integer.MIN_VALUE, 12}, Utils.readNumbersStartEnd("2147483648,12"));
    assertArrayEquals(new int[]{6, 1, 3, 5, 5, 7, 9}, Utils.readNumbers("1-3,5,7-9"));
  }

  @Test
  public void roundingRetainsTruncationAndObservedX86Overflow() {
    assertEquals(-2, Utils.roundInt(-2.5));
    assertEquals(-1, Utils.roundInt(-1.5));
    assertEquals(0, Utils.roundInt(-.5));
    assertEquals(-2, Utils.roundToInf(-2.9));
    assertEquals(-1, Utils.roundToSup(-2.5));
    assertEquals(Integer.MIN_VALUE, Utils.roundToInf(Double.NaN));
    assertEquals(Integer.MIN_VALUE, Utils.roundToSup(2147483647.0));
    assertEquals(Integer.MIN_VALUE, Utils.roundInt(Double.POSITIVE_INFINITY));
  }

  @Test
  public void expfDiffersFromDoubleThenFloatAtActualOracleInputs() {
    int[][] golden = {
        {1042494126, 1066802561}, {1069920486, 1083567916}, {(int) 3234109872L, 990687722},
        {(int) 3246199761L, 873532611}, {1113197007, 1724354578}, {(int) 3261385200L, 372425331}
    };
    for (int[] g : golden) {
      float x = Float.intBitsToFloat(g[0]);
      assertEquals(g[1], Float.floatToRawIntBits(Utils.expf(x)));
      assertNotEquals(g[1], Float.floatToRawIntBits((float) Math.exp(x)));
    }
    assertEquals(1f, Utils.expf(0));
    assertEquals(0f, Utils.expf(Float.NEGATIVE_INFINITY));
    assertEquals(Float.POSITIVE_INFINITY, Utils.expf(Float.POSITIVE_INFINITY));
    assertEquals(0xffc01234, Float.floatToRawIntBits(Utils.expf(Float.intBitsToFloat(0xffc01234))));
  }

  @Test
  public void log10fMatchesFloatOverloadAndPreservesNaNPayload() {
    int[][] golden = {
        {1093551566, 1065664512}, {1068797687, 1041825377}, {1086697870, 1061844136},
        {1062733066, (int) 3180791931L}, {852212712, (int) 3237216102L}, {1048612640, (int) 3206128784L}
    };
    for (int[] g : golden) {
      float x = Float.intBitsToFloat(g[0]);
      assertEquals(g[1], Float.floatToRawIntBits(Utils.log10f(x)));
      assertNotEquals(g[1], Float.floatToRawIntBits((float) Math.log10(x)));
    }
    assertEquals(0f, Utils.log10f(1));
    assertEquals(Float.NEGATIVE_INFINITY, Utils.log10f(-0f));
    assertEquals(0x7fc00000, Float.floatToRawIntBits(Utils.log10f(-1)));
    assertEquals(0xffc01234, Float.floatToRawIntBits(Utils.log10f(Float.intBitsToFloat(0xffc01234))));
  }
}
