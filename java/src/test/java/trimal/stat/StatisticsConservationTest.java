package trimal.stat;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import trimal.util.Utils;

public class StatisticsConservationTest {
  private static int[] bits(float[] values) {
    int[] out = new int[values.length];
    for (int i = 0; i < values.length; ++i) out[i] = Float.floatToRawIntBits(values[i]);
    return out;
  }

  @Test
  public void identityConservationAndReflectedWindowMatchOracleBits() throws Exception {
    String[] rows = {"AAC-", "AGC-"};
    StatisticsGaps gaps = new StatisticsGaps(rows, 2, 4, Utils.DNAType);
    assertTrue(gaps.applyWindow(1));
    SimilarityMatrix matrix = new SimilarityMatrix(); matrix.defaultNTSimMatrix();
    StatisticsConservation stats = new StatisticsConservation(rows, 2, 4, Utils.DNAType);
    assertTrue(stats.setSimilarityMatrix(matrix));
    assertTrue(stats.calculateVectors(rows, gaps.getGapsWindow()));
    assertEquals(1107645781, Float.floatToRawIntBits(stats.matrixIdentity[0][1]));
    assertArrayEquals(new int[]{0, 1068827891, 0, 0}, bits(stats.Q));
    assertArrayEquals(new int[]{1065353216, 1048114072, 1065353216, 1065353216}, bits(stats.getMdkVector()));
    assertTrue(stats.applyWindow(1));
    assertArrayEquals(new int[]{1056810632, 1061120419, 1061120419, 1065353216}, bits(stats.getMdkwVector()));
    assertEquals(4602596153937100800L, Double.doubleToRawLongBits(stats.calcCutPoint(0f, .5f)));
    assertFalse(stats.applyWindow(1));
    // Upstream does not clear MDK_Window when the window size changes.
    assertTrue(stats.applyWindow(0));
    assertArrayEquals(new int[]{1069509026, 1065199241, 1071625426, 1073741824}, bits(stats.getMdkwVector()));
  }

  @Test
  public void gapThresholdUsesFloatRatioPromotedToDouble() throws Exception {
    String[] rows = {"-", "-", "-", "-", "A"};
    SimilarityMatrix matrix = new SimilarityMatrix(); matrix.defaultNTSimMatrix();
    StatisticsConservation stats = new StatisticsConservation(rows, 5, 1, Utils.DNAType);
    stats.setSimilarityMatrix(matrix);
    stats.calculateVectors(rows, new int[]{4});
    assertEquals(0f, stats.getMdkVector()[0]);
    // 4/5 is a float slightly ABOVE the double literal 0.8.
    assertFalse(stats.applyWindow(1));
  }

  @Test
  public void emptyPairIdentityIsNanAndMissingSimilarityReturnsFalse() throws Exception {
    StatisticsConservation stats = new StatisticsConservation(new String[]{"--", "NN"}, 2, 2, Utils.DNAType);
    assertTrue(Float.isNaN(stats.matrixIdentity[0][1]));
    assertFalse(stats.calculateVectors(new String[]{"--", "NN"}, null));
    assertFalse(stats.setSimilarityMatrix(null));
  }

  @Test
  public void allAboveThresholdUsesExplicitUbCompatibilityFallback() {
    StatisticsConservation stats = new StatisticsConservation(new String[]{"AAAA", "AAAA"}, 2, 4, Utils.DNAType);
    java.util.Arrays.fill(stats.getMdkwVector(), 1f);
    // Real g++ ordinary-build output is zero in all tested runs; ASan confirms
    // vectAux[-1] UB. The test specifies our fallback, not portable C++ behavior.
    assertEquals(0.0, stats.calcCutPoint(0, 0));
  }
}
