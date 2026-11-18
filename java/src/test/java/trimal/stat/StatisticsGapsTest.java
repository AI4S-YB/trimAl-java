package trimal.stat;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.Test;
import trimal.util.Utils;

public class StatisticsGapsTest {
  @Test
  public void reflectedWindowRoundingAndCutPointMatchOracle() throws Exception {
    String[] rows = {"AC--AC--", "A---ACGT", "AC-NACGT"};
    StatisticsGaps gaps = new StatisticsGaps(rows, 3, 8, Utils.DNAType);
    assertArrayEquals(new int[]{0, 1, 3, 2, 0, 0, 1, 1}, gaps.getGapsWindow());
    assertEquals(1, gaps.aminosXInColumn[3]);
    assertTrue(gaps.applyWindow(1));
    assertArrayEquals(new int[]{1, 1, 2, 2, 1, 0, 1, 1}, gaps.getGapsWindow());
    assertEquals(4603579539312869376L, Double.doubleToRawLongBits(gaps.calcCutPoint(50f, .2f)));
    assertEquals(0, gaps.calcCutPoint2ndSlope());
    assertFalse(gaps.applyWindow(1));
    assertFalse(gaps.applyWindow(3));
  }

  @Test
  public void emptyHistogramDoesNotExposeJavaArrayException() throws Exception {
    StatisticsGaps gaps = new StatisticsGaps(new String[]{"", ""}, 2, 0, Utils.DNAType);
    // Upstream reads outside the empty histogram and can print arbitrary slope
    // values; this is a deliberate safe result, not a claimed original match.
    assertEquals(0, gaps.calcCutPoint2ndSlope());
  }

  @Test
  public void allGapColumnsRemainIncludedInHistogram() throws Exception {
    StatisticsGaps gaps = new StatisticsGaps(new String[]{"----", "----"}, 2, 4, Utils.DNAType);
    assertArrayEquals(new int[]{2, 2, 2, 2}, gaps.getGapsWindow());
    assertArrayEquals(new int[]{0, 0, 4}, gaps.numColumnsWithGaps);
    assertEquals(2, gaps.calcCutPoint2ndSlope());
    // The source omits the full-gap histogram bucket from its acum loop and
    // extrapolates using that bucket afterward: 100% yields 2 - (-4/4) == 3.
    assertEquals(3.0, gaps.calcCutPoint(100f, 0f));
  }
}
