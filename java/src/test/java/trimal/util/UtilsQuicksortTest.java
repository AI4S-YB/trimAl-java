package trimal.util;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

/**
 * Quicksort order tests with golden index permutations captured from the compiled
 * reference C++ (verify/oracle/qs_oracle.cpp, g++ -O2). The int** variant is
 * unstable and its tie order is exchange-history dependent — Java must reproduce
 * it element for element.
 */
public class UtilsQuicksortTest {

  private static final int[] TIE_LENS = {5,3,5,5,1,5,3,3,5,1,1,5,3,5,1,3,5,5,3,1,5,3,1,5};
  private static final int[] GOLDEN_TIES = {
      10, 4, 19, 22, 9, 14, 18, 12, 15, 21, 1, 7, 6, 17, 16, 20, 23, 2, 3, 0, 13, 11, 5, 8};
  private static final int[] GOLDEN_REVERSE = {
      23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0};
  private static final int[] GOLDEN_ALL_EQUAL = {
      15, 16, 17, 14, 13, 12, 21, 22, 18, 20, 19, 23, 2, 3, 4, 1, 0, 11, 9, 10, 5, 6, 8, 7};

  private int[][] makeIntPairs(int[] lens) {
    int[][] seqs = new int[lens.length][];
    for (int i = 0; i < lens.length; i++) {
      seqs[i] = new int[2];
      seqs[i][0] = lens[i];
      seqs[i][1] = i;
    }
    return seqs;
  }

  private int[] indexOf(int[][] seqs) {
    int[] out = new int[seqs.length];
    for (int i = 0; i < seqs.length; i++) out[i] = seqs[i][1];
    return out;
  }

  @Test
  public void testIntPairsTies() {
    int[][] seqs = makeIntPairs(TIE_LENS);
    Utils.quicksort(seqs, 0, seqs.length - 1);
    assertArrayEquals(GOLDEN_TIES, indexOf(seqs));
  }

  @Test
  public void testIntPairsReverse() {
    int[] lens = new int[24];
    for (int i = 0; i < 24; i++) lens[i] = 24 - i;
    int[][] seqs = makeIntPairs(lens);
    Utils.quicksort(seqs, 0, seqs.length - 1);
    assertArrayEquals(GOLDEN_REVERSE, indexOf(seqs));
  }

  @Test
  public void testIntPairsAllEqual() {
    int[] lens = new int[24];
    java.util.Arrays.fill(lens, 7);
    int[][] seqs = makeIntPairs(lens);
    Utils.quicksort(seqs, 0, seqs.length - 1);
    assertArrayEquals(GOLDEN_ALL_EQUAL, indexOf(seqs));
  }

  @Test
  public void testFloatSortValues() {
    float[] v = {0.5f, -1f, 0.25f, 0.5f, 0f, -1f};
    Utils.quicksort(v, 0, v.length - 1);
    assertArrayEquals(new float[]{-1f, -1f, 0f, 0.25f, 0.5f, 0.5f}, v);
  }

  @Test
  public void testFloatSortNaNEndures() {
    // IEEE semantics shared by C++ and Java: the transcribed loop turns
    // {1, NaN, 0} into {NaN, 0, 1} (NaN comparisons are always false).
    float[] v = {1f, Float.NaN, 0f};
    Utils.quicksort(v, 0, v.length - 1);
    assertArrayEquals(new float[]{Float.NaN, 0f, 1f}, v);
  }

  @Test
  public void testRound() {
    assertEquals(2, Utils.roundInt(1.5));
    assertEquals(1, Utils.roundInt(1.4));
    assertEquals(1, Utils.roundInt(0.5)); // (int)(0.5 + 0.5)
    assertEquals(0, Utils.roundInt(0.4));
    assertEquals(3, Utils.roundToSup(2.1)); // (int)(2.1 + 1.0)
    assertEquals(2, Utils.roundToInf(2.9)); // plain truncation
  }
}
