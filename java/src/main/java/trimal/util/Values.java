package trimal.util;

/**
 * Faithful transcription of trimAl v1.5.1 values.h (guarded by SIMMatrix there and
 * compiled once into similarityMatrix.cpp). All matrix values are bit-identical
 * literals; the 1/4. and 1/6. C double divisions are rounded to float on storage.
 */
public final class Values {

  private Values() {}

  public static final String listNTSym = "ACGTU";
  public static final String listAASym = "ARNDCQEGHILKMFPSTWYV";
  public static final String listNTDegenerateSym = "ACGTURYKMSWBDHV";

  public static final String proteinWildcards = "BX";
  public static final String proteinAlternativeAminoacids = "UO";

  /** Default Identity Matrix for Canonical Nucleotides. */
  public static final float[][] defaultNTMatrix = {
    {1, 0, 0, 0, 0},
    {0, 1, 0, 0, 0},
    {0, 0, 1, 0, 0},
    {0, 0, 0, 1, 0},
    {0, 0, 0, 0, 1},
  };

  public static final float[][] defaultNTDegeneratedMatrix = buildNTDegenerated();

  private static float[][] buildNTDegenerated() {
    float q4 = (float) (1 / 4.);
    float q6 = (float) (1 / 6.);
    return new float[][] {
      { 1.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { 0.0f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { 0.0f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { 0.0f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { 0.0f,  0.0f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { q4,    0.0f,  q4,    0.0f,  0.0f,  q4,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { 0.0f,  q4,    0.0f,  q4,    q4,    0.0f,  q4,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { 0.0f,  0.0f,  q4,    q4,    q4,    0.0f,  0.0f,  q4,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { q4,    q4,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  q4,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { 0.0f,  q4,    q4,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  q4,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f},
      { q4,    0.0f,  0.0f,  q4,    q4,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  q4,    0.0f,  0.0f,  0.0f,  0.0f},
      { 0.0f,  q6,    q6,    q6,    q6,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  q6,    0.0f,  0.0f,  0.0f},
      { q6,    0.0f,  q6,    q6,    q6,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  q6,    0.0f,  0.0f},
      { q6,    q6,    0.0f,  q6,    q6,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  q6,    0.0f},
      { q6,    q6,    q6,    0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  0.0f,  q6},
    };
  }

  /** BLOSUM62. */
  public static final float[][] defaultAAMatrix = {
    {  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
    { -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
    { -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
    { -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
    {  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
    { -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
    { -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
    {  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
    { -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
    { -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
    { -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
    { -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
    { -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
    { -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
    { -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
    {  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
    {  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
    { -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
    { -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
    {  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4},
  };

  /** Alternative matrices: 15x15 identity (alternative_1). */
  public static final float[][] alternative1NTDegeneratedMatrix = buildAltIdentity();

  private static float[][] buildAltIdentity() {
    float[][] m = new float[15][15];
    for (int i = 0; i < 15; i++) m[i][i] = 1;
    return m;
  }
}
