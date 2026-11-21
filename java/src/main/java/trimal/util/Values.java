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
