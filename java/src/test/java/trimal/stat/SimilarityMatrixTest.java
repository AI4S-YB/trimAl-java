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

import static org.junit.jupiter.api.Assertions.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import trimal.util.Utils;

/** Matrix parsing/distance golden values taken from the linked real C++ oracle. */
public class SimilarityMatrixTest {
  @TempDir Path temporary;

  private SimilarityMatrix load(String text, boolean expected) throws Exception {
    Path file = Files.createTempFile(temporary, "matrix", ".txt");
    Files.write(file, text.getBytes(StandardCharsets.ISO_8859_1));
    SimilarityMatrix matrix = new SimilarityMatrix();
    assertEquals(expected, matrix.loadSimMatrix(file.toString()));
    return matrix;
  }

  @Test
  public void lowerCaseHeaderAndReorderedRowsAreSupported() throws Exception {
    SimilarityMatrix matrix = load("a c\nC 0 1\nA 1 0\n", true);
    assertEquals(1068827891, Float.floatToRawIntBits(matrix.getDistance('a', 'c')));
    SimilarityMatrix noRows = load("a c\n1 0\n0 1\n", true);
    assertEquals(1068827891, Float.floatToRawIntBits(noRows.getDistance('A', 'C')));
  }

  @Test
  public void rawCarriageReturnInHeaderIsRejectedByReference() throws Exception {
    load("A C\r\nA 1 0\r\nC 0 1\r\n", false);
    load("A A\nA 1 0\nA 0 1\n", false);
    load("A 1\nA 1 0\nC 0 1\n", false);
    load("\u00df A\nA 1 0\nC 0 1\n", false);
  }

  @Test
  public void numericExtractionFailStateIsNotAtOfTokenization() throws Exception {
    for (String invalid : new String[]{"x", "NaN", "1e"}) {
      SimilarityMatrix matrix = load("A C\nA 1 " + invalid + "\nC 3 4\n", true);
      assertEquals(1065353216, Float.floatToRawIntBits(matrix.getDistance('A', 'C')));
    }
    load("A C\nA 1 2junk\nC 3 4\n", false);
    load("A C\nA 1 0x1p2\nC 3 4\n", false);
    SimilarityMatrix noRow = load("A C\n1 2\na 3\n", true);
    assertEquals(1073741824, Float.floatToRawIntBits(noRow.getDistance('A', 'C')));
  }

  @Test
  public void asymmetricAndDefaultDistancesUseFloatOrder() throws Exception {
    SimilarityMatrix matrix = load("A C\nA 1 2\nC 4 3\n", true);
    assertEquals(2f, matrix.getDistance('A', 'C'));
    matrix.defaultNTSimMatrix();
    assertEquals(1068827891, Float.floatToRawIntBits(matrix.getDistance('A', 'C')));
    assertEquals(0f, matrix.getDistance('G', 'G'));
  }

  @Test
  public void aaAlternativePreservesMeasuredReferenceLayoutWithoutJavaOob() throws Exception {
    SimilarityMatrix matrix = new SimilarityMatrix();
    matrix.alternativeSimilarityMatrices(1, Utils.AAType);
    java.lang.reflect.Field field = SimilarityMatrix.class.getDeclaredField("simMat");
    field.setAccessible(true);
    float[][] sim = (float[][]) field.get(matrix);
    assertEquals(20, matrix.getNumPositions());
    assertEquals(1f, sim[0][16]); // C++ [0][16] spills into next 15-wide row
    assertEquals(0f, sim[15][6]); // seven padding words
    assertEquals(4f, sim[15][7]); // defaultAAMatrix follows padding
    assertEquals(-3f, sim[19][19]);
    assertTrue(Float.isFinite(matrix.getDistance('A', 'V')));
  }
}
