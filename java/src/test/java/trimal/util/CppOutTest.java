package trimal.util;

import static org.junit.jupiter.api.Assertions.*;
import java.io.ByteArrayOutputStream;
import java.nio.charset.StandardCharsets;
import org.junit.jupiter.api.Test;

public class CppOutTest {
  @Test
  void widthIsOneShotAndDoesNotTruncate() throws Exception {
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    try (CppOut out = new CppOut(bytes)) {
      out.setw(4).p("x").p("y").left().setw(3).p("z").setw(1).p("long").endl();
    }
    assertEquals("   xyz  long\n", bytes.toString(StandardCharsets.ISO_8859_1));
  }

  @Test
  void precisionAndAlignmentAreSticky() throws Exception {
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    try (CppOut out = new CppOut(bytes)) {
      out.precision(4).fixed().left().setw(9).p(0.5).p("|").setw(9).p(1.0).endl();
    }
    assertEquals("0.5000   |1.0000   \n", bytes.toString(StandardCharsets.ISO_8859_1));
  }

  @Test
  void standardStreamStateBelongsToTheStream() throws Exception {
    CppOut first = CppOut.toStdout();
    int precision = first.precision();
    boolean left = first.isLeft();
    first.precision(10).left();
    first.close();
    CppOut second = CppOut.toStdout();
    assertSame(first, second);
    assertEquals(10, second.precision());
    assertTrue(second.isLeft());
    second.precision(precision);
    if (!left) second.right();
  }
}
