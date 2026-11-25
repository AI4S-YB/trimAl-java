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

package trimal;

import static org.junit.jupiter.api.Assertions.*;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.DynamicTest;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.TestFactory;
import trimal.util.CppOut;

/**
 * In-process gate against recorded reference behaviour.
 *
 * verify/record_golden.py runs the reference trimal / readAl / statal on the
 * arguments of the differential manifests and stores what they printed and wrote
 * (src/test/resources/golden). This test runs the port on the same arguments, in this
 * JVM, and requires the same exit status, stdout, stderr and produced files byte for
 * byte. Nothing here needs the reference binaries, so it runs wherever the project
 * builds: any operating system, any JDK that can compile it.
 *
 * The working directory is the golden root, a copy of the input files at their
 * repository-relative paths (Maven sets it). That keeps the arguments exactly as
 * recorded: relative paths, `$OUT`, `$HTML` and `$IN`, with no path rewritten and
 * nothing normalised.
 */
public class GoldenCliTest {

  private static final class Case {
    final String id;
    final String tool;
    final int rc;
    final List<String> args = new ArrayList<>();
    final List<String> files = new ArrayList<>();

    Case(String id, String tool, int rc) {
      this.id = id;
      this.tool = tool;
      this.rc = rc;
    }
  }

  private static List<Case> cases;
  private static Map<String, byte[]> expected;

  @BeforeAll
  static void loadRecording() throws IOException {
    cases = new ArrayList<>();
    try (BufferedReader reader = new BufferedReader(new InputStreamReader(
        resource("/golden/cases.txt"), StandardCharsets.UTF_8))) {
      Case current = null;
      for (String line; (line = reader.readLine()) != null; ) {
        if (line.isEmpty() || line.startsWith("#"))
          continue;
        if (line.startsWith("case ")) {
          String[] head = line.split(" ", 4);
          current = new Case(head[1], head[2], Integer.parseInt(head[3]));
          cases.add(current);
        } else if (line.startsWith("arg ")) {
          current.args.add(unescape(line.substring(4)));
        } else if (line.startsWith("file ")) {
          current.files.add(line.substring(5));
        }
      }
    }
    expected = new HashMap<>();
    try (ZipInputStream zip = new ZipInputStream(resource("/golden/expected.zip"))) {
      for (ZipEntry entry; (entry = zip.getNextEntry()) != null; )
        expected.put(entry.getName(), zip.readAllBytes());
    }
  }

  @TestFactory
  Stream<DynamicTest> matchesTheRecordedReferenceOutput() {
    return cases.stream().map(item -> DynamicTest.dynamicTest(item.id, () -> runOne(item)));
  }

  private static void runOne(Case item) throws Exception {
    Path root = Paths.get("").toAbsolutePath().normalize();
    assertTrue(Files.isDirectory(root.resolve("upstream/dataset")),
        "the working directory must be the golden root, src/test/resources/golden/root "
        + "(Maven runs the tests there); it is " + root);
    Path out = Files.createTempDirectory("trimal-golden-");
    try {
      ByteArrayOutputStream stdout = new ByteArrayOutputStream();
      ByteArrayOutputStream stderr = new ByteArrayOutputStream();
      int rc;
      String thrown = null;
      CppOut.redirectStandard(stdout, stderr);
      try {
        String[] argv = new String[item.args.size() + 1];
        argv[0] = item.tool;
        for (int i = 0; i < item.args.size(); i++)
          argv[i + 1] = item.args.get(i).replace("$OUT", out.resolve("o.aln").toString())
              .replace("$HTML", out.resolve("o.html").toString()).replace("$IN", root.toString());
        int code = item.tool.equals("readal") ? ReadAl.run(argv)
            : item.tool.equals("statal") ? StatAl.run(argv) : TrimAl.run(argv);
        rc = code < 0 ? code + 256 : code;   // the exit status main() leaves
      } catch (IOException e) {
        rc = 255;                            // main() exits with 255 on an I/O error
      } catch (RuntimeException | Error e) {
        rc = -1;
        thrown = e.toString();
      } finally {
        CppOut.restoreStandard();
      }

      // Everything the command flushed by the time run() returned is compared as is:
      // the test adds no flush of its own, so output a real process would lose before
      // it could leave the JVM is lost here too, and shows up as a difference.
      assertNull(thrown, item.id + ": the command threw " + thrown
          + " (a real run prints a stack trace and exits with status 1)");
      assertEquals(item.rc, rc, item.id + ": exit status");
      assertArrayEquals(content(item.id + "/stderr"), stderr.toByteArray(), item.id + ": stderr");
      assertArrayEquals(content(item.id + "/stdout"), stdout.toByteArray(), item.id + ": stdout");
      TreeSet<String> produced;
      try (Stream<Path> written = Files.list(out)) {
        produced = written.map(path -> path.getFileName().toString()).collect(Collectors.toCollection(TreeSet::new));
      }
      assertEquals(new TreeSet<>(item.files), produced, item.id + ": produced files");
      for (String name : item.files)
        assertArrayEquals(content(item.id + "/" + name), Files.readAllBytes(out.resolve(name)),
            item.id + ": file " + name);
    } finally {
      try (Stream<Path> walk = Files.walk(out)) {
        for (Path path : walk.sorted(Comparator.reverseOrder()).collect(Collectors.toList()))
          Files.delete(path);
      }
    }
  }

  private static byte[] content(String entry) {
    byte[] data = expected.get(entry);
    return data == null ? new byte[0] : data;
  }

  private static InputStream resource(String name) throws IOException {
    InputStream stream = GoldenCliTest.class.getResourceAsStream(name);
    if (stream == null)
      throw new IOException("missing test resource " + name + " (run verify/record_golden.py)");
    return stream;
  }

  private static String unescape(String text) {
    StringBuilder result = new StringBuilder(text.length());
    for (int i = 0; i < text.length(); i++) {
      char c = text.charAt(i);
      if (c == '\\' && i + 1 < text.length()) {
        char next = text.charAt(++i);
        result.append(next == 'n' ? '\n' : next == 'r' ? '\r' : next);
      } else {
        result.append(c);
      }
    }
    return result.toString();
  }

  /** The seam itself: two commands in a row are formatted like two processes. */
  @Test
  void redirectingTheStandardStreamsStartsFromPristineState() throws Exception {
    CppOut real = CppOut.toStdout();
    ByteArrayOutputStream first = new ByteArrayOutputStream(), second = new ByteArrayOutputStream();
    ByteArrayOutputStream ignored = new ByteArrayOutputStream();
    try {
      CppOut.redirectStandard(first, ignored);
      try (CppOut out = CppOut.toStdout()) {
        out.precision(9).left().fixed().setw(12).p(0.5).endl();
      }
      CppOut.redirectStandard(second, ignored);
      try (CppOut out = CppOut.toStdout()) {
        out.setw(4).p(0.5).endl();
      }
    } finally {
      CppOut.restoreStandard();
    }
    // first: the flags set on the first stream are still sticky within that stream
    // second: a fresh stream is right-aligned, general format, precision 6
    assertEquals("0.500000000 \n", first.toString(StandardCharsets.ISO_8859_1));
    assertEquals(" 0.5\n", second.toString(StandardCharsets.ISO_8859_1));
    assertNotSame(real, CppOut.toStdout());
    assertEquals(6, CppOut.toStdout().precision());
  }
}
