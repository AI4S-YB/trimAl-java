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

import java.util.Arrays;
import java.util.Locale;
import trimal.io.Tok;
import trimal.stat.SequencesMatrix;

/** Public/direct-function probe for semantics not exposed by trimal's CLI. */
public final class AuditProbe {
  public static void main(String[] args) throws Exception {
    switch (args[0]) {
      case "identity": {
        trimal.util.CppOut.toStdout().precision(9);
        for (int i = 1; i < args.length; i++) {
          Alignment alignment = new Alignment();
          if (!alignment.loadAlignment(args[i])) throw new AssertionError("input not loaded");
          alignment.printColumnsIdentity_DescriptiveStats();
        }
        break;
      }
      case "matrix": {
        SequencesMatrix matrix = new SequencesMatrix(new String[]{"-A-C", "AC--", "--AC"},
            new String[]{"alpha", "beta", "gamma"}, 3, 4);
        matrix.setOrder(new int[]{2, 0, 1});
        int[] row = new int[4];
        for (String name : new String[]{"alpha", "beta", "gamma"}) {
          if (!matrix.getSequence(name, row)) throw new AssertionError(name);
          System.out.print(name);
          for (int cell : row) System.out.print(" " + cell);
          System.out.println();
        }
        matrix.setOrder(new int[]{1, 1, 2});
        for (String name : new String[]{"alpha", "beta", "gamma"}) {
          matrix.getSequence(name, row);
          System.out.print(name);
          for (int cell : row) System.out.print(" " + cell);
          System.out.println();
        }
        break;
      }
      case "locale": {
        Locale.setDefault(new Locale("tr", "TR"));
        Alignment alignment = new Alignment();
        System.out.println(alignment.loadAlignment(args[1]) ? "loaded" : "not loaded");
        alignment.setOutputFormat(8, false);
        alignment.printAlignment();
        break;
      }
      case "contracts": {
        Tok shortHeader = Tok.strtok(">a", " ");
        if (shortHeader.at(3) != '\0') throw new AssertionError("short header UB fallback");
        Tok rawHeader = Tok.strtok(">a ;desc", " ");
        if (rawHeader.at(3) != ';') throw new AssertionError("raw token remainder");
        Alignment.Cursor withLf = new Alignment.Cursor("abc\n".getBytes("ISO-8859-1"));
        if (!"abc".equals(withLf.readLine()) || withLf.atEof()) throw new AssertionError("EOF before failed read");
        if (withLf.readLine() != null || !withLf.atEof()) throw new AssertionError("EOF after failed read");
        Alignment.Cursor withoutLf = new Alignment.Cursor("abc".getBytes("ISO-8859-1"));
        if (!"abc".equals(withoutLf.readLine()) || !withoutLf.atEof()) throw new AssertionError("missing LF EOF");
        Alignment alignment = new Alignment();
        if (alignment.formatInputAlignment(args[1]) != 0) throw new AssertionError("missing file format");
        System.out.println("contracts passed");
        break;
      }
      default: throw new IllegalArgumentException(args[0]);
    }
  }
}
