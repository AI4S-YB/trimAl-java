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

package trimal.io;

import trimal.util.Utils;

/**
 * Emulation of the C char-buffer + strtok environment used by the trimAl loaders.
 * The raw line is kept in a mutable char[] (ISO-8859-1 chars, NUL-terminated) so that
 * pointer-level quirks of the original code (strtok replacing delimiters with NUL,
 * `str + 1` walking past a terminator, firstWord[3] reading past the token end into
 * the raw buffer) reproduce byte for byte.
 */
public final class Tok {

  public final char[] buf;
  public int start;

  public Tok(String line) {
    buf = new char[line.length() + 1];
    for (int i = 0; i < line.length(); i++) buf[i] = line.charAt(i);
    buf[line.length()] = '\0';
    start = 0;
  }

  private Tok(char[] buf, int start) {
    this.buf = buf;
    this.start = start;
  }

  /** strlen: chars until NUL. */
  public int strlen() {
    int p = start;
    while (buf[p] != '\0') p++;
    return p - start;
  }

  public boolean isEmpty() {
    return buf[start] == '\0';
  }

  /** s[i] may read past a token terminator, but not the allocated line buffer. */
  public char at(int i) {
    int p = start + i;
    // C++ firstWord[3] on a short FASTA header can be outside its allocation (UB).
    // Keep raw in-buffer reads; choose NUL for the otherwise undefined read.
    return p >= 0 && p < buf.length ? buf[p] : '\0';
  }

  /** The token text up to its NUL (C string, not the raw remainder). */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    for (int p = start; buf[p] != '\0'; p++) sb.append(buf[p]);
    return sb.toString();
  }

  /** Copy of the raw buffer from start (used where C code reads beyond tokens). */
  public String raw() {
    return toString();
  }

  /** str = str + 1 */
  public Tok advance() {
    return new Tok(buf, start + 1);
  }

  private static boolean isDelim(char c, String delims) {
    return c != '\0' && delims.indexOf(c) >= 0;
  }

  /**
   * strtok(tok, delims): NULL continues the scan of the same underlying buffer.
   * Returns null when no token remains. Delimiters are replaced by NUL in the buffer.
   */
  public static Tok strtok(Tok t, String delims) {
    if (t == null) return null;
    int p = t.start;
    while (isDelim(t.buf[p], delims)) p++;
    if (t.buf[p] == '\0') {
      // no more tokens: subsequent strtok(NULL) on the same buffer also return null
      t.start = p;
      return null;
    }
    int s = p;
    while (t.buf[p] != '\0' && !isDelim(t.buf[p], delims)) p++;
    Tok res = new Tok(t.buf, s);
    if (t.buf[p] != '\0') {
      t.buf[p] = '\0';
      t.start = p + 1;
    } else {
      t.start = p; // ended at the string terminator: scan position stays on the NUL
    }
    return res;
  }

  /** Convenience: first token of a fresh line. */
  public static Tok strtok(String line, String delims) {
    return strtok(new Tok(line), delims);
  }

  /** utils::removeCharacter('-', seq) — removes every '-' (bytes 0-255 safe). */
  public static String removeCharacter(char c, String line) {
    return Utils.removeCharacter(c, line);
  }
}
