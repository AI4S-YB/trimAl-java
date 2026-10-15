package trimal.util;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.io.FileDescriptor;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;

/**
 * Emulation of std::cout state machine: sticky precision / fixed / left-right flags,
 * one-shot setw field width, endl that flushes. Bytes are written verbatim with an
 * ISO-8859-1 mapping so every byte 0-255 stays transparent (C++ writes bytes).
 */
public final class CppOut implements AutoCloseable {

  private static final Charset LATIN1 = StandardCharsets.ISO_8859_1;

  private final OutputStream out;
  private final boolean stdStream;
  private boolean flushOnEndl = true;
  private int precision = 6;
  private boolean fixed = false;
  private boolean left = false;
  private int width = 0;

  private static final class StandardStreams {
    private static final CppOut STDOUT = new CppOut(
        new BufferedOutputStream(new FileOutputStream(FileDescriptor.out), 1 << 16), true);
    private static final CppOut STDERR = new CppOut(new FileOutputStream(FileDescriptor.err), true);
  }

  /** Formatting flags persist across callers, just as they do on std::cout. */
  public static CppOut toStdout() throws IOException {
    return StandardStreams.STDOUT;
  }

  /** Closing a standard stream wrapper flushes it without closing its descriptor. */
  public static CppOut toStderr() throws IOException {
    return StandardStreams.STDERR;
  }

  /**
   * File-backed writer that accumulates lines and flushes once on close(). C++ endl
   * flushes per line, which costs one write syscall per line on slow filesystems;
   * the byte sequence is identical either way, so file output batches instead.
   */
  public static CppOut toFile(OutputStream out) {
    CppOut fileOut = new CppOut(out, false);
    fileOut.flushOnEndl = false;
    return fileOut;
  }

  public CppOut(OutputStream out) {
    this(out, false);
  }

  private CppOut(OutputStream out, boolean stdStream) {
    this.out = out;
    this.stdStream = stdStream;
  }

  public CppOut precision(int p) { this.precision = p; return this; }
  public int precision() { return precision; }

  public CppOut fixed() { this.fixed = true; return this; }
  public CppOut unsetfixed() { this.fixed = false; return this; }

  public CppOut left() { this.left = true; return this; }
  public CppOut right() { this.left = false; return this; }
  public boolean isLeft() { return left; }

  public CppOut setw(int w) { this.width = w; return this; }

  public CppOut endl() throws IOException {
    out.write('\n');
    if (flushOnEndl) out.flush();
    return this;
  }

  public void flush() throws IOException { out.flush(); }

  public CppOut p(String s) throws IOException {
    emit(s);
    return this;
  }

  /** operator<<(char) */
  public CppOut p(char c) throws IOException {
    emit(Character.toString((char) (c & 0xFFFF)));
    return this;
  }

  public CppOut p(int v) throws IOException {
    emit(Integer.toString(v));
    return this;
  }

  public CppOut p(long v) throws IOException {
    emit(Long.toString(v));
    return this;
  }

  public CppOut p(float v) throws IOException {
    emit(fixed ? CppFormat.formatFixed(v, precision) : CppFormat.format(v, precision));
    return this;
  }

  public CppOut p(double v) throws IOException {
    emit(fixed ? CppFormat.formatFixed(v, precision) : CppFormat.format(v, precision));
    return this;
  }

  private void emit(String s) throws IOException {
    if (width > s.length()) {
      StringBuilder pad = new StringBuilder();
      int n = width - s.length();
      for (int i = 0; i < n; ++i) pad.append(' ');
      if (left) {
        out.write(s.getBytes(LATIN1));
        out.write(pad.toString().getBytes(LATIN1));
      } else {
        out.write(pad.toString().getBytes(LATIN1));
        out.write(s.getBytes(LATIN1));
      }
    } else {
      out.write(s.getBytes(LATIN1));
    }
    width = 0; // setw applies to the next field only
  }

  @Override
  public void close() throws IOException {
    if (stdStream) {
      out.flush(); // never close stdout/stderr: the process writes through them again
    } else {
      out.close();
    }
  }
}
