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

import java.io.IOException;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Minimal range-parallel runner for trimAl's O(n^2) statistics loops.
 *
 * Units must be independent: each index writes only distinct result slots and keeps
 * its own float accumulation order internally, so results are bit-identical to the
 * serial transcription regardless of scheduling. Workers are daemons and the calling
 * thread participates, so behavior is unchanged when only one core is available and
 * System.exit is never blocked.
 *
 * -Dtrimal.threads=N overrides the worker count (default: availableProcessors);
 * 1 forces the plain serial loop.
 */
public final class Parallel {

  @FunctionalInterface
  public interface IntBody {
    void accept(int index) throws IOException;
  }

  private Parallel() {}

  public static int threads() {
    int threads = Integer.getInteger("trimal.threads", Runtime.getRuntime().availableProcessors());
    return Math.max(1, threads);
  }

  /** Runs body(i) for every i in [begin, end). */
  public static void forIndex(int begin, int end, IntBody body) throws IOException {
    int units = end - begin;
    int threads = Math.min(threads(), units);
    if (threads <= 1) {
      for (int i = begin; i < end; i++)
        body.accept(i);
      return;
    }

    AtomicInteger next = new AtomicInteger(begin);
    IOException[] failure = new IOException[1];
    Runnable work = () -> {
      try {
        while (true) {
          int i = next.getAndIncrement();
          if (i >= end) return;
          body.accept(i);
        }
      } catch (IOException e) {
        if (failure[0] == null) failure[0] = e;
        next.set(end); // stop the remaining workers
      }
    };

    Thread[] workers = new Thread[threads - 1];
    for (int t = 0; t < workers.length; t++) {
      workers[t] = new Thread(work, "trimal-parallel-" + t);
      workers[t].setDaemon(true);
      workers[t].start();
    }
    work.run(); // the calling thread participates
    for (Thread worker : workers)
      while (worker.isAlive())
        try {
          worker.join();
        } catch (InterruptedException e) {
          // keep waiting for a clean join
        }
    if (failure[0] != null) throw failure[0];
  }
}
