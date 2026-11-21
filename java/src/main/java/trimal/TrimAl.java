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

import java.io.IOException;
import java.util.ArrayList;

import trimal.io.Tok;
import trimal.stat.SimilarityMatrix;
import trimal.stat.SequencesMatrix;
import trimal.util.CppOut;
import trimal.util.Utils;

/**
 * Faithful transcription of trimAl v1.5.1 main.cpp: argument parsing with its exact
 * error paths, post-processing checks, similarity matrix setup, statistics printing,
 * trimming method dispatch, html/backtranslation and output.
 */
public final class TrimAl {

  // defines.h
  static final String BUILD = "2025-11-25";
  static final double VERSION = 1.5;
  static final int REVISION = 1;
  // AUTHORS with UTF-8 bytes preserved through the Latin-1 writer
  static final String AUTHORS = "2009-2020. Victor Fern\u00C3\u00A1ndez-Rodr\u00C3\u00ADguez, "
      + "Toni Gabald\u00C3\u00B3n, and Salvador Capella-Guti\u00C3\u00A9rrez";

  public static void main(String[] args) {
    try {
      // C argv[0] is the program name; java args[] start at argv[1]. Reintroduce the
      // offset so the transcription keeps the original indexing.
      String[] argv = new String[args.length + 1];
      argv[0] = "trimal";
      System.arraycopy(args, 0, argv, 1, args.length);
      // C `return -1` from main exits with status 255 on Linux; Java on Windows
      // would translate -1 differently, so normalize like the shell does.
      int rc = run(argv);
      System.exit(rc < 0 ? rc + 256 : rc);
    } catch (IOException e) {
      System.exit(255);
    }
  }

  static int run(String[] args) throws IOException {
    int argc = args.length;

    boolean appearErrors = false, complementary = false, colnumbering = false, nogaps = false,
        noallgaps = false, gappyout = false, strict = false, strictplus = false, automated1 = false,
        sgc = false, sgt = false, scc = false, sct = false, sfc = false, sft = false, sident = false,
        soverlap = false, selectSeqs = false, selectCols = false, shortNames = false,
        splitbystop = false, terminal = false, keepSeqs = false, keepHeader = false,
        ignorestop = false;

    float conserve = -1, gapThreshold = -1, simThreshold = -1, comThreshold = -1,
        resOverlap = -1, seqOverlap = -1, maxIdentity = -1;

    int outformat = -1, compareset = -1, stats = 0, windowSize = -1, gapWindow = -1,
        simWindow = -1, conWindow = -1, blockSize = -1, clusters = -1, alternative_matrix = -1,
        alignDataType = -1;

    float[] compareVect = null;
    Alignment[] compAlig = null;
    String[] seqNames = null;
    SequencesMatrix seqMatrix = null;
    SimilarityMatrix similMatrix = null;
    Alignment origAlig = null, intermediateAlig = null, singleAlig = null, backtranslation = null;

    int i = 1, num = 0, maxAminos = 0, numfiles = 0, referFile = 0;
    int[] delColumns = null, delSequences = null, seqLengths = null, boundaries = null;
    String forceFile = null, infile = null, backtransFile = null, outfile = null, outhtml = null,
        matrix = null;
    String[] filesToCompare = null;

    if (argc == 1) {
      menu();
      return 0;
    }

    if (args[1].equals("-h") && (2 == argc)) {
      menu();
      examples();
      return 0;
    }

    if (args[1].equals("--version") && (2 == argc)) {
      try (CppOut cout = CppOut.toStdout()) {
        cout.endl().p("trimAl v").p(VERSION).p(".rev").p(REVISION)
            .p(" build[").p(BUILD).p("]").endl().endl();
      }
      return 0;
    }

    origAlig = new Alignment();

    while (i < argc) {

      if (args[i].equals("-in") && (i + 1 != argc) && (infile == null)) {
        if ((sfc) || (sft) || (comThreshold != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of file comparision.").endl().endl();
          }
          appearErrors = true;
          i++;
        } else if ((compareset == -1) || (forceFile != null)) {
          infile = args[++i];
          if (!origAlig.loadAlignment(infile)) {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: Alignment not loaded: \"").p(infile)
                  .p("\" Check the file's content.").endl().endl();
            }
            appearErrors = true;
          }
        } else {
          if (compareset != -1) {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: Option \"").p(args[i])
                  .p("\" not valid. A reference file exists with alignments to compare.").endl().endl();
            }
          }
          if (forceFile != null) {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: Option \"").p(args[i])
                  .p("\" not valid. A alignment file has been setting up to be compare with a set of alignmets.").endl().endl();
            }
          }
          appearErrors = true;
          i++;
        }
      } else if ((args[i].equals("-out")) && (i + 1 != argc) && (outfile == null)) {
        outfile = args[++i];
      } else if ((args[i].equals("-htmlout")) && (i + 1 != argc) && (outhtml == null)) {
        outhtml = args[++i];
      } else if (args[i].equals("-clustal") && (outformat == -1))
        outformat = 1;
      else if (args[i].equals("-fasta") && (outformat == -1))
        outformat = 8;
      else if (args[i].equals("-fasta_m10") && (outformat == -1)) {
        outformat = 8;
        shortNames = true;
      } else if (args[i].equals("-nbrf") && (outformat == -1))
        outformat = 3;
      else if (args[i].equals("-nexus") && (outformat == -1))
        outformat = 17;
      else if (args[i].equals("-mega") && (outformat == -1))
        outformat = 21;
      else if (args[i].equals("-phylip3.2") && (outformat == -1))
        outformat = 11;
      else if (args[i].equals("-phylip3.2_m10") && (outformat == -1)) {
        outformat = 11;
        shortNames = true;
      } else if (args[i].equals("-phylip") && (outformat == -1))
        outformat = 12;
      else if (args[i].equals("-phylip_m10") && (outformat == -1)) {
        outformat = 12;
        shortNames = true;
      } else if (args[i].equals("-phylip_paml") && (outformat == -1))
        outformat = 13;
      else if (args[i].equals("-phylip_paml_m10") && (outformat == -1)) {
        outformat = 13;
        shortNames = true;
      } else if (args[i].equals("-matrix") && (i + 1 != argc) && (matrix == null)) {
        matrix = args[++i];
      } else if (args[i].equals("--alternative_matrix") && (i + 1 != argc) && (alternative_matrix == -1)) {
        i++;
        if (args[i].equals("degenerated_nt_identity"))
          alternative_matrix = 1;
        else {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Alternative not recognized \"").p(args[i]).p("\"").endl().endl();
          }
          appearErrors = true;
        }
      } else if (args[i].equals("-compareset") && (i + 1 != argc) && (compareset == -1)) {
        if (infile == null) {
          java.nio.file.Path p = java.nio.file.Path.of(args[++i]);
          if (!java.nio.file.Files.isRegularFile(p)) {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: Check the reference file with the alignments to compare.").endl().endl();
            }
            appearErrors = true;
          } else {
            byte[] raw = java.nio.file.Files.readAllBytes(p);
            numfiles = countGetlines(raw);
          }
          compareset = i;
        } else {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Option \"").p(args[i])
                .p("\" not valid. A single alignment file has been set by the user.").endl().endl();
          }
          appearErrors = true;
          i++;
        }
      } else if (args[i].equals("-forceselect") && (i + 1 != argc) && (forceFile == null)) {
        if (infile == null) {
          forceFile = args[++i];
          if (!origAlig.loadAlignment(forceFile)) {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: Alignment not loaded: \"").p(forceFile)
                  .p("\" Check the file's content.").endl().endl();
            }
            appearErrors = true;
          }
        } else {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Option \"").p(args[i])
                .p("\" not valid. A single alignment file has been setting it up").endl().endl();
          }
          appearErrors = true;
          i++;
        }
      } else if (args[i].equals("-backtrans") && (i + 1 != argc) && (backtransFile == null)) {
        backtransFile = args[++i];
        backtranslation = new Alignment();
        if (!backtranslation.loadAlignment(backtransFile)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Alignment not loaded: \"").p(backtransFile)
                .p("\" Check the file's content.").endl().endl();
          }
          appearErrors = true;
        }
      } else if ((args[i].equals("-gapthreshold") || args[i].equals("-gt")) && (i + 1 != argc) && (gapThreshold == -1)) {
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[++i])) {
            gapThreshold = (float) (1 - Utils.atof(args[i]));
            if ((gapThreshold < 0) || (gapThreshold > 1)) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The gap threshold value should be between 0 and 1.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The gap threshold value should be a positive real number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-simthreshold") || args[i].equals("-st")) && (i + 1 != argc) && (simThreshold == -1)) {
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[++i])) {
            simThreshold = (float) Utils.atof(args[i]);
            if ((simThreshold < 0) || (simThreshold > 1)) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The similarity threshold value should be between 0 and 1.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The similarity threshold value should be a positive real number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-conthreshold") || args[i].equals("-ct")) && (i + 1 != argc) && (comThreshold == -1)) {
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns.").endl().endl();
          }
          appearErrors = true;
        } else if (infile != null) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination with -in option.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[++i])) {
            comThreshold = (float) Utils.atof(args[i]);
            if ((comThreshold < 0) || (comThreshold > 1)) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The consistency threshold value should be between 0 and 1.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The consistency threshold value should be a positive real number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-cons")) && (i + 1 != argc) && (conserve == -1)) {
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of other manual methods such as manual selection of sequences/columns.").endl().endl();
          }
          appearErrors = true;
        } else if (blockSize != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of column block size value.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[++i])) {
            conserve = (float) Utils.atof(args[i]);
            if ((conserve < 0) || (conserve > 100)) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The minimal positions value should be between 0 and 100.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The minimal positions value should be a positive real number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-selectcols")) && (selectCols == false) && ((i + 3) < argc)
          && (args[++i].equals("{")) && (args[i + 2].equals("}"))) {
        if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed.").endl().endl();
          }
          appearErrors = true;
        } else if (blockSize != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of column block size value.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (comThreshold != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of other manual methods.").endl().endl();
          }
          appearErrors = true;
        } else if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: It's imposible to use this windows size in combination of selection method.").endl().endl();
          }
          appearErrors = true;
        } else if ((delColumns = Utils.readNumbers(args[++i])) == null) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Impossible to parser the sequences number").endl().endl();
          }
          appearErrors = true;
        } else
          selectCols = true;
        i++;
      } else if (args[i].equals("-nogaps") && (!nogaps)) {
        if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of window values.").endl().endl();
          }
          appearErrors = true;
        } else if (blockSize != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of column block size value.").endl().endl();
          }
          appearErrors = true;
        } else if ((noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic methods are not allowed.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else
          nogaps = true;
      } else if (args[i].equals("-noallgaps") && (!noallgaps)) {
        if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of window values.").endl().endl();
          }
          appearErrors = true;
        } else if (blockSize != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of column block size value.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic methods are not allowed.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else
          noallgaps = true;
      } else if (args[i].equals("-keepseqs") && (!keepSeqs)) {
        keepSeqs = true;
      } else if (args[i].equals("-keepheader") && (!keepHeader)) {
        keepHeader = true;
      } else if (args[i].equals("-gappyout") && (!strict)) {
        if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of window values.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic methods are not allowed.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else
          gappyout = true;
      } else if (args[i].equals("-strict") && (!strict)) {
        if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of window values.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic methods are not allowed.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else
          strict = true;
      } else if ((args[i].equals("-strictplus")) && (!strictplus)) {
        if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination with this window value.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic methods are not allowed.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else
          strictplus = true;
      } else if ((args[i].equals("-automated1")) && (!automated1)) {
        if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination with this window value.").endl().endl();
          }
          appearErrors = true;
        } else if (blockSize != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of column block size value.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (comThreshold != -1) || (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic methods are not allowed.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else
          automated1 = true;
      } else if ((args[i].equals("-resoverlap")) && (i + 1 != argc) && (resOverlap == -1)) {
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of methods such as manual selection of sequences/columns.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[++i])) {
            resOverlap = (float) Utils.atof(args[i]);
            if ((resOverlap < 0) || (resOverlap > 1)) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The residue overlap value should be between 0 and 1.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The residue overlap value should be a positive real number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-seqoverlap")) && (i + 1 != argc) && (seqOverlap == -1)) {
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of methods such as manual selection of sequences/columns.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[++i])) {
            seqOverlap = (float) Utils.atof(args[i]);
            if ((seqOverlap < 0) || (seqOverlap > 100)) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The sequences overlap value should be between 0 and 100.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The minimal positions value should be a positive real number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-selectseqs")) && (selectSeqs == false) && ((i + 3) < argc)
          && (args[++i].equals("{")) && (args[i + 2].equals("}"))) {
        if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed.").endl().endl();
          }
          appearErrors = true;
        } else if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (comThreshold != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of other manual methods.").endl().endl();
          }
          appearErrors = true;
        } else if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: It's imposible to use this windows size in combination of selection method.").endl().endl();
          }
          appearErrors = true;
        } else if ((clusters != -1) || (maxIdentity != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Only one method to chose sequences can be applied.").endl().endl();
          }
          appearErrors = true;
        } else if ((delSequences = Utils.readNumbers(args[++i])) == null) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Impossible to parser the sequences number").endl().endl();
          }
          appearErrors = true;
        } else
          selectSeqs = true;
        i++;
      } else if ((args[i].equals("-maxidentity")) && (i + 1 != argc) && (maxIdentity == -1)) {
        if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (comThreshold != -1) || (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of other manual methods such as manual ")
                .p("selection of sequences/columns.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination with window values.").endl().endl();
          }
          appearErrors = true;
        } else if (clusters != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Only one method to chose representative sequences can be applied.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[++i])) {
            maxIdentity = (float) Utils.atof(args[i]);
            if ((maxIdentity < 0) || (maxIdentity > 1)) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The maximum identity threshold should be between 0 and 1.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The minimal positions value should be a positive real number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-clusters")) && (i + 1 != argc) && (clusters == -1)) {
        if ((gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) ||
            (comThreshold != -1) || (selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of other manual methods such as manual ")
                .p("selection of sequences/columns.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Combinations between automatic and manual methods are not allowed").endl().endl();
          }
          appearErrors = true;
        } else if ((windowSize != -1) || (gapWindow != -1) || (simWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination with window values.").endl().endl();
          }
          appearErrors = true;
        } else if (maxIdentity != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Only one method to chose representative sequences can be applied.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[++i])) {
            clusters = Utils.atoi(args[i]);
            if (clusters < 1) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: There is a problem with the given clusters number.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The clusters number should be a positive integer number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-terminalonly")) && (!terminal)) {
        terminal = true;
      } else if ((args[i].equals("--set_boundaries")) && (!terminal) && ((i + 3) < argc)
          && (args[++i].equals("{")) && (args[i + 2].equals("}"))) {
        if ((boundaries = Utils.readNumbersStartEnd(args[++i])) == null) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Impossible to parser the sequences number").endl().endl();
          }
          appearErrors = true;
        }
        terminal = true;
        i++;
      } else if (args[i].equals("-w") && (i + 1 != argc) && (windowSize == -1)) {
        if ((gapWindow != -1) || (simWindow != -1) || (conWindow != -1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination with this specific window value.").endl().endl();
          }
          appearErrors = true;
        }
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: It's imposible to use this windows size in combination of manual selection method.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of automatic methods.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[i + 1])) {
            windowSize = Utils.atoi(args[++i]);
            if (windowSize <= 0) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The window value should be a positive integer number.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The window value should be a number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if (args[i].equals("-gw") && (i + 1 != argc) && (gapWindow == -1)) {
        if (windowSize != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of general window value.").endl().endl();
          }
          appearErrors = true;
        }
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: It's imposible to use this windows size in combination of manual selection method.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of automatic methods.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[i + 1])) {
            gapWindow = Utils.atoi(args[++i]);
            if (gapWindow <= 0) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The window value should be a positive integer number.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The window value should be a number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if (args[i].equals("-sw") && (i + 1 != argc) && (simWindow == -1)) {
        if (windowSize != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of general window value.").endl().endl();
          }
          appearErrors = true;
        }
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: It's imposible to use this windows size in combination of manual selection method.").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of automatic methods.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[i + 1])) {
            simWindow = Utils.atoi(args[++i]);
            if (simWindow <= 0) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The window value should be a positive integer number.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The window value should be a number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if (args[i].equals("-cw") && (i + 1 != argc) && (conWindow == -1)) {
        if (windowSize != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of general window value.").endl().endl();
          }
          appearErrors = true;
        }
        if ((selectCols) || (selectSeqs)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: It's imposible to use this windows size in combination of manual selection method.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[i + 1])) {
            conWindow = Utils.atoi(args[++i]);
            if (conWindow <= 0) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The window value should be a positive integer number.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The window value should be a number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if (args[i].equals("-block") && (i + 1 != argc) && (blockSize == -1)) {
        if (selectCols) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: It's imposible to set a block size value in combination with a column manual selection").endl().endl();
          }
          appearErrors = true;
        } else if (conserve != -1) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: It's imposible to ask for a minimum percentage of the input alignment in combination with column block size").endl().endl();
          }
          appearErrors = true;
        } else if ((nogaps) || (noallgaps)) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination of automatic methods.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (Utils.isNumber(args[i + 1])) {
            blockSize = Utils.atoi(args[++i]);
            if (blockSize <= 0) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The block size value should be a positive integer number.").endl().endl();
              }
              appearErrors = true;
            }
          } else {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The block size value should be a number.").endl().endl();
            }
            appearErrors = true;
          }
        }
      } else if ((args[i].equals("-sgc")) && (!sgc)) {
        sgc = true;
        stats--;
      } else if ((args[i].equals("-sgt")) && (!sgt)) {
        sgt = true;
        stats--;
      } else if ((args[i].equals("-ssc")) && (!scc)) {
        scc = true;
        stats--;
      } else if ((args[i].equals("-sst")) && (!sct)) {
        sct = true;
        stats--;
      } else if ((args[i].equals("-sident")) && (!sident)) {
        sident = true;
        stats--;
      } else if ((args[i].equals("-soverlap")) && (!soverlap)) {
        soverlap = true;
        stats--;
      } else if ((args[i].equals("-sfc")) && (!sfc)) {
        if (infile != null) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination with -in option.").endl().endl();
          }
          appearErrors = true;
          i++;
        } else {
          sfc = true;
          stats--;
        }
      } else if ((args[i].equals("-sft")) && (!sft)) {
        if (infile != null) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("ERROR: Not allowed in combination with -in option.").endl().endl();
          }
          appearErrors = true;
          i++;
        } else {
          sft = true;
          stats--;
        }
      } else if ((args[i].equals("-complementary")) && (complementary == false)) {
        complementary = true;
      } else if ((args[i].equals("-colnumbering")) && (colnumbering == false)) {
        colnumbering = true;
      } else if ((args[i].equals("-splitbystopcodon")) && (splitbystop == false)) {
        splitbystop = true;
      } else if ((args[i].equals("-ignorestopcodon")) && (ignorestop == false)) {
        ignorestop = true;
      } else {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: Parameter \"").p(args[i]).p("\" not valid.").endl().endl();
        }
        appearErrors = true;
      }
      i++;

      if (appearErrors)
        break;
    }

    if ((!appearErrors) && (infile != null) && (forceFile != null)) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: You can not use a single alignmet at the same ")
            .p("time that you force the alignment selection.").endl().endl();
      }
      appearErrors = true;
    }
    if ((!appearErrors) && (compareset == -1) && (forceFile != null)) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: You can not force the alignment selection without set")
            .p(" an alignment dataset against to compare it.").endl().endl();
      }
      appearErrors = true;
    }
    if ((!appearErrors) && (infile == null) && (compareset == -1) && (forceFile == null) && (backtransFile != null)) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: It is impossible to use a Coding Sequences file to apply the back translation method")
            .p(" without define an input alignment.").endl().endl();
      }
      appearErrors = true;
    }
    if ((!appearErrors) && (infile != null)) {
      if (((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1) ||
          (gapThreshold != -1) || (conserve != -1) || (simThreshold != -1) || (selectCols) || (selectSeqs) ||
          (resOverlap != -1) || (seqOverlap != -1) || (stats < 0)) &&
          (!origAlig.isFileAligned())) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The sequences in the input alignment should be aligned in order to use trimming method.").endl().endl();
        }
        appearErrors = true;
      }
    }
    if ((!appearErrors) && (windowSize != -1) && (compareset != -1)) {
      try (CppOut err = CppOut.toStderr()) {
        err.p("INFO: Try with specific comparison file window value. parameter -cw.").endl().endl();
      }
    }

    if ((matrix != null) && (!appearErrors)) {
      if ((!strict) && (!strictplus) && (!automated1) && (simThreshold == -1.0) && (!scc) && (!sct)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The Similarity Matrix can only be used with methods that use this matrix.").endl().endl();
        }
        appearErrors = true;
      }

      if ((gapWindow != -1) || ((compareset == -1) && (conWindow != -1))) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The Similarity Matrix can only be used with general/similarity windows size.").endl().endl();
        }
        appearErrors = true;
      }
    }

    if ((complementary) && (!appearErrors))
      if ((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
          && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (!selectCols) && (!selectSeqs)
          && (resOverlap == -1) && (seqOverlap == -1) && (maxIdentity == -1) && (clusters == -1)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: This parameter can only be used with either an automatic or a manual method.").endl().endl();
        }
        appearErrors = true;
      }

    if ((terminal) && (boundaries != null) && (!appearErrors)) {
      num = origAlig.getNumAminos();

      if ((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
          && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (!selectCols) && (!selectSeqs)
          && (resOverlap == -1) && (seqOverlap == -1) && (maxIdentity == -1) && (clusters == -1)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: This parameter '--set_boundaries' can only be used with either an automatic or a manual method.").endl().endl();
        }
        appearErrors = true;
      } else if (boundaries[1] >= num) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR:  \"--set_boundaries\" parameter only accepts ")
              .p("integer numbers between 0 and the number of positions (").p(num)
              .p(") - 1.").endl().endl();
        }
        appearErrors = true;
      }
    }

    if ((terminal) && (boundaries == null) && (!appearErrors))
      if ((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
          && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (!selectCols) && (!selectSeqs)
          && (resOverlap == -1) && (seqOverlap == -1) && (maxIdentity == -1) && (clusters == -1)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: This parameter '-terminalonly' can only be used with either an automatic or a manual method.").endl().endl();
        }
        appearErrors = true;
      }

    if ((colnumbering) && (!appearErrors)) {
      if ((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1)
          && (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (comThreshold == -1) && (!selectCols) && (!selectSeqs)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: This parameter can only be used with any trimming method.").endl().endl();
        }
        appearErrors = true;
      } else if (stats < 0) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: This parameter is not valid when statistics' parameters are defined.").endl().endl();
        }
        appearErrors = true;
      }
    }

    if ((outhtml != null) && (outfile != null) && (!appearErrors)) {
      if (outhtml.equals(outfile)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The output and html files should not be the same.").endl().endl();
        }
        appearErrors = true;
      }
    }

    if ((outhtml != null) && (!appearErrors)) {
      if ((!nogaps) && (!noallgaps) && (!gappyout) && (!strict) && (!strictplus) && (!automated1) &&
          (gapThreshold == -1) && (conserve == -1) && (simThreshold == -1) && (comThreshold == -1) &&
          (!selectCols) && (!selectSeqs) && (resOverlap == -1) && (seqOverlap == -1) && (maxIdentity == -1) &&
          (clusters == -1)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: This parameter can only be used with any trimming method.").endl().endl();
        }
        appearErrors = true;
      }
    }

    if (((resOverlap != -1) || (seqOverlap != -1)) && (!appearErrors)) {
      if ((resOverlap != -1) && (seqOverlap == -1)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The sequence overlap value should be defined.").endl().endl();
        }
        appearErrors = true;
      } else if ((resOverlap == -1) && (seqOverlap != -1)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: The residue overlap value should be defined.").endl().endl();
        }
        appearErrors = true;
      }
    }

    if ((stats < 0) && (!appearErrors)) {
      stats--;

      if (((nogaps) || (noallgaps) || (gappyout) || (strict) || (strictplus) || (automated1)
          || (gapThreshold != -1) || (conserve != -1) || (simThreshold != -1)) && (outfile == null)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: An output file should be defined in order to get the alignment's statistics.").endl().endl();
        }
        appearErrors = true;
      }
    }

    if ((comThreshold != -1) && (conserve != -1) && (!appearErrors)) {
      if ((gapThreshold != -1) || (simThreshold != -1)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: Combinations among thresholds are not allowed.").endl().endl();
        }
        appearErrors = true;
      }
    }

    /* Files Comparison Methods */
    if ((compareset != -1) && (!appearErrors)) {

      compAlig = new Alignment[numfiles];
      filesToCompare = new String[numfiles];

      ArrayList<String> lines = readGetlines(java.nio.file.Path.of(args[compareset]), numfiles);

      for (i = 0; (i < numfiles) && (!appearErrors); i++) {

        filesToCompare[i] = lines.get(i);

        compAlig[i] = new Alignment();
        if (!compAlig[i].loadAlignment(filesToCompare[i])) {
          try (CppOut err = CppOut.toStderr()) {
            err.endl().p("Alignment not loaded: \"").p(filesToCompare[i])
                .p("\" Check the file's content.").endl().endl();
          }
          appearErrors = true;
        } else {
          if (!compAlig[i].isFileAligned()) {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: The sequences in the input alignment should be aligned in order to use this method.").endl().endl();
            }
            appearErrors = true;
          } else {
            compAlig[i].sequenMatrix();

            if (compAlig[i].getNumAminos() > maxAminos)
              maxAminos = compAlig[i].getNumAminos();

            if ((compAlig[i].getTypeAlignment() != alignDataType) && (alignDataType != -1)) {
              try (CppOut err = CppOut.toStderr()) {
                err.endl().p("ERROR: The alignments' datatypes are different. Check your dataset.").endl().endl();
              }
              appearErrors = true;
            } else
              alignDataType = compAlig[i].getTypeAlignment();
          }
        }
      }

      if ((!appearErrors) && (forceFile == null)) {

        compareVect = new float[maxAminos];
        if ((stats >= 0) && (outfile != null))
          referFile = CompareFiles.algorithm(compAlig, filesToCompare, compareVect, numfiles, true);
        else
          referFile = CompareFiles.algorithm(compAlig, filesToCompare, compareVect, numfiles, false);

        if (windowSize != -1)
          CompareFiles.applyWindow(compAlig[referFile].getNumAminos(), windowSize, compareVect);
        else if (conWindow != -1)
          CompareFiles.applyWindow(compAlig[referFile].getNumAminos(), conWindow, compareVect);

        // C++ ignores algorithm()'s -1 and indexes filesToCompare[-1] (OOB garbage
        // pointer); the observed downstream behavior is a silently failed reload with
        // origAlig left unloaded, so replicate just that.
        if (referFile >= 0)
          origAlig.loadAlignment(filesToCompare[referFile]);

      } else if ((!appearErrors) && (forceFile != null)) {

        compareVect = new float[origAlig.getNumAminos()];
        appearErrors = !(CompareFiles.forceComparison(compAlig, numfiles, origAlig, compareVect));

        if ((windowSize != -1) && (!appearErrors))
          CompareFiles.applyWindow(origAlig.getNumAminos(), windowSize, compareVect);
        else if ((conWindow != -1) && (!appearErrors))
          CompareFiles.applyWindow(origAlig.getNumAminos(), conWindow, compareVect);
      }
    }

    if ((!appearErrors) && (origAlig.getNumAminos() < (blockSize / 4))) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: The block size value is too big. Please, choose another one smaller than residues number / 4.").endl().endl();
      }
      appearErrors = true;
    }

    if ((!appearErrors) && (backtransFile != null)
        && (backtranslation.getTypeAlignment() != Utils.DNAType
            && backtranslation.getTypeAlignment() != Utils.DNADeg)) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: Check your Coding sequences file. It has been detected other kind of biological sequences.").endl().endl();
      }
      appearErrors = true;
    }
    if ((!appearErrors) && (origAlig.isFileAligned() != true) && (backtransFile != null)) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: The input protein file has to be aligned to carry out the backtranslation process").endl().endl();
      }
      appearErrors = true;
    }
    if ((!appearErrors) && (backtransFile == null) && (splitbystop)) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: The -splitbystopcodon parameter can be only set up with backtranslation functionality.").endl().endl();
      }
      appearErrors = true;
    }
    if ((!appearErrors) && (backtransFile == null) && (ignorestop)) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: The -ignorestopcodon parameter can be only set up with backtranslation functionality.").endl().endl();
      }
      appearErrors = true;
    }
    if ((!appearErrors) && (ignorestop) && (splitbystop)) {
      try (CppOut err = CppOut.toStderr()) {
        err.endl().p("ERROR: Incompatibility of -ignorestopcodon & -splitbystopcodon parameters. Choose one.").endl().endl();
      }
      appearErrors = true;
    }
    if ((!appearErrors) && (backtransFile != null)
        && (backtranslation.prepareCodingSequence(splitbystop, ignorestop, origAlig) != true))
      appearErrors = true;

    if ((!appearErrors) && (backtransFile != null)) {
      seqNames = new String[backtranslation.getNumSpecies()];
      seqLengths = new int[backtranslation.getNumSpecies()];
      backtranslation.getSequences(seqNames, seqLengths);

      if (origAlig.checkCorrespondence(seqNames, seqLengths, backtranslation.getNumSpecies(), 3) != true)
        appearErrors = true;
    }

    /* Errors Control */
    if (appearErrors) {
      return -1;
    }

    if (conserve == -1)
      conserve = 0;

    origAlig.trimTerminalGaps(terminal, boundaries);
    origAlig.setKeepSequencesFlag(keepSeqs);
    origAlig.setKeepSeqsHeaderFlag(keepHeader);

    if (windowSize != -1) {
      gapWindow = windowSize;
      simWindow = windowSize;
    } else {
      if (gapWindow == -1)
        gapWindow = 0;
      if (simWindow == -1)
        simWindow = 0;
    }
    origAlig.setWindowsSize(gapWindow, simWindow);

    if (blockSize != -1)
      origAlig.setBlockSize(blockSize);

    if (outformat != -1)
      origAlig.setOutputFormat(outformat, shortNames);

    if ((strict) || (strictplus) || (automated1) || (simThreshold != -1.0) || (scc == true) || (sct == true)) {
      similMatrix = new SimilarityMatrix();

      if (matrix != null)
        similMatrix.loadSimMatrix(matrix);
      else if (alternative_matrix != -1) {
        alignDataType = origAlig.getTypeAlignment();
        similMatrix.alternativeSimilarityMatrices(alternative_matrix, alignDataType);
      } else {
        alignDataType = origAlig.getTypeAlignment();
        if (alignDataType == Utils.AAType)
          similMatrix.defaultAASimMatrix();
        else if ((alignDataType == Utils.DNAType) || (alignDataType == Utils.RNAType))
          similMatrix.defaultNTSimMatrix();
        else if ((alignDataType == Utils.DNADeg) || (alignDataType == Utils.RNADeg))
          similMatrix.defaultNTDegeneratedSimMatrix();
      }

      if (!origAlig.setSimilarityMatrix(similMatrix)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: It's imposible to proccess the Similarity Matrix.").endl().endl();
        }
        return -1;
      }
    }

    if (sgc) {
      origAlig.printStatisticsGapsColumns();
      stats++;
      if (stats < -1) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.endl();
        }
      }
    }

    if (sgt) {
      origAlig.printStatisticsGapsTotal();
      stats++;
      if (stats < -1) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.endl();
        }
      }
    }

    if (scc) {
      origAlig.printStatisticsConservationColumns();
      stats++;
      if (stats < -1) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.endl();
        }
      }
    }

    if (sct) {
      origAlig.printStatisticsConservationTotal();
      stats++;
      if (stats < -1) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.endl();
        }
      }
    }

    if (sident) {
      origAlig.printSeqIdentity();
      stats++;
      if (stats < -1) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.endl();
        }
      }
    }

    if (soverlap) {
      origAlig.printSeqOverlap();
      stats++;
      if (stats < -1) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.endl();
        }
      }
    }

    if (compareset != -1) {
      if (sfc)
        CompareFiles.printStatisticsFileColumns(origAlig.getNumAminos(), compareVect);
      if (sft)
        CompareFiles.printStatisticsFileAcl(origAlig.getNumAminos(), compareVect);
    }

    if (backtransFile != null)
      seqMatrix = origAlig.getSeqMatrix();

    if (nogaps)
      singleAlig = origAlig.cleanGaps(0, 0, complementary);
    else if (noallgaps)
      singleAlig = origAlig.cleanNoAllGaps(complementary);
    else if (gappyout)
      singleAlig = origAlig.clean2ndSlope(complementary);
    else if (strict)
      singleAlig = origAlig.cleanCombMethods(complementary, false);
    else if (strictplus)
      singleAlig = origAlig.cleanCombMethods(complementary, true);
    else if (automated1) {
      if (origAlig.selectMethod() == Utils.GAPPYOUT)
        singleAlig = origAlig.clean2ndSlope(complementary);
      else
        singleAlig = origAlig.cleanCombMethods(complementary, false);
    }

    if (comThreshold != -1)
      singleAlig = origAlig.cleanCompareFile(comThreshold, conserve, compareVect, complementary);

    if ((resOverlap != -1) && (seqOverlap != -1)) {
      intermediateAlig = origAlig.cleanSpuriousSeq(resOverlap, (seqOverlap / 100), complementary);
      singleAlig = intermediateAlig.cleanNoAllGaps(false);
    }

    if (simThreshold != -1.0) {
      if (gapThreshold != -1.0)
        singleAlig = origAlig.clean(conserve, gapThreshold, simThreshold, complementary);
      else
        singleAlig = origAlig.cleanConservation(conserve, simThreshold, complementary);
    } else if (gapThreshold != -1.0)
      singleAlig = origAlig.cleanGaps(conserve, gapThreshold, complementary);

    if ((selectCols) || (selectSeqs)) {

      if (delColumns != null) {
        num = origAlig.getNumAminos();

        for (i = 1; i < delColumns[0] + 1; i++)
          if (delColumns[i] >= num) {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: This option only accepts integer numbers ")
                  .p("between 0 and the number of columns - 1.").endl().endl();
            }
            appearErrors = true;
            break;
          }

        if (!appearErrors)
          singleAlig = origAlig.removeColumns(delColumns, 1, delColumns[0], complementary);
      }

      if (delSequences != null) {
        num = origAlig.getNumSpecies();

        for (i = 1; i < delSequences[0] + 1; i++)
          if (delSequences[i] >= num) {
            try (CppOut err = CppOut.toStderr()) {
              err.endl().p("ERROR: This option only accepts integer numbers ")
                  .p("between 0 and the number of sequences - 1.").endl().endl();
            }
            appearErrors = true;
            break;
          }

        if (!appearErrors) {
          intermediateAlig = origAlig.removeSequences(delSequences, 1, delSequences[0], complementary);
          singleAlig = intermediateAlig.cleanNoAllGaps(false);
        }
      }
    }

    if (maxIdentity != -1) {
      intermediateAlig = origAlig.getClustering(maxIdentity);
      singleAlig = intermediateAlig.cleanNoAllGaps(false);
    } else if (clusters != -1) {
      if (clusters > origAlig.getNumSpecies()) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR:The number of clusters from the alignment can not be larger than the number of sequences from that alignment.").endl().endl();
        }
        appearErrors = true;
      } else {
        intermediateAlig = origAlig.getClustering(origAlig.getCutPointClusters(clusters));
        singleAlig = intermediateAlig.cleanNoAllGaps(false);
      }
    }

    if (singleAlig == null) {
      singleAlig = origAlig;
      origAlig = null;
    }

    if ((outhtml != null) && (!appearErrors))
      if (!origAlig.alignmentSummaryHTML(outhtml, singleAlig.getNumAminos(), singleAlig.getNumSpecies(),
          singleAlig.getCorrespResidues(), singleAlig.getCorrespSequences(), compareVect)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: It's imposible to generate the HTML output file.").endl().endl();
        }
        appearErrors = true;
      }

    if (backtransFile != null) {
      seqNames = new String[singleAlig.getNumSpecies()];
      singleAlig.getSequencesNames(seqNames);

      singleAlig = backtranslation.getTranslationCDS(singleAlig.getNumAminos(), singleAlig.getNumSpecies(),
          singleAlig.getCorrespResidues(), seqNames, seqMatrix, singleAlig);
    }

    if ((outfile != null) && (!appearErrors)) {
      if (!singleAlig.saveAlignment(outfile)) {
        try (CppOut err = CppOut.toStderr()) {
          err.endl().p("ERROR: It's imposible to generate the output file.").endl().endl();
        }
        appearErrors = true;
      }
    } else if ((stats >= 0) && (!appearErrors))
      singleAlig.printAlignment();

    if ((colnumbering) && (!appearErrors))
      singleAlig.printCorrespondence();

    return 0;
  }

  /** while(compare.getline(line, 256)) numfiles++; — counts lines incl. empty ones. */
  private static int countGetlines(byte[] raw) {
    int count = 0;
    int pos = 0;
    boolean any = raw.length > 0;
    while (any) {
      // getline consumes up to '\n' (or EOF); returns false only at EOF with no data
      if (pos >= raw.length) {
        break;
      }
      while (pos < raw.length && raw[pos] != '\n')
        pos++;
      if (pos < raw.length)
        pos++; // consume '\n'
      count++;
      if (pos >= raw.length)
        break;
    }
    return count;
  }

  /**
   * The main reading loop (main.cpp:1342): raw char reads until '\n' or EOF — no
   * stripping at all (unlike utils::readLine). Produces up to numfiles entries.
   */
  private static ArrayList<String> readGetlines(java.nio.file.Path p, int numfiles) throws IOException {
    byte[] raw = java.nio.file.Files.readAllBytes(p);
    ArrayList<String> res = new ArrayList<>();
    StringBuilder nline = new StringBuilder();
    int pos = 0;
    boolean first = true;
    while (pos < raw.length && res.size() < numfiles) {
      if (!first)
        pos++; // consume the '\n' of the previous line
      nline.setLength(0);
      while (pos < raw.length && raw[pos] != '\n')
        nline.append((char) (raw[pos++] & 0xFF));
      res.add(nline.toString());
      first = false;
    }
    return res;
  }

  static void menu() throws IOException {
    try (CppOut cout = CppOut.toStdout()) {
      cout.endl();
      cout.p("trimAl v").p(VERSION).p(".rev").p(REVISION).p(" build[").p(BUILD)
          .p("]. ").p(AUTHORS).endl().endl();

      cout.p("trimAl webpage: https://trimal.readthedocs.io").endl().endl();

      cout.p("This program is free software: you can redistribute it and/or modify ").endl()
          .p("it under the terms of the GNU General Public License as published by ").endl()
          .p("the Free Software Foundation, the last available version.").endl().endl();

      cout.p("Please cite:").endl()
          .p("\t\ttrimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses.")
          .p("\n\t\tSalvador Capella-Gutierrez; Jose M. Silla-Martinez; Toni Gabaldon.")
          .p("\n\t\tBioinformatics 2009, 25:1972-1973.").endl().endl();

      cout.p("Basic usage").endl()
          .p("\ttrimal -in <inputfile> -out <outputfile> -(other options).").endl().endl();

      cout.p("Common options (for a complete list please see the User Guide or visit https://trimal.readthedocs.io):").endl().endl();
      cout.p("    -h                          ").p("Print this information and show some examples.").endl();
      cout.p("    --version                   ").p("Print the trimAl version.").endl().endl();

      cout.p("    -in <inputfile>             ").p("Input file in several formats (clustal, fasta, NBRF/PIR, nexus, phylip3.2, phylip).").endl().endl();

      cout.p("    -compareset <inputfile>     ").p("Input list of paths for the files containing the alignments to compare.").endl();
      cout.p("    -forceselect <inputfile>    ").p("Force selection of the given input file in the files comparison method.").endl().endl();

      cout.p("    -backtrans <inputfile>      ").p("Use a Coding Sequences file to get a backtranslation for a given AA alignment").endl();
      cout.p("    -ignorestopcodon            ").p("Ignore stop codons in the input coding sequences").endl();
      cout.p("    -splitbystopcodon           ").p("Split input coding sequences up to first stop codon appearance").endl().endl();

      cout.p("    -matrix <inpufile>          ").p("Input file for user-defined similarity matrix (default is Blosum62).").endl();
      cout.p("    --alternative_matrix <name> ").p("Select an alternative similarity matrix already loaded. ").endl()
          .p("                                Only available 'degenerated_nt_identity'").endl().endl();

      cout.p("    -out <outputfile>           ").p("Output alignment in the same input format (default stdout). (default input format)").endl();
      cout.p("    -htmlout <outputfile>       ").p("Get a summary of trimal's work in an HTML file.").endl().endl();

      cout.p("    -keepheader                 ").p("Keep original sequence header including non-alphanumeric characters.").endl();
      cout.p("                                ").p("Only available for input FASTA format files. (future versions will extend this feature)").endl().endl();

      cout.p("    -nbrf                       ").p("Output file in NBRF/PIR format").endl();
      cout.p("    -mega                       ").p("Output file in MEGA format").endl();
      cout.p("    -nexus                      ").p("Output file in NEXUS format").endl();
      cout.p("    -clustal                    ").p("Output file in CLUSTAL format").endl().endl();

      cout.p("    -fasta                      ").p("Output file in FASTA format").endl();
      cout.p("    -fasta_m10                  ").p("Output file in FASTA format. Sequences name length up to 10 characters.").endl().endl();

      cout.p("    -phylip                     ").p("Output file in PHYLIP/PHYLIP4 format").endl();
      cout.p("    -phylip_m10                 ").p("Output file in PHYLIP/PHYLIP4 format. Sequences name length up to 10 characters.").endl();
      cout.p("    -phylip_paml                ").p("Output file in PHYLIP format compatible with PAML").endl();
      cout.p("    -phylip_paml_m10            ").p("Output file in PHYLIP format compatible with PAML. Sequences name length up to 10 characters.").endl();
      cout.p("    -phylip3.2                  ").p("Output file in PHYLIP3.2 format").endl();
      cout.p("    -phylip3.2_m10              ").p("Output file in PHYLIP3.2 format. Sequences name length up to 10 characters.").endl().endl();

      cout.p("    -complementary              ").p("Get the complementary alignment.").endl();
      cout.p("    -colnumbering               ").p("Get the relationship between the columns in the old and new alignment.").endl().endl();

      cout.p("    -selectcols { n,l,m-k }     ").p("Selection of columns to be removed from the alignment. Range: [0 - (Number of Columns - 1)]. (see User Guide).").endl();
      cout.p("    -selectseqs { n,l,m-k }     ").p("Selection of sequences to be removed from the alignment. Range: [0 - (Number of Sequences - 1)]. (see User Guide).").endl().endl();

      cout.p("    -gt -gapthreshold <n>       ").p("1 - (fraction of sequences with a gap allowed). Range: [0 - 1]").endl();
      cout.p("    -st -simthreshold <n>       ").p("Minimum average similarity allowed. Range: [0 - 1]").endl();
      cout.p("    -ct -conthreshold <n>       ").p("Minimum consistency value allowed.Range: [0 - 1]").endl();
      cout.p("    -cons <n>                   ").p("Minimum percentage of the positions in the original alignment to conserve. Range: [0 - 100]").endl().endl();

      cout.p("    -nogaps                     ").p("Remove all positions with gaps in the alignment.").endl();
      cout.p("    -noallgaps                  ").p("Remove columns composed only by gaps.").endl();
      cout.p("    -keepseqs                   ").p("Keep sequences even if they are composed only by gaps.").endl().endl();

      cout.p("    -gappyout                   ").p("Use automated selection on \"gappyout\" mode. This method only uses ")
          .p("information based on gaps' distribution. (see User Guide).").endl();
      cout.p("    -strict                     ").p("Use automated selection on \"strict\" mode. (see User Guide).").endl();
      cout.p("    -strictplus                 ").p("Use automated selection on \"strictplus\" mode. (see User Guide).").endl();
      cout.p("                               ").p("(Optimized for Neighbour Joining phylogenetic tree reconstruction).").endl().endl();

      cout.p("    -automated1                 ").p("Use a heuristic selection of the automatic method based on similarity statistics. ")
          .p("(see User Guide). (Optimized for Maximum Likelihood phylogenetic tree reconstruction).")
          .endl().endl();

      cout.p("    -terminalonly               ").p("Only columns out of internal boundaries (first and last column without gaps) are ").endl();
      cout.p("                                ").p("candidates to be trimmed depending on the selected method").endl();

      cout.p("    --set_boundaries { l,r }    ").p("Set manually left (l) and right (r) boundaries - only columns out of these boundaries are ").endl();
      cout.p("                                ").p("candidates to be trimmed depending on the selected method. Range: [0 - (Number of Columns - 1)]").endl();

      cout.p("    -block <n>                  ").p("Minimum column block size to be kept in the trimmed alignment. Available with manual")
          .p(" and automatic (gappyout) methods").endl().endl();

      cout.p("    -resoverlap                 ").p("Minimum overlap of a positions with other positions in the column to be considered a ")
          .p("\"good position\". Range: [0 - 1]. (see User Guide).").endl();
      cout.p("    -seqoverlap                 ").p("Minimum percentage of \"good positions\" that a sequence must have in order to be conserved. Range: [0 - 100]")
          .p("(see User Guide).").endl().endl();

      cout.p("    -clusters <n>               ").p("Get the most Nth representatives sequences from a given alignment. Range: [1 - (Number of sequences)]").endl();
      cout.p("    -maxidentity <n>            ").p("Get the representatives sequences for a given identity threshold. Range: [0 - 1].").endl().endl();

      cout.p("    -w <n>                      ").p("(half) Window size, score of position i is the average of the window (i - n) to (i + n).")
          .endl();
      cout.p("    -gw <n>                     ").p("(half) Window size only applies to statistics/methods based on Gaps.").endl();
      cout.p("    -sw <n>                     ").p("(half) Window size only applies to statistics/methods based on Similarity.").endl();
      cout.p("    -cw <n>                     ").p("(half) Window size only applies to statistics/methods based on Consistency.").endl().endl();

      cout.p("    -sgc                        ").p("Print gap scores for each column in the input alignment.").endl();
      cout.p("    -sgt                        ").p("Print accumulated gap scores for the input alignment.").endl();
      cout.p("    -ssc                        ").p("Print similarity scores for each column in the input alignment.").endl();
      cout.p("    -sst                        ").p("Print accumulated similarity scores for the input alignment.").endl();
      cout.p("    -sfc                        ").p("Print sum-of-pairs scores for each column from the selected alignment")
          .endl();
      cout.p("    -sft                        ").p("Print accumulated sum-of-pairs scores for the selected alignment")
          .endl();
      cout.p("    -sident                     ").p("Print identity scores matrix for all sequences in the input alignment. (see User Guide).")
          .endl();
      cout.p("    -soverlap                   ").p("Print overlap scores matrix for all sequences in the input alignment. (see User Guide).")
          .endl().endl();
    }
  }

  static void examples() throws IOException {
    try (CppOut cout = CppOut.toStdout()) {
      cout.p("Some Examples:").endl().endl();

      cout.p("1) Removes all positions in the alignment with gaps in 10% or more of").endl()
          .p("   the sequences, unless this leaves less than 60% of original alignment. ").endl()
          .p("   In such case, print the 60% best (with less gaps) positions.").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60").endl().endl();

      cout.p("2) As above but, the gap score is averaged over a window starting").endl()
          .p("   3 positions before and ending 3 positions after each column.").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -gt 0.9 -cons 60 -w 3").endl().endl();

      cout.p("3) Use an automatic method to decide optimal thresholds, based in the gap scores").endl()
          .p("   from input alignment. (see User Guide for details).").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -gappyout").endl().endl();

      cout.p("4) Use automatic methods to decide optimal thresholds, based on the combination ").endl()
          .p("   of gap and similarity scores. (see User Guide for details).").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -strictplus").endl().endl();

      cout.p("5) Use an heuristic to decide the optimal method for trimming the alignment. ").endl()
          .p("   (see User Guide for details).").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -automated1").endl().endl();

      cout.p("6) Use residues and sequences overlap thresholds to delete some sequences from the ").endl()
          .p("   alignemnt. (see User Guide for details).").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -resoverlap 0.8 -seqoverlap 75").endl().endl();

      cout.p("7) Selection of columns to be deleted from the alignment. The selection can ").endl()
          .p("   be a column number or a column number interval. Start from 0").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -selectcols { 0,2,3,10,45-60,68,70-78 }").endl().endl();

      cout.p("8) Get the complementary alignment from the alignment previously trimmed.").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -selectcols { 0,2,3,10,45-60,68,70-78 } -complementary").endl().endl();

      cout.p("9) Selection of sequences to be deleted from the alignment. Start in 0").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -selectseqs { 2,4,8-12 } ").endl().endl();

      cout.p("10) Select the 5 most representative sequences from the alignment").endl().endl();

      cout.p("   trimal -in <inputfile> -out <outputfile> -clusters 5 ").endl().endl();
    }
  }
}
