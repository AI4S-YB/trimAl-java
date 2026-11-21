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
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import trimal.stat.SimilarityMatrix;
import trimal.util.CppOut;
import trimal.util.Utils;

/** Faithful transcription of trimAl v1.5.1 source/statAl.cpp. */
public final class StatAl {

  private StatAl() {}

  public static void main(String[] args) {
    String[] argv = new String[args.length + 1];
    argv[0] = "statal";
    System.arraycopy(args, 0, argv, 1, args.length);
    try {
      int rc = run(argv);
      System.exit(rc < 0 ? rc + 256 : rc);
    } catch (IOException e) {
      System.exit(255);
    }
  }

  /** args includes the C argv[0] program-name slot. */
  static int run(String[] args) throws IOException {
    int argc = args.length;
    String inFile = null, forceFile = null, setAlignments = null, matrix = null;
    int windowSize = -1, gapWindow = -1, simWindow = -1, conWindow = -1,
        alternativeMatrix = -1;
    boolean statsGapsColumns = false, statsGapsDist = false, statsSimilColumns = false,
        statsSimilDist = false, statsSeqsIdent = false, statsColIdentGen = false,
        statsFileColumns = false, statsFileDist = false;
    int i = 1, numFiles = 0, maxResidues = 0, alignDataType = -1;
    SimilarityMatrix similMatrix = null;
    boolean appearErrors = false;
    float[] compareVect = null;
    List<String> alignmentPaths = null;

    if (argc == 1 || ((args[1].equals("-h") || args[1].equals("--help")) && argc == 2)) {
      showMenu();
      showExamples();
      return 0;
    }
    if ((args[1].equals("-v") || args[1].equals("--version")) && argc == 2) {
      try (CppOut cout = CppOut.toStdout()) {
        cout.endl().p("statAl v").p(TrimAl.VERSION).p(".rev").p(TrimAl.REVISION)
            .p(" build[").p(TrimAl.BUILD).p("]").endl().endl();
      }
      return 0;
    }

    Alignment origAlig = new Alignment();
    while (i < argc) {
      if ((args[i].equals("-i") || args[i].equals("-in")) && i + 1 != argc
          && inFile == null) {
        inFile = args[++i];
        if (!origAlig.loadAlignment(inFile)) {
          error("\nERROR: Alignment not loaded: \"" + inFile
              + "\" Check input aligmment content.\n\n");
          appearErrors = true;
        }
      } else if (args[i].equals("-forceselect") && i + 1 != argc && forceFile == null) {
        forceFile = args[++i];
        if (!origAlig.loadAlignment(forceFile)) {
          error("\nERROR: Alignment not loaded: \"" + forceFile
              + "\" Check input aligmment content.\n\n");
          appearErrors = true;
        }
      } else if (args[i].equals("-compareset") && i + 1 != argc && setAlignments == null) {
        setAlignments = args[++i];
      } else if (args[i].equals("-matrix") && i + 1 != argc && matrix == null) {
        matrix = args[++i];
        similMatrix = new SimilarityMatrix();
        if (!similMatrix.loadSimMatrix(matrix)) {
          error("\nERROR: Similarity Matrix not loaded: \"" + matrix
              + "\" Check input file content.\n\n");
          appearErrors = true;
        }
      } else if (args[i].equals("--alternative_matrix") && i + 1 != argc
          && alternativeMatrix == -1) {
        i++;
        if (args[i].equals("degenerated_nt_identity")) {
          alternativeMatrix = 1;
        } else {
          error("\nERROR: Alternative not recognized \"" + args[i] + "\"\n\n");
          appearErrors = true;
        }
      } else if (args[i].equals("-sgc") && !statsGapsColumns) {
        statsGapsColumns = true;
      } else if (args[i].equals("-sgt") && !statsGapsDist) {
        statsGapsDist = true;
      } else if (args[i].equals("-ssc") && !statsSimilColumns) {
        statsSimilColumns = true;
      } else if (args[i].equals("-sst") && !statsSimilDist) {
        statsSimilDist = true;
      } else if (args[i].equals("-sident") && !statsSeqsIdent) {
        statsSeqsIdent = true;
      } else if (args[i].equals("-scolidentt") && !statsColIdentGen) {
        statsColIdentGen = true;
      } else if (args[i].equals("-sfc") && !statsFileColumns) {
        statsFileColumns = true;
      } else if (args[i].equals("-sft") && !statsFileDist) {
        statsFileDist = true;
      } else if (args[i].equals("-w") && i + 1 != argc && windowSize == -1) {
        if (Utils.isNumber(args[i + 1])) {
          windowSize = Utils.atoi(args[++i]);
        } else {
          error("\nERROR: Window size value should be a number\n\n.");
          appearErrors = true;
        }
        if (windowSize < 1) {
          error("\nERROR: Windows size value should be equal or greater "
              + "than 1. Check your command-line parameter\n\n");
          appearErrors = true;
        }
      } else if (args[i].equals("-gw") && i + 1 != argc && gapWindow == -1) {
        if (Utils.isNumber(args[i + 1])) {
          gapWindow = Utils.atoi(args[++i]);
        } else {
          error("\nERROR: Gap Window size value should be a number\n\n.");
          appearErrors = true;
        }
        if (gapWindow < 1) {
          error("\nERROR: Gap Windows size should be equal or greater "
              + "than 1. Check your command-line parameter\n\n");
          appearErrors = true;
        }
      } else if (args[i].equals("-sw") && i + 1 != argc && simWindow == -1) {
        if (Utils.isNumber(args[i + 1])) {
          simWindow = Utils.atoi(args[++i]);
        } else {
          error("\nERROR: Similarity Window size should be a number\n\n.");
          appearErrors = true;
        }
        if (simWindow < 1) {
          error("\nERROR: Similarity Windows size should be equal or "
              + "gretaer than 1. Check your command-line parameter\n\n");
          appearErrors = true;
        }
      } else if (args[i].equals("-cw") && i + 1 != argc && conWindow == -1) {
        if (Utils.isNumber(args[i + 1])) {
          conWindow = Utils.atoi(args[++i]);
        } else {
          error("\nERROR: Consistency Window size should be a number\n\n");
          appearErrors = true;
        }
        if (conWindow < 1) {
          error("\nERROR: Consistency Windows size should be equal or "
              + "gretaer than 1. Check your command-line parameter\n\n");
          appearErrors = true;
        }
      } else {
        error("\nERROR: Parameter \"" + args[i] + "\" not valid\n\n.");
        appearErrors = true;
      }
      i++;
      if (appearErrors)
        break;
    }

    if (!appearErrors && inFile != null) {
      if (statsFileColumns || statsFileDist) {
        error("\nERROR: No stats about alignments comparison can be "
            + "obtained from Single Input Alignment\n\n");
        appearErrors = true;
      }
      // Both messages are emitted when both incompatibilities occur.
      if (forceFile != null || setAlignments != null) {
        error("\nERROR: Incompatibilities detected: Sigle Alignment "
            + "provided as well as Set of Alignments for being compared\n\n");
        appearErrors = true;
      }
    }
    if (!appearErrors && forceFile != null && setAlignments == null) {
      error("\nERROR: An alignmnet has been set as reference but it has"
          + "not been defined a set to compare with\n\n");
      appearErrors = true;
    }
    if (!appearErrors && (statsFileColumns || statsFileDist) && setAlignments == null) {
      error("\nERROR: A set of alignments should be provided to compute"
          + "stats about its comparison\n\n");
      appearErrors = true;
    }
    if (!appearErrors && windowSize != -1
        && (gapWindow != -1 || simWindow != -1 || conWindow != -1)) {
      error("\nERROR: General and specific windows size specified. Just"
          + " select one of them\n\n");
      appearErrors = true;
    }

    if (!appearErrors && setAlignments != null) {
      try {
        alignmentPaths = readGetlines(setAlignments);
        numFiles = alignmentPaths.size();
      } catch (IOException e) {
        error("\nERROR: Check input file \"" + setAlignments
            + "\" containing alignments paths\n\n");
        appearErrors = true;
      }
    }

    if (!appearErrors && numFiles != 0) {
      String[] filesToCompare = new String[numFiles];
      Alignment[] compAlig = new Alignment[numFiles];
      for (i = 0; i < numFiles && !appearErrors; i++) {
        String line = alignmentPaths.get(i);
        if (line.length() == 0)
          continue;
        filesToCompare[i] = line;
        compAlig[i] = new Alignment();
        if (!compAlig[i].loadAlignment(filesToCompare[i])) {
          error("\nERROR: Alignment not loaded: \"" + filesToCompare[i]
              + "\". Check input file content\n\n");
          appearErrors = true;
          break;
        }
        if (!compAlig[i].isFileAligned()) {
          error("\nERROR: Input alignment \"" + filesToCompare[i] + "\""
              + "should have sequences aligned for computing any stats from them.\n\n");
          appearErrors = true;
          break;
        }
        if (alignDataType == -1) {
          alignDataType = compAlig[i].getTypeAlignment();
        } else if (compAlig[i].getTypeAlignment() != alignDataType) {
          error("\nERROR: Check input alignment data-type. Alignment \""
              + filesToCompare[0] + "\" contains a different data-type different "
              + "than \"" + filesToCompare[i] + "\"\n\n");
          appearErrors = true;
          break;
        }
        compAlig[i].sequenMatrix();
        if (compAlig[i].getNumAminos() > maxResidues)
          maxResidues = compAlig[i].getNumAminos();
      }

      if (!appearErrors && forceFile == null) {
        compareVect = new float[maxResidues];
        int referFile = CompareFiles.algorithm(compAlig, filesToCompare, compareVect,
            numFiles, true);
        // C++ fails to check algorithm()=-1 before indexing filesToCompare. Avoid
        // Java's stack trace for that undefined native path; final checks report
        // the still-unloaded reference, as TrimAl's transcription also does.
        if (referFile >= 0)
          origAlig.loadAlignment(filesToCompare[referFile]);
      } else if (!appearErrors && forceFile != null) {
        compareVect = new float[origAlig.getNumAminos()];
        appearErrors = !CompareFiles.forceComparison(compAlig, numFiles, origAlig,
            compareVect);
      }
      conWindow = conWindow != -1 ? conWindow : windowSize;
      if (conWindow != -1 && !appearErrors)
        CompareFiles.applyWindow(origAlig.getNumAminos(), conWindow, compareVect);
    }

    boolean anyStats = statsGapsColumns || statsGapsDist || statsSimilColumns
        || statsSimilDist || statsSeqsIdent || statsColIdentGen || statsFileColumns
        || statsFileDist;
    if (!appearErrors && !origAlig.isFileAligned() && anyStats) {
      error("\nERROR: Reference alignment should aligned to compute any "
          + "stats using it\n\n");
      appearErrors = true;
    }
    if (!appearErrors && !anyStats) {
      error("\nERROR: An option has to be chosen.\n\n");
      appearErrors = true;
    }
    if (!appearErrors) {
      gapWindow = windowSize != -1 ? windowSize : gapWindow != -1 ? gapWindow : 0;
      simWindow = windowSize != -1 ? windowSize : simWindow != -1 ? simWindow : 0;
      origAlig.setWindowsSize(gapWindow, simWindow);
    }
    if (!appearErrors && (statsSimilColumns || statsSimilDist)) {
      if (matrix == null) {
        similMatrix = new SimilarityMatrix();
        alignDataType = origAlig.getTypeAlignment();
        if (alternativeMatrix != -1) {
          similMatrix.alternativeSimilarityMatrices(alternativeMatrix, alignDataType);
        } else if (alignDataType == Utils.AAType) {
          similMatrix.defaultAASimMatrix();
        } else if (alignDataType == Utils.DNAType || alignDataType == Utils.RNAType) {
          similMatrix.defaultNTSimMatrix();
        } else if (alignDataType == Utils.DNADeg || alignDataType == Utils.RNADeg) {
          similMatrix.defaultNTDegeneratedSimMatrix();
        }
      }
      if (!origAlig.setSimilarityMatrix(similMatrix)) {
        error("\nERROR: Impossible to associate similarity matrix to "
            + "selected alignment. Check input matrix file content\n\n");
        appearErrors = true;
      }
    }

    if (!appearErrors) {
      if (statsGapsColumns) {
        heading("## Gaps scores per column");
        origAlig.printStatisticsGapsColumns();
      }
      if (statsGapsDist) {
        heading("## Gaps scores distribution");
        origAlig.printStatisticsGapsTotal();
      }
      if (statsSimilColumns) {
        heading("## Similarity values per column");
        origAlig.printStatisticsConservationColumns();
      }
      if (statsSimilDist) {
        heading("## Similarity values distribution");
        origAlig.printStatisticsConservationTotal();
      }
      if (statsSeqsIdent) {
        heading("## Sequence Identity Values");
        origAlig.printSeqIdentity();
      }
      if (statsColIdentGen) {
        heading("## Columns Identity Descriptive Statistics");
        origAlig.printColumnsIdentity_DescriptiveStats();
      }
      if (statsFileColumns) {
        heading("## Consistency values per column for selected align");
        CompareFiles.printStatisticsFileColumns(origAlig.getNumAminos(), compareVect);
      }
      if (statsFileDist) {
        heading("## Consistency values distribution for selected align");
        CompareFiles.printStatisticsFileAcl(origAlig.getNumAminos(), compareVect);
      }
    }
    return appearErrors ? -1 : 0;
  }

  /** std::getline strips LF only, retaining CR and not adding an EOF empty line. */
  private static List<String> readGetlines(String filename) throws IOException {
    String text = new String(Files.readAllBytes(Path.of(filename)),
        StandardCharsets.ISO_8859_1);
    List<String> lines = new ArrayList<>();
    int start = 0;
    while (start < text.length()) {
      int end = text.indexOf('\n', start);
      if (end < 0) {
        lines.add(text.substring(start));
        break;
      }
      lines.add(text.substring(start, end));
      start = end + 1;
    }
    return lines;
  }

  private static void error(String message) throws IOException {
    try (CppOut cerr = CppOut.toStderr()) {
      cerr.p(message);
    }
  }

  private static void heading(String text) throws IOException {
    try (CppOut cout = CppOut.toStdout()) {
      cout.endl().p(text).endl();
    }
  }

  static void showMenu() throws IOException {
    try (CppOut cout = CppOut.toStdout()) {
      cout.endl().p("statAl v").p(TrimAl.VERSION).p(".rev").p(TrimAl.REVISION)
          .p(" build[").p(TrimAl.BUILD).p("]. ").p(TrimAl.AUTHORS).endl().endl()
          .p("statAl webpage: https://trimal.readthedocs.io\n\n")
          .p("This program is free software: you can redistribute it and/or modify")
          .p("\nit under the terms of the GNU General Public License as published by the")
          .p("\nFree Software Foundation, the last available version.\n\n")
          .p("Please cite:\n")
          .p("\ttrimAl: a tool for automated alignment trimming in large-scale ")
          .p("phylogenetic analyses.\n\tSalvador Capella-Gutierrez; Jose M. Silla-")
          .p("Martinez; Toni Gabaldon. Bioinformatics 2009, 25:1972-1973.\n\n")
          .p("Basic usage:\n\tstatal -in <inputfile> (options).\n\n")
          .p("Available options:\n\n")
          .p("    -h --help                Print this information and show some examples.\n")
          .p("    -v --version             Print the trimAl version.\n\n")
          .p("    -i -in <inputfile>       Input file in several formats ")
          .p("(clustal, fasta, nexus, phylip, etc).\n\n")
          .p("    -compareset <inputfile>     Input list of paths for the alignments to compare.\n")
          .p("    -forceselect <inputfile>    Force selection of a given file ")
          .p("as reference for being compare with others.\n\n")
          .p("    -matrix <inpufile>          Input file for user-defined similarity ")
          .p("matrix (default: Blosum62).\n")
          .p("    --alternative_matrix <name> Select an alternative similarity matrix already loaded. \n")
          .p("                                Only available 'degenerated_nt_identity'\n\n")
          .p("    -sgc                        Print gap score per column from input alignment.\n")
          .p("    -sgt                        Print accumulated gap scores distribution ")
          .p("from input alignment.\n\n")
          .p("    -ssc                        Print similarity score per column from input alignment.\n")
          .p("    -sst                        Print accumulated similarity scores ")
          .p("distribution for input alignment.\n\n")
          .p("    -sfc                        Print sum-of-pairs score per column ")
          .p("for the selected alignment\n")
          .p("    -sft                        Print accumulated sum-of-pairs scores ")
          .p("distribution for the selected alignment\n\n")
          .p("    -sident                     Print identity scores for sequences in the alignemnt.\n")
          .p("    -scolidentt                 Print general descriptive statistics ")
          .p("for column identity scores from input alignemnt.\n\n")
          .p("    -w <n>                      (half) Window size, score of position i ")
          .p("is the average of the window (i - n) to (i + n).\n")
          .p("    -gw <n>                     (half) Window size only applies to statistics based on Gaps.\n")
          .p("    -sw <n>                     (half) Window size only applies to statistics based on Similarity.\n")
          .p("    -cw <n>                     (half) Window size only applies to statistics based on Consistency.\n\n");
    }
  }

  static void showExamples() throws IOException {
    try (CppOut cout = CppOut.toStdout()) {
      cout.p("Some Examples:\n\n")
          .p("1) Get information about gaps distribution for input alignment\n")
          .p("   statal -in <inputfile> -sgt\n\n")
          .p("2) Get information about consistency score per column for the most ")
          .p("consistent input alignment, if more than one is provided\n")
          .p("   statal -in <inputfile> -sfc\n\n")
          .p("3) Get general descriptive statistics for columns identity\n")
          .p("   statal -in <inputfile> -scolidentt\n\n")
          .p("4) Change the windows size for computing similarity score per column\n")
          .p("   statal -in <inputfile> -sw 3 -ssc\n\n");
    }
  }
}
