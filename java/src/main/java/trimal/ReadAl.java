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

import trimal.util.CppOut;
import trimal.util.Utils;

/** Faithful transcription of trimAl v1.5.1 source/readAl.cpp. */
public final class ReadAl {

  private ReadAl() {}

  public static void main(String[] args) {
    String[] argv = new String[args.length + 1];
    argv[0] = "readal";
    System.arraycopy(args, 0, argv, 1, args.length);
    try {
      int rc = run(argv);
      // C main's -1 is status 255 on Linux, including when launched on Windows.
      System.exit(rc < 0 ? rc + 256 : rc);
    } catch (IOException e) {
      System.exit(255);
    }
  }

  /** args includes the C argv[0] program-name slot. */
  static int run(String[] args) throws IOException {
    int argc = args.length;
    Alignment inAlig = new Alignment();
    int outformat = -1;
    String infile = null, outfile = null;
    boolean errors = false, reverse = false, shortNames = false, format = false,
        type = false, info = false;

    if (argc == 1) {
      menu();
      return 0;
    }
    if (args[1].equals("-h") && argc == 2) {
      menu();
      return 0;
    }
    if (args[1].equals("--version") && argc == 2) {
      try (CppOut cout = CppOut.toStdout()) {
        cout.endl().p("readAl v").p(TrimAl.VERSION).p(".rev").p(TrimAl.REVISION)
            .p(" build[").p(TrimAl.BUILD).p("]").endl().endl();
      }
      return 0;
    }

    int i = 1;
    while (i < argc) {
      if (args[i].equals("-in") && i + 1 != argc && infile == null) {
        infile = args[++i];
        if (!inAlig.loadAlignment(infile)) {
          error("\nERROR: Alignment not loaded: \"" + infile
              + "\" Check the file's content.\n\n");
          errors = true;
        }
      } else if (args[i].equals("-out") && i + 1 != argc && outfile == null) {
        outfile = args[++i];
      } else if (args[i].equals("-format") && !format) {
        format = true;
      } else if (args[i].equals("-type") && !type) {
        type = true;
      } else if (args[i].equals("-info") && !info) {
        info = true;
      } else if (args[i].equals("-reverse") && !reverse) {
        reverse = true;
      } else if (args[i].equals("-clustal") && outformat == -1) {
        outformat = 1;
      } else if (args[i].equals("-fasta") && outformat == -1) {
        outformat = 8;
      } else if (args[i].equals("-fasta_m10") && outformat == -1) {
        outformat = 8;
        shortNames = true;
      } else if (args[i].equals("-nbrf") && outformat == -1) {
        outformat = 3;
      } else if (args[i].equals("-nexus") && outformat == -1) {
        outformat = 17;
      } else if (args[i].equals("-mega") && outformat == -1) {
        outformat = 21;
      } else if (args[i].equals("-phylip3.2") && outformat == -1) {
        outformat = 11;
      } else if (args[i].equals("-phylip3.2_m10") && outformat == -1) {
        outformat = 11;
        shortNames = true;
      } else if (args[i].equals("-phylip") && outformat == -1) {
        outformat = 12;
      } else if (args[i].equals("-phylip_m10") && outformat == -1) {
        outformat = 12;
        shortNames = true;
      } else if (args[i].equals("-phylip_paml") && outformat == -1) {
        outformat = 13;
      } else if (args[i].equals("-phylip_paml_m10") && outformat == -1) {
        outformat = 13;
        shortNames = true;
      } else if (args[i].equals("-html") && outformat == -1) {
        outformat = 100;
      } else if (args[i].equals("-onlyseqs") && outformat == -1) {
        outformat = 99;
      } else {
        error("\nERROR: Parameter \"" + args[i] + "\" not valid.\n\n");
        errors = true;
      }
      i++;
      if (errors)
        break;
    }

    if (infile == null && !errors) {
      error("\nERROR: An input file has to be defined.\n\n");
      errors = true;
    }
    if (outformat == -1 && !reverse && !format && !type && !info && !errors) {
      error("\nERROR: An option has to be chosen.\n\n");
      errors = true;
    }
    if (outfile == null && (outformat != -1 || reverse) && (format || type || info)
        && !errors) {
      error("\nERROR: Only one option can be selected: either an output "
          + "format or get information about input file when an output file is "
          + "not defined\n\n");
      errors = true;
    }
    if (outfile != null && outformat == -1 && !reverse && (format || type || info)
        && !errors) {
      error("\nERROR: An output file should not be provided when only "
          + "information about input alignment is requested\n\n");
      errors = true;
    }

    if (!errors) {
      if (format || type || info) {
        try (CppOut cout = CppOut.toStdout()) {
          cout.p("## Input filename\t'").p(infile).p("'").endl();
          if (format) {
            String alignFormat;
            switch (inAlig.getInputFormat()) {
              case 1: alignFormat = "clustal"; break;
              case 3: alignFormat = "nbrf/pir"; break;
              case 8: alignFormat = "fasta"; break;
              case 11: alignFormat = "phylip3.2"; break;
              case 12: alignFormat = "phylip"; break;
              case 17: alignFormat = "nexus"; break;
              case 21: alignFormat = "mega_interleaved"; break;
              case 22: alignFormat = "mega_sequential"; break;
              default: alignFormat = "unknown"; break;
            }
            cout.p("## Input file format\t").p(alignFormat).endl()
                .p("## Input file aligned\t").p(inAlig.isFileAligned() ? "YES" : "NO")
                .endl();
          }
          if (type) {
            String dataType;
            switch (inAlig.getTypeAlignment()) {
              case Utils.DNAType: dataType = "nucleotides:dna"; break;
              case Utils.DNADeg: dataType = "nucleotides:dna_degenerate_codes"; break;
              case Utils.RNAType: dataType = "nucleotides:rna"; break;
              case Utils.RNADeg: dataType = "nucleotides:rna_degenerate_codes"; break;
              case Utils.AAType: dataType = "amino-acids"; break;
              default: dataType = "unknown"; break;
            }
            cout.p("## Input file datatype\t").p(dataType).endl();
          }
          if (info)
            inAlig.printAlignmentInfo(cout);
        }
      }

      if (outfile != null || outformat != -1 || reverse || shortNames) {
        if (outformat != -1 || shortNames)
          inAlig.setOutputFormat(outformat, shortNames);
        if (reverse)
          inAlig.setReverse();
        if (outfile != null) {
          boolean saved;
          try {
            saved = inAlig.saveAlignment(outfile);
          } catch (IOException e) {
            // C++ ofstream reports an open/write failure as saveAlignment=false.
            saved = false;
          }
          if (!saved) {
            error("\nERROR: Impossible to generate OUTPUT file.\n\n");
            return -1;
          }
        } else {
          inAlig.printAlignment();
        }
      }
    }
    return errors ? -1 : 0;
  }

  private static void error(String message) throws IOException {
    try (CppOut cerr = CppOut.toStderr()) {
      cerr.p(message);
    }
  }

  static void menu() throws IOException {
    try (CppOut cout = CppOut.toStdout()) {
      cout.endl().p("readAl v").p(TrimAl.VERSION).p(".rev").p(TrimAl.REVISION)
          .p(" build[").p(TrimAl.BUILD).p("]. ").p(TrimAl.AUTHORS).endl().endl()
          .p("readAl webpage: https://trimal.readthedocs.io\n\n")
          .p("This program is free software: you can redistribute it and/or modify \n")
          .p("it under the terms of the GNU General Public License as published by \n")
          .p("the Free Software Foundation, the last available version.\n\n")
          .p("Basic usage\n")
          .p("\treadal -in <inputfile> -out <outputfile> [options].\n\n")
          .p("\t-h                   Show this information.\n")
          .p("\t--version            Show readAl version.\n\n")
          .p("\t-in <inputfile>      Input file in several formats.\n")
          .p("\t-out <outputfile>    Output file name (default STDOUT).\n\n")
          .p("\t-format              Print information about input file format ")
          .p("and if sequences are aligned or not.\n")
          .p("\t-type                Print information about biological ")
          .p("sequences datatype (e.g. nucleotides:dna, nucleotides:rna, aminoacids, etc)\n")
          .p("\t-info                Print information about sequences number, ")
          .p("average sequence length, max & min sequence length\n\n")
          .p("\t-onlyseqs            Generate output with only residues from input file\n\n")
          .p("\t-html                Output residues colored according their ")
          .p("physicochemical properties. HTML file.\n\n")
          .p("\t-reverse             Output the reverse of sequences in input file.\n\n")
          .p("\t-nbrf                Output file in NBRF/PIR format\n")
          .p("\t-mega                Output file in MEGA format\n")
          .p("\t-nexus               Output file in NEXUS format\n")
          .p("\t-clustal             Output file in CLUSTAL format\n\n")
          .p("\t-fasta               Output file in FASTA format\n")
          .p("\t-fasta_m10           Output file in FASTA format. Sequences ")
          .p("name up to 10 characters.\n\n")
          .p("\t-phylip              Output file in PHYLIP/PHYLIP4 format\n")
          .p("\t-phylip_m10          Output file in PHYLIP/PHYLIP4 format. ")
          .p("Sequences name up to 10 characters.\n")
          .p("\t-phylip_paml         Output file in PHYLIP format compatible with PAML\n")
          .p("\t-phylip_paml_m10     Output file in PHYLIP format compatible ")
          .p("with PAML. Sequences name up to 10 characters.\n")
          .p("\t-phylip3.2           Output file in PHYLIP3.2 format\n")
          .p("\t-phylip3.2_m10       Output file in PHYLIP3.2 format. Sequences")
          .p(" name up to 10 characters.\n\n");
    }
  }
}
