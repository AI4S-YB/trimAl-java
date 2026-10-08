package trimal;

import java.io.IOException;
import java.util.ArrayDeque;

import trimal.io.Tok;
import trimal.stat.SequencesMatrix;
import trimal.stat.SimilarityMatrix;
import trimal.stat.StatisticsConservation;
import trimal.stat.StatisticsGaps;
import trimal.util.CppFormat;
import trimal.util.CppOut;
import trimal.util.Parallel;
import trimal.util.Utils;

/**
 * Faithful transcription of trimAl v1.5.1 alignment.cpp + rwAlignment.cpp +
 * autAlignment.cpp. All strings are ISO-8859-1 (one char = one byte); writer output
 * goes through CppOut (std::ostream emulation). Quirks and UB-adjacent behavior are
 * preserved and marked with their C++ line references.
 */
public class Alignment {

  int sequenNumber;
  int residNumber;

  boolean isAligned;
  boolean reverse;

  boolean terminalGapOnly;
  int left_boundary;
  int right_boundary;

  int iformat;
  int oformat;
  boolean shortNames;

  boolean keepSequences;
  boolean keepHeader;

  String filename;
  String aligInfo;

  int dataType;

  int ghWindow;
  int shWindow;

  int blockSize;

  String[] sequences;
  String[] seqsName;
  String[] seqsInfo;

  StatisticsGaps sgaps;
  StatisticsConservation scons;
  SequencesMatrix seqMatrix;

  float[][] identities;
  float[][] overlaps;

  int[] residuesNumber;
  int[] saveResidues;
  int[] saveSequences;

  // Load-time accumulation buffers (C++ std::string::append is amortized O(1);
  // Java String += is quadratic, so loaders append here and fillMatrices materializes).
  private StringBuilder[] seqAcc;
  private StringBuilder[] nameAcc;
  private StringBuilder[] infoAcc;

  private void accInit(int n, boolean names, boolean infos) {
    seqAcc = new StringBuilder[n];
    nameAcc = names ? new StringBuilder[n] : null;
    infoAcc = infos ? new StringBuilder[n] : null;
    for (int z = 0; z < n; z++) {
      seqAcc[z] = new StringBuilder();
      if (names) nameAcc[z] = new StringBuilder();
      if (infos) infoAcc[z] = new StringBuilder();
    }
  }

  private void accMaterialize(int n) {
    if (seqAcc == null)
      return;
    for (int z = 0; z < n; z++) {
      sequences[z] = seqAcc[z].toString();
      if (nameAcc != null) seqsName[z] = nameAcc[z].toString();
      if (infoAcc != null) seqsInfo[z] = infoAcc[z].toString();
    }
    seqAcc = null;
    nameAcc = null;
    infoAcc = null;
  }

  public Alignment() {
    sequenNumber = 0;
    residNumber = 0;
    isAligned = false;
    reverse = false;
    terminalGapOnly = false;
    left_boundary = -1;
    right_boundary = -1;
    iformat = 0;
    oformat = 0;
    shortNames = false;
    keepSequences = false;
    keepHeader = false;
    filename = "";
    aligInfo = "";
    dataType = 0;
    ghWindow = 0;
    shWindow = 0;
    blockSize = 0;
    residuesNumber = null;
    saveResidues = null;
    saveSequences = null;
    sequences = null;
    seqsName = null;
    seqsInfo = null;
    sgaps = null;
    scons = null;
    seqMatrix = null;
    identities = null;
    overlaps = null;
  }

  /**
   * Copy constructor used by every clean* method (alignment.cpp:128-274). C++ quirk
   * preserved: identities is built from o_identities and then overwritten with NULL
   * at the end of the constructor (alignment.cpp:272), so the new object always
   * recomputes identities; overlaps are kept.
   */
  public Alignment(String o_filename, String o_aligInfo, String[] o_sequences, String[] o_seqsName,
      String[] o_seqsInfo, int o_sequenNumber, int o_residNumber, int o_iformat, int o_oformat,
      boolean o_shortNames, int o_dataType, boolean o_isAligned, boolean o_reverse,
      boolean o_terminalGapOnly, int o_left_boundary, int o_right_boundary,
      boolean o_keepSeqs, boolean o_keepHeader, int OldSequences, int OldResidues,
      int[] o_residuesNumber, int[] o_saveResidues, int[] o_saveSequences,
      int o_ghWindow, int o_shWindow, int o_blockSize, float[][] o_identities, float[][] o_overlaps) {

    int i, j;

    sequenNumber = o_sequenNumber;
    residNumber = o_residNumber;

    iformat = o_iformat;
    oformat = o_oformat;
    shortNames = o_shortNames;

    dataType = o_dataType;

    ghWindow = o_ghWindow;
    shWindow = o_shWindow;

    blockSize = o_blockSize;

    isAligned = o_isAligned;
    reverse = o_reverse;

    terminalGapOnly = o_terminalGapOnly;
    right_boundary = o_right_boundary;
    left_boundary = o_left_boundary;

    filename = o_filename;
    aligInfo = o_aligInfo;

    keepSequences = o_keepSeqs;
    keepHeader = o_keepHeader;

    sequences = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      sequences[i] = o_sequences[i];

    seqsName = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      seqsName[i] = o_seqsName[i];

    residuesNumber = new int[sequenNumber];
    if ((isAligned) || (o_residuesNumber != null)) {
      for (i = 0; i < sequenNumber; i++)
        residuesNumber[i] = residNumber;
    } else {
      for (i = 0; i < sequenNumber; i++)
        residuesNumber[i] = o_residuesNumber[i];
    }

    if (o_seqsInfo != null) {
      seqsInfo = new String[sequenNumber];
      for (i = 0; i < sequenNumber; i++)
        seqsInfo[i] = o_seqsInfo[i];
    } else seqsInfo = null;

    saveResidues = null;
    if (o_saveResidues != null) {
      saveResidues = new int[residNumber];
      for (i = 0, j = 0; i < OldResidues; i++)
        if (o_saveResidues[i] != -1) {
          saveResidues[j] = o_saveResidues[i];
          j++;
        }
    }

    saveSequences = null;
    if (o_saveSequences != null) {
      saveSequences = new int[sequenNumber];
      for (i = 0, j = 0; i < OldSequences; i++)
        if (o_saveSequences[i] != -1) {
          saveSequences[j] = o_saveSequences[i];
          j++;
        }
    }

    // C++ builds identities from o_identities here and then NULLs it (line 272).
    identities = null;

    overlaps = null;
    if (o_overlaps != null) {
      overlaps = new float[sequenNumber][];
      for (i = 0, j = 0; i < OldSequences; i++) {
        if (o_saveSequences[i] != -1) {
          float[] row = new float[sequenNumber];
          for (int k = 0, ll = 0; k < OldSequences; k++) {
            if (o_saveSequences[k] != -1) {
              row[ll] = o_overlaps[i][k];
              ll++;
            }
          }
          overlaps[j] = row;
          j++;
        }
      }
    }

    sgaps = null;
    scons = null;
    seqMatrix = null;
  }

  // ------------------------------------------------------------------ io infra

  /**
   * utils::readLine emulation over an in-memory buffer: one artificial leading space
   * is prepended (the C++ for-loop writes its initial c=' ' first), '\r' is dropped,
   * then all leading spaces are stripped followed by all leading tabs; empty -> null.
   */
  static final class Cursor {
    final String content;
    int pos = 0;
    boolean eofFlag = false;

    Cursor(byte[] bytes) {
      // ISO-8859-1 maps every byte 0-255 to the identical char: this is a plain copy,
      // so it matches the previous per-byte loop while running at memcpy speed.
      content = new String(bytes, java.nio.charset.StandardCharsets.ISO_8859_1);
    }

    /** ifstream::read(&c, 1): at end leaves c unchanged and raises eof. */
    int read() {
      if (pos < content.length()) return content.charAt(pos++);
      eofFlag = true;
      return -1;
    }

    boolean atEof() {
      // ifstream::eof() is set by a failed read, not by consuming the final byte.
      return eofFlag;
    }

    String readLine() {
      if (eofFlag) return null;
      StringBuilder nline = new StringBuilder();
      nline.append(' '); // C++ initial c=' ' append
      while (true) {
        int c = read();
        if (c == -1 || c == '\n') break;
        if (c != '\r') nline.append((char) c);
      }
      int s = 0;
      while (s < nline.length() && nline.charAt(s) == ' ') s++;
      int t = s;
      while (t < nline.length() && nline.charAt(t) == '\t') t++;
      if (t == nline.length()) return null;
      return nline.substring(t);
    }
  }

  static Cursor open(String file) throws IOException {
    byte[] raw = java.nio.file.Files.readAllBytes(java.nio.file.Path.of(file));
    return new Cursor(raw);
  }

  static CppOut cerr() throws IOException {
    return CppOut.toStderr();
  }

  private static boolean isalpha(char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
  }

  private static boolean ispunct(char c) {
    return (c >= 33 && c <= 47) || (c >= 58 && c <= 64) || (c >= 91 && c <= 96) || (c >= 123 && c <= 126);
  }

  private static char toupper(char c) {
    return c >= 'a' && c <= 'z' ? (char) (c - ('a' - 'A')) : c;
  }

  private static String toupper(String s) {
    StringBuilder upper = new StringBuilder(s.length());
    for (int i = 0; i < s.length(); i++) upper.append(toupper(s.charAt(i)));
    return upper.toString();
  }

  /** C++ string::substr(pos, len) with clamping. */
  static String sub(String s, int pos, int len) {
    if (pos >= s.length()) return "";
    int end = (len >= s.length() - pos) ? s.length() : pos + len;
    return s.substring(pos, end);
  }

  /** appends the whole C-string of t to sb */
  private static void append(StringBuilder sb, Tok t) {
    if (t != null) sb.append(t.toString());
  }

  public boolean fillMatrices(boolean aligned) throws IOException {
    int i, j;

    accMaterialize(sequenNumber);

    residuesNumber = new int[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      residuesNumber[i] = sequences[i].length();

    for (i = 0; i < sequenNumber; i++)
      for (j = 0; j < residuesNumber[i]; j++) {
        char c = sequences[i].charAt(j);
        if (!(isalpha(c)) && !(ispunct(c))) {
          try (CppOut err = cerr()) {
            err.endl().p("ERROR: The sequence \"").p(seqsName[i]).p("\" has an ")
                .p("unknown (").p(c).p(") character.").endl();
          }
          return false;
        }
      }

    for (i = 1; i < sequenNumber; i++)
      if (residuesNumber[i] != residuesNumber[i - 1])
        break;
    isAligned = (i != sequenNumber) ? false : true;

    if (aligned && !isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences should be aligned (all with same length) ")
            .p("and there are not. Check your input alignment").endl();
      }
      return false;
    }

    if (residNumber == 0)
      residNumber = residuesNumber[0];

    for (i = 0; (i < sequenNumber) && (aligned); i++) {
      if (residuesNumber[i] != residNumber) {
        try (CppOut err = cerr()) {
          err.endl().p("ERROR: The sequence \"").p(seqsName[i]).p("\" (")
              .p(residuesNumber[i]).p(") does not have the same number of residues ")
              .p("fixed by the alignment (").p(residNumber).p(").").endl();
        }
        return false;
      }
    }

    if ((aligned) || (isAligned)) {
      saveResidues = new int[residNumber];
      for (i = 0; i < residNumber; i++)
        saveResidues[i] = i;

      saveSequences = new int[sequenNumber];
      for (i = 0; i < sequenNumber; i++)
        saveSequences[i] = i;
    }

    return true;
  }

  // ------------------------------------------------------------------ detection

  public int formatInputAlignment(String alignmentFile) throws IOException {
    int format = 0, blocks = 0;

    Cursor file;
    try {
      file = open(alignmentFile);
    } catch (IOException e) {
      return 0;
    }

    String line = null;
    do {
      line = file.readLine();
    } while ((line == null) && (!file.atEof()));

    if (file.atEof() || line == null)
      return 0;

    Tok lw = new Tok(line);
    Tok firstWord = Tok.strtok(lw, Utils.OTHDELIMITERS);

    if (firstWord == null)
      return 0;

    String fw = firstWord.toString();

    if (fw.equals("CLUSTAL") || fw.equals("clustal"))
      format = 1;
    else if (firstWord.at(0) == '>' && firstWord.at(3) == ';')
      format = 3;
    else if (firstWord.at(0) == '>')
      format = 8;
    else if (fw.equals("#NEXUS") || fw.equals("#nexus"))
      format = 17;
    else if (fw.equals("#MEGA") || fw.equals("#mega")) {
      // byte-wise scan counting '#' line blocks (rwAlignment.cpp:191-202)
      int c;
      blocks = 0;
      do {
        c = file.read();
      } while ((c != '#') && (!file.eofFlag));

      do {
        while ((c != '\n') && (!file.eofFlag))
          c = file.read();
        c = file.read();
        if (c == '#')
          blocks++;
      } while ((c != '\n') && (!file.eofFlag));

      format = (blocks == 0) ? 22 : 21;
    } else {
      sequenNumber = Utils.atoi(fw);
      Tok t2 = Tok.strtok(lw, Utils.DELIMITERS);
      residNumber = 0;
      if (t2 != null)
        residNumber = Utils.atoi(t2.toString());

      if ((sequenNumber == 1) && (residNumber != 0))
        format = 12;
      else if ((sequenNumber != 0) && (residNumber != 0)) {
        blocks = 0;

        do {
          line = file.readLine();
        } while ((line == null) && (!file.atEof()));

        if (file.atEof() || line == null)
          return 0;

        Tok lw2 = new Tok(line);
        Tok w = Tok.strtok(lw2, Utils.DELIMITERS);
        while (w != null) {
          blocks++;
          w = Tok.strtok(lw2, Utils.DELIMITERS);
        }

        do {
          line = file.readLine();
        } while ((line == null) && (!file.atEof()));

        if (line == null)
          return 0;

        lw2 = new Tok(line);
        w = Tok.strtok(lw2, Utils.DELIMITERS);
        while (w != null) {
          blocks--;
          w = Tok.strtok(lw2, Utils.DELIMITERS);
        }

        if (file.atEof())
          return 0;

        format = (blocks == 0) ? 12 : 11;
      }
    }

    return format;
  }

  public boolean loadAlignment(String alignmentFile) {
    // C++ ifstream::open failure yields checkFile()==false -> format 0 -> false;
    // an unreadable file maps to the same observable path.
    try {
      iformat = formatInputAlignment(alignmentFile);
    } catch (IOException e) {
      iformat = 0;
    }
    oformat = iformat;

    try {
      switch (iformat) {
        case 1:
          return loadClustalAlignment(alignmentFile);
        case 3:
          return loadNBRF_PirAlignment(alignmentFile);
        case 8:
          return loadFastaAlignment(alignmentFile);
        case 11:
          return loadPhylip3_2Alignment(alignmentFile);
        case 12:
          return loadPhylipAlignment(alignmentFile);
        case 17:
          return loadNexusAlignment(alignmentFile);
        case 21:
          return loadMegaInterleavedAlignment(alignmentFile);
        case 22:
          return loadMegaNonInterleavedAlignment(alignmentFile);
        default:
          return false;
      }
    } catch (IOException e) {
      return false;
    }
  }

  // ------------------------------------------------------------------ loaders

  /** PHYLIP interleaved parser (rwAlignment.cpp:276-381). */
  public boolean loadPhylipAlignment(String alignmentFile) throws IOException {
    Cursor file = open(alignmentFile);

    filename = "!Title " + alignmentFile + ";";

    String line = null;
    do {
      line = file.readLine();
    } while ((line == null) && (!file.atEof()));

    if (file.atEof() || line == null)
      return false;

    Tok lw = new Tok(line);
    Tok str = Tok.strtok(lw, Utils.DELIMITERS);
    sequenNumber = 0;
    if (str != null)
      sequenNumber = Utils.atoi(str.toString());

    str = Tok.strtok(lw, Utils.DELIMITERS);
    residNumber = 0;
    if (str != null)
      residNumber = Utils.atoi(str.toString());

    if ((sequenNumber == 0) || (residNumber == 0))
      return false;

    sequences = new String[sequenNumber];
    seqsName = new String[sequenNumber];
    accInit(sequenNumber, true, false);
    for (int z = 0; z < sequenNumber; z++) {
      sequences[z] = "";
      seqsName[z] = "";
    }

    int i = 0;
    while ((i < sequenNumber) && (!file.atEof())) {
      line = file.readLine();
      if (line == null)
        continue;

      lw = new Tok(line);
      str = Tok.strtok(lw, Utils.DELIMITERS);
      if (str != null) {
        nameAcc[i].append(str.toString());
        str = Tok.strtok(lw, Utils.DELIMITERS);
        while (str != null) {
          seqAcc[i].append(str.toString());
          str = Tok.strtok(lw, Utils.DELIMITERS);
        }
      }
      i++;
    }

    while (!file.atEof()) {
      i = 0;
      while ((i < sequenNumber) && (!file.atEof())) {
        line = file.readLine();
        if (line == null)
          continue;

        lw = new Tok(line);
        str = Tok.strtok(lw, Utils.DELIMITERS);
        while (str != null) {
          seqAcc[i].append(str.toString());
          str = Tok.strtok(lw, Utils.DELIMITERS);
        }
        i++;
      }
    }

    return fillMatrices(true);
  }

  /** PHYLIP 3.2 sequential parser (rwAlignment.cpp:383-496). */
  public boolean loadPhylip3_2Alignment(String alignmentFile) throws IOException {
    Cursor file = open(alignmentFile);

    filename = "!Title " + alignmentFile + ";";

    String line = null;
    do {
      line = file.readLine();
    } while ((line == null) && (!file.atEof()));

    if (file.atEof() || line == null)
      return false;

    Tok lw = new Tok(line);
    Tok str = Tok.strtok(lw, Utils.DELIMITERS);
    sequenNumber = 0;
    if (str != null)
      sequenNumber = Utils.atoi(str.toString());

    str = Tok.strtok(lw, Utils.DELIMITERS);
    residNumber = 0;
    if (str != null)
      residNumber = Utils.atoi(str.toString());

    if ((sequenNumber == 0) || (residNumber == 0))
      return false;

    sequences = new String[sequenNumber];
    seqsName = new String[sequenNumber];
    accInit(sequenNumber, true, false);
    for (int z = 0; z < sequenNumber; z++) {
      sequences[z] = "";
      seqsName[z] = "";
    }

    int i = 0;
    int blocksFirstLine = 0;
    int firstLine = 1;

    do {
      line = file.readLine();
      if (line == null)
        continue;

      lw = new Tok(line);
      str = Tok.strtok(lw, Utils.DELIMITERS);
      if (firstLine != 0) {
        if (str != null) nameAcc[i].append(str.toString());
        str = Tok.strtok(lw, Utils.OTHDELIMITERS);
        firstLine = 1;
      }

      while (str != null) {
        seqAcc[i].append(str.toString());
        str = Tok.strtok(lw, Utils.OTHDELIMITERS);
        if (firstLine != 0)
          firstLine += 1;
      }

      if ((blocksFirstLine == 0) && (firstLine != 0))
        blocksFirstLine = firstLine;

      if ((firstLine != 0) && (firstLine != blocksFirstLine)) {
        if (i > 0) {
          seqAcc[i - 1].append(nameAcc[i]);
          nameAcc[i].setLength(0);
          seqAcc[i - 1].append(seqAcc[i]);
          seqAcc[i].setLength(0);
          i--;
        }
      }

      firstLine = 0;
      if (seqAcc[i].length() == residNumber) {
        firstLine = 1;
        i++;
      }
    } while (!file.atEof());

    return fillMatrices(true);
  }

  /** CLUSTAL parser (rwAlignment.cpp:498-671). */
  public boolean loadClustalAlignment(String alignmentFile) throws IOException {
    Cursor file = open(alignmentFile);

    filename = "!Title " + alignmentFile + ";";

    String line = null;
    do {
      line = file.readLine();
    } while ((line == null) && (!file.atEof()));

    if (file.atEof() || line == null)
      return false;

    // count sequences from the first block
    while (!file.atEof()) {
      line = file.readLine();
      if (line != null)
        break;
    }

    sequenNumber = 0;
    while (!file.atEof()) {
      if (line == null)
        break;

      int seqLength = line.length();
      int pos;
      for (pos = 0; pos < seqLength; pos++)
        if (isalpha(line.charAt(pos)) || (line.charAt(pos) == '-'))
          break;

      if (pos == seqLength)
        break;
      sequenNumber++;

      line = file.readLine();
    }

    // second pass from the beginning (C++: file.clear(); file.seekg(0);)
    file = open(alignmentFile);
    line = file.readLine();
    if (line == null)
      return false;
    aligInfo = line;

    while (!file.atEof()) {
      line = file.readLine();
      if (line != null)
        break;
    }

    seqsName = new String[sequenNumber];
    sequences = new String[sequenNumber];
    accInit(sequenNumber, true, false);
    for (int z = 0; z < sequenNumber; z++) {
      sequences[z] = "";
      seqsName[z] = "";
    }

    int i = 0;
    boolean firstBlock = true;

    while (!file.atEof()) {
      if (line == null) {
        if (i == 0)
          firstBlock = false;
        line = file.readLine();
        continue;
      }

      int seqLength = line.length();
      int pos;
      for (pos = 0; pos < seqLength; pos++)
        if (isalpha(line.charAt(pos)) || (line.charAt(pos) == '-'))
          break;

      if (pos == seqLength) {
        firstBlock = false;
        line = file.readLine();
        continue;
      }

      Tok lw = new Tok(line);
      Tok str = Tok.strtok(lw, Utils.OTHDELIMITERS);
      if (str != null) {
        if (firstBlock)
          nameAcc[i].append(str.toString());
        str = Tok.strtok(lw, Utils.OTHDELIMITERS);
        if (str != null)
          seqAcc[i].append(str.toString());
        i = (i + 1) % sequenNumber;
      }

      line = file.readLine();
    }

    return fillMatrices(true);
  }

  /** FASTA parser (rwAlignment.cpp:673-770). */
  public boolean loadFastaAlignment(String alignmentFile) throws IOException {
    Cursor file = open(alignmentFile);

    filename = "!Title " + alignmentFile + ";";

    String line = null;
    sequenNumber = 0;
    while (!file.atEof()) {
      line = file.readLine();
      if (line == null)
        continue;

      Tok lw = new Tok(line);
      Tok str = Tok.strtok(lw, Utils.DELIMITERS);
      if (str == null)
        continue;

      if (str.at(0) == '>')
        sequenNumber++;
    }

    // second pass (C++: file.clear(); file.seekg(0);)
    file = open(alignmentFile);
    seqsName = new String[sequenNumber];
    sequences = new String[sequenNumber];
    seqsInfo = new String[sequenNumber];
    accInit(sequenNumber, true, true);
    for (int z = 0; z < sequenNumber; z++) {
      sequences[z] = "";
      seqsName[z] = "";
      seqsInfo[z] = "";
    }

    int i = -1;
    while ((i < sequenNumber) && (!file.atEof())) {
      line = file.readLine();
      if (line == null)
        continue;

      if (line.charAt(0) == '>')
        infoAcc[i + 1].append(line, 1, line.length());

      Tok lw = new Tok(line);
      Tok str = Tok.strtok(lw, Utils.OTHDELIMITERS);
      if (str == null)
        continue;

      if (str.at(0) == '>') {
        // do { str = str + 1; } while(strlen(str) == 0);
        str = str.advance();
        while (str.isEmpty() && str.start + 1 < str.buf.length)
          str = str.advance();
        i++;
        if (i >= sequenNumber) break;
        nameAcc[i].append(str.toString());
        continue;
      }

      while (str != null) {
        seqAcc[i].append(str.toString());
        str = Tok.strtok(lw, Utils.DELIMITERS);
      }
    }

    return fillMatrices(false);
  }

  /** NEXUS parser (rwAlignment.cpp:772-913). */
  public boolean loadNexusAlignment(String alignmentFile) throws IOException {
    Cursor file = open(alignmentFile);

    filename = "!Title " + alignmentFile + ";";

    int i, pos, state = 0;
    boolean firstBlock;
    String line = null;
    boolean sawMatrix = false;
    String lastToken = null;

    do {
      line = file.readLine();
      if (line == null)
        continue;

      Tok lw = new Tok(line);
      Tok str = Tok.strtok(lw, Utils.DELIMITERS);
      if (str == null)
        continue;

      String up = toupper(str.toString());
      lastToken = up;

      if (up.equals("BEGIN"))
        state = 1;
      else if (up.equals("MATRIX")) {
        sawMatrix = true;
        break;
      } else if (up.equals("FORMAT")) {
        str = Tok.strtok(lw, Utils.DELIMITERS);
        while (str != null) {
          aligInfo += str.toString();
          aligInfo += " ";
          str = Tok.strtok(lw, Utils.DELIMITERS);
        }
      } else if (up.equals("DIMENSIONS") && state != 0) {
        str = Tok.strtok(lw, Utils.DELIMITERS);
        Tok frag = Tok.strtok(lw, Utils.DELIMITERS);
        // str = strtok(str, "=;"); sequenNumber = atoi(strtok(NULL, "=;"));
        String s1 = str != null ? str.toString() : "";
        String[] parts1 = splitOn(s1, "=;");
        sequenNumber = parts1.length > 1 ? Utils.atoi(parts1[1]) : 0;
        String f1 = frag != null ? frag.toString() : "";
        String[] parts2 = splitOn(f1, "=;");
        residNumber = parts2.length > 1 ? Utils.atoi(parts2[1]) : 0;
      }
    } while (!file.atEof());

    if ((!sawMatrix) || (sequenNumber == 0) || (residNumber == 0))
      return false;

    seqsName = new String[sequenNumber];
    sequences = new String[sequenNumber];
    accInit(sequenNumber, true, false);
    for (int z = 0; z < sequenNumber; z++) {
      sequences[z] = "";
      seqsName[z] = "";
    }

    pos = 0;
    state = 0;
    firstBlock = true;

    while (!file.atEof()) {
      line = file.readLine();
      if (line == null)
        continue;

      // comment skip: any line containing '[' or ']' is skipped entirely
      int idx;
      for (idx = 0; idx < line.length(); idx++) {
        if (line.charAt(idx) == '[')
          state = 1;
        else if ((line.charAt(idx) == ']') && (state != 0)) {
          state = 0;
          break;
        }
      }
      if ((state != 0) || (state == 0 && idx != line.length()))
        continue;

      if (line.startsWith("end;") || line.startsWith("END;"))
        break;

      Tok lw = new Tok(line);
      Tok str = Tok.strtok(lw, Utils.OTH2DELIMITERS);
      if (str == null)
        continue;

      if (firstBlock)
        nameAcc[pos].append(str.toString());

      str = Tok.strtok(lw, Utils.OTH2DELIMITERS);
      while (str != null) {
        seqAcc[pos].append(str.toString());
        str = Tok.strtok(lw, Utils.OTH2DELIMITERS);
      }

      pos = (pos + 1) % sequenNumber;
      if (pos == 0)
        firstBlock = false;
    }

    return fillMatrices(true);
  }

  /** "a=b" style split honoring C strtok(str, "=;") + strtok(NULL, "=;") on one token. */
  private static String[] splitOn(String s, String delims) {
    // C: strtok(str, "=;") returns the prefix up to the first delimiter;
    // strtok(NULL) continues after it. On "NTAX=6": "NTAX", "6".
    java.util.ArrayList<String> parts = new java.util.ArrayList<>();
    int start = 0;
    // skip leading delimiters
    while (start < s.length() && delims.indexOf(s.charAt(start)) >= 0) start++;
    while (start < s.length()) {
      int end = start;
      while (end < s.length() && delims.indexOf(s.charAt(end)) < 0) end++;
      parts.add(s.substring(start, end));
      if (end >= s.length()) break;
      start = end + 1;
      while (start < s.length() && delims.indexOf(s.charAt(start)) >= 0) start++;
    }
    return parts.toArray(new String[0]);
  }

  /** MEGA sequential parser (rwAlignment.cpp:915-1112). */
  public boolean loadMegaNonInterleavedAlignment(String alignmentFile) throws IOException {
    Cursor file = open(alignmentFile);

    filename = "!Title " + alignmentFile + ";";

    String line = null;
    do {
      line = file.readLine();
    } while ((line == null) && (!file.atEof()));

    if (file.atEof() || line == null)
      return false;

    // gather TITLE/FORMAT information
    while (!file.atEof()) {
      line = file.readLine();
      if (line == null)
        continue;

      if (line.startsWith("#"))
        break;

      String frag = line;
      // str = strtok(frag, "!: "); then uppercase each char
      Tok lw = new Tok(frag);
      Tok str = Tok.strtok(lw, "!: ");
      if (str == null) continue;
      String up = toupper(str.toString());

      if (up.equals("TITLE")) {
        filename = "";
        if (!line.startsWith("!"))
          filename += "!";
        filename += line;
      } else if (up.equals("FORMAT"))
        aligInfo += line;
    }

    // count sequences (lines starting with '#')
    do {
      if (line == null) {
        line = file.readLine();
        continue;
      }
      if (line.startsWith("#"))
        sequenNumber++;

      line = file.readLine();
    } while (!file.atEof());

    // restart
    file = open(alignmentFile);
    seqsName = new String[sequenNumber];
    sequences = new String[sequenNumber];
    accInit(sequenNumber, true, false);
    for (int z = 0; z < sequenNumber; z++) {
      sequences[z] = "";
      seqsName[z] = "";
    }

    line = file.readLine(); // skip first line

    while (!file.atEof()) {
      line = file.readLine();
      if (line == null)
        continue;
      if (line.startsWith("#"))
        break;
    }

    int i = -1;

    while (!file.atEof()) {
      if (line == null) {
        line = file.readLine();
        continue;
      }

      if (line.startsWith("!")) {
        line = file.readLine();
        continue;
      }

      String frag = Utils.trimLine(line);
      if (frag == null) {
        line = file.readLine();
        continue;
      }

      Tok lw = new Tok(frag);
      Tok str = Tok.strtok(lw, " #\n");

      if (line.startsWith("#")) {
        i += 1;
        if (str != null) nameAcc[i].append(str.toString());
        str = Tok.strtok(lw, " #\n");
      }

      while (str != null) {
        seqAcc[i].append(str.toString());
        str = Tok.strtok(lw, " \n");
      }

      line = file.readLine();
    }

    return fillMatrices(true);
  }

  /** MEGA interleaved parser (rwAlignment.cpp:1114-1295). */
  public boolean loadMegaInterleavedAlignment(String alignmentFile) throws IOException {
    Cursor file = open(alignmentFile);

    filename = "!Title " + alignmentFile + ";";

    String line = null;
    do {
      line = file.readLine();
    } while ((line == null) && (!file.atEof()));

    if (file.atEof() || line == null)
      return false;

    while (!file.atEof()) {
      line = file.readLine();
      if (line == null)
        continue;

      if (line.startsWith("#"))
        break;

      Tok lw = new Tok(line);
      Tok str = Tok.strtok(lw, "!: ");
      if (str == null) continue;
      String up = toupper(str.toString());

      if (up.equals("TITLE")) {
        filename = "";
        if (!line.startsWith("!"))
          filename += "!";
        filename += line;
      } else if (up.equals("FORMAT"))
        aligInfo += line;
    }

    // count sequences in the first block
    while (!file.atEof()) {
      if (line != null && line.startsWith("#"))
        sequenNumber++;

      line = file.readLine();

      if (line == null)
        break;
    }

    // restart
    file = open(alignmentFile);
    seqsName = new String[sequenNumber];
    sequences = new String[sequenNumber];
    accInit(sequenNumber, true, false);
    for (int z = 0; z < sequenNumber; z++) {
      sequences[z] = "";
      seqsName[z] = "";
    }

    line = file.readLine(); // skip first line

    while (!file.atEof()) {
      line = file.readLine();
      if (line != null && line.startsWith("#"))
        break;
    }

    int i = 0;
    boolean firstBlock = true;

    while (!file.atEof()) {
      if (line == null) {
        line = file.readLine();
        continue;
      }

      if (line.startsWith("!")) {
        line = file.readLine();
        continue;
      }

      String frag = Utils.trimLine(line);
      Tok lw = frag != null ? new Tok(frag) : null;
      Tok str = lw != null ? Tok.strtok(lw, " #\n") : null;

      if (str == null) {
        line = file.readLine();
        continue;
      }

      if (firstBlock)
        nameAcc[i].append(str.toString());

      str = Tok.strtok(lw, " \n");
      while (str != null) {
        seqAcc[i].append(str.toString());
        str = Tok.strtok(lw, " \n");
      }

      line = file.readLine();

      i = (i + 1) % sequenNumber;
      if (i == 0)
        firstBlock = false;
    }

    return fillMatrices(true);
  }

  /** NBRF/PIR parser (rwAlignment.cpp:1297-1419). */
  public boolean loadNBRF_PirAlignment(String alignmentFile) throws IOException {
    Cursor file = open(alignmentFile);

    filename = "!Title " + alignmentFile + ";";

    String line = null;
    sequenNumber = 0;
    while (!file.atEof()) {
      line = file.readLine();
      if (line == null)
        continue;

      Tok lw = new Tok(line);
      Tok str = Tok.strtok(lw, Utils.DELIMITERS);
      if (str == null)
        continue;

      if (str.at(0) == '>')
        sequenNumber++;
    }

    // restart
    file = open(alignmentFile);
    sequences = new String[sequenNumber];
    seqsName = new String[sequenNumber];
    seqsInfo = new String[sequenNumber];
    accInit(sequenNumber, true, true);
    for (int z = 0; z < sequenNumber; z++) {
      sequences[z] = "";
      seqsName[z] = "";
      seqsInfo[z] = "";
    }

    boolean seqIdLine = true;
    boolean seqLines = false;
    int i = -1;

    while (!file.atEof()) {
      line = file.readLine();
      if (line == null)
        continue;

      if ((line.charAt(0) == '>') && (line.length() > 3 && line.charAt(3) == ';') && (seqIdLine)) {
        seqIdLine = false;
        i += 1;

        Tok lw = new Tok(line);
        Tok str = Tok.strtok(lw, ">;");
        if (str != null) infoAcc[i].append(str.toString());
        str = Tok.strtok(lw, ">;");
        if (str != null) nameAcc[i].append(str.toString());
      } else if ((!seqIdLine) && (!seqLines)) {
        seqLines = true;
        infoAcc[i].append(line);
      } else if (seqLines) {
        if (line.charAt(line.length() - 1) == '*') {
          seqLines = false;
          seqIdLine = true;
        }

        Tok lw = new Tok(line);
        Tok str = Tok.strtok(lw, Utils.OTHDELIMITERS);
        while (str != null) {
          seqAcc[i].append(str.toString());
          str = Tok.strtok(lw, Utils.OTHDELIMITERS);
        }

        if (seqAcc[i].charAt(seqAcc[i].length() - 1) == '*')
          seqAcc[i].setLength(seqAcc[i].length() - 1);
      }
    }

    return fillMatrices(true);
  }

  // ------------------------------------------------------------------ output

  public boolean printAlignment() throws IOException {
    if (sequences == null)
      return false;

    if ((residNumber == 0) || (sequenNumber == 0)) {
      try (CppOut err = cerr()) {
        err.endl().p("WARNING: Output alignment has not been generated. ")
            .p("It is empty").endl().endl();
      }
      return true;
    }

    try (CppOut cout = CppOut.toStdout()) {
      switch (oformat) {
        case 1:
          alignmentClustalToFile(cout);
          break;
        case 3:
          alignmentNBRF_PirToFile(cout);
          break;
        case 8:
          alignmentFastaToFile(cout);
          break;
        case 11:
          alignmentPhylip3_2ToFile(cout);
          break;
        case 12:
          alignmentPhylipToFile(cout);
          break;
        case 13:
          alignmentPhylip_PamlToFile(cout);
          break;
        case 17:
          alignmentNexusToFile(cout);
          break;
        case 21:
        case 22:
          alignmentMegaToFile(cout);
          break;
        case 99:
          getSequences(cout);
          break;
        case 100:
          alignmentColourHTML(cout);
          break;
        default:
          return false;
      }
    }
    return true;
  }

  public boolean saveAlignment(String destFile) throws IOException {
    if (sequences == null)
      return false;

    if ((residNumber == 0) || (sequenNumber == 0)) {
      try (CppOut err = cerr()) {
        err.endl().p("WARNING: Output alignment has not been generated. ")
            .p("It is empty.").endl().endl();
      }
      return true;
    }

    switch (oformat) {
      case 1:
      case 11:
      case 12:
      case 13:
      case 17:
      case 21:
      case 22:
        if (!isAligned) {
          // C++ omits endl after this message (alignment.cpp:632-634)
          try (CppOut err = cerr()) {
            err.endl().p("ERROR: Sequences are not aligned. Output format is ")
                .p("not compatible with unaligned sequences.");
          }
          return false;
        }
        break;
      default:
        break;
    }

    java.io.OutputStream destination;
    try {
      destination = new java.io.FileOutputStream(destFile);
    } catch (IOException e) {
      // ofstream::open reports failure to the caller without throwing.
      return false;
    }
    try (CppOut file = new CppOut(new java.io.BufferedOutputStream(destination, 1 << 16))) {
      switch (oformat) {
        case 1:
          alignmentClustalToFile(file);
          break;
        case 3:
          alignmentNBRF_PirToFile(file);
          break;
        case 8:
          alignmentFastaToFile(file);
          break;
        case 11:
          alignmentPhylip3_2ToFile(file);
          break;
        case 12:
          alignmentPhylipToFile(file);
          break;
        case 13:
          alignmentPhylip_PamlToFile(file);
          break;
        case 17:
          alignmentNexusToFile(file);
          break;
        case 21:
        case 22:
          alignmentMegaToFile(file);
          break;
        case 99:
          getSequences(file);
          break;
        case 100:
          alignmentColourHTML(file);
          break;
        default:
          return false;
      }
    }
    return true;
  }

  /** -onlyseqs writer is getSequences(ostream&) — see below in writers section. */

  // ------------------------------------------------------------------ writers

  /** -onlyseqs writer (rwAlignment.cpp:115-142): gaps removed, blank line per record. */
  public void getSequences(CppOut file) throws IOException {
    int i, j;
    String[] tmpMatrix = new String[sequenNumber];

    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? Utils.removeCharacter('-', sequences[i])
          : Utils.removeCharacter('-', Utils.getReverse(sequences[i]));

    for (i = 0; i < sequenNumber; i++) {
      file.p(">").p(seqsName[i]).endl();
      for (j = 0; j < tmpMatrix[i].length(); j += 60)
        file.p(sub(tmpMatrix[i], j, 60)).endl();
      file.endl();
    }
  }

  public void alignmentPhylipToFile(CppOut file) throws IOException {
    int i, j, maxLongName;

    if (!isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences are not aligned. Format (PHYLIP) ")
            .p("not compatible with unaligned sequences.").endl().endl();
      }
      return;
    }

    String[] tmpMatrix = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? sequences[i] : Utils.getReverse(sequences[i]);

    maxLongName = Utils.PHYLIPDISTANCE;
    for (i = 0; (i < sequenNumber) && (!shortNames); i++)
      maxLongName = Utils.max(maxLongName, seqsName[i].length());

    file.p(" ").p(sequenNumber).p(" ").p(residNumber).endl();

    for (i = 0; i < sequenNumber; i++)
      file.setw(maxLongName + 3).left().p(sub(seqsName[i], 0, maxLongName))
          .p(sub(tmpMatrix[i], 0, 60)).endl();
    file.endl();

    for (i = 60; i < residNumber; i += 60) {
      for (j = 0; j < sequenNumber; j++)
        file.p(sub(tmpMatrix[j], i, 60)).endl();
      file.endl();
    }
    file.endl();
  }

  public void alignmentPhylip3_2ToFile(CppOut file) throws IOException {
    int i, j, k, maxLongName;

    if (!isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences are not aligned. Format (PHYLIP) ")
            .p("not compatible with unaligned sequences.").endl().endl();
      }
      return;
    }

    String[] tmpMatrix = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? sequences[i] : Utils.getReverse(sequences[i]);

    maxLongName = Utils.PHYLIPDISTANCE;
    for (i = 0; (i < sequenNumber) && (!shortNames); i++)
      maxLongName = Utils.max(maxLongName, seqsName[i].length());

    file.p(" ").p(sequenNumber).p(" ").p(residNumber).endl();

    for (i = 0; i < sequenNumber; i++) {
      file.setw(maxLongName + 3).left().p(sub(seqsName[i], 0, maxLongName));
      for (j = 0; j < residNumber; j += 50) {
        // C++ bug kept: body reads sequences[i], not the (possibly reversed) tmpMatrix
        for (k = j; (k < residNumber) && (k < (j + 50)); k += 10)
          file.p(sub(sequences[i], k, 10)).p(" ");
        file.endl();
        if ((j + 50) < residNumber)
          file.setw(maxLongName + 3).p(" ");
      }
      file.endl();
    }
  }

  public void alignmentPhylip_PamlToFile(CppOut file) throws IOException {
    int i, maxLongName;

    if (!isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences are not aligned. Format (PHYLIP) ")
            .p("not compatible with unaligned sequences.").endl().endl();
      }
      return;
    }

    String[] tmpMatrix = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? sequences[i] : Utils.getReverse(sequences[i]);

    maxLongName = Utils.PHYLIPDISTANCE;
    for (i = 0; (i < sequenNumber) && (!shortNames); i++)
      maxLongName = Utils.max(maxLongName, seqsName[i].length());

    file.p(" ").p(sequenNumber).p(" ").p(residNumber).endl();

    for (i = 0; i < sequenNumber; i++)
      // C++ bug kept: sequences[i], not tmpMatrix[i]
      file.setw(maxLongName + 3).left().p(sub(seqsName[i], 0, maxLongName))
          .p(sequences[i]).endl();
    file.endl();
  }

  public void alignmentClustalToFile(CppOut file) throws IOException {
    int i, j, maxLongName = 0;

    if (!isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences are not aligned. Format (CLUSTAL) ")
            .p("not compatible with unaligned sequences.").endl().endl();
      }
      return;
    }

    String[] tmpMatrix = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? sequences[i] : Utils.getReverse(sequences[i]);

    for (i = 0; (i < sequenNumber) && (!shortNames); i++)
      maxLongName = Utils.max(maxLongName, seqsName[i].length());

    if ((aligInfo.length() != 0) && (iformat == oformat))
      file.p(aligInfo).endl().endl();
    else
      file.p("CLUSTAL multiple sequence alignment").endl().endl();

    for (j = 0; j < residNumber; j += 60) {
      for (i = 0; i < sequenNumber; i++)
        file.setw(maxLongName + 5).left().p(seqsName[i])
            .p(sub(tmpMatrix[i], j, 60)).endl();
      file.endl().endl();
    }
  }

  public void alignmentFastaToFile(CppOut file) throws IOException {
    int i, j, maxLongName;

    String[] tmpMatrix = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? sequences[i] : Utils.getReverse(sequences[i]);

    maxLongName = 0;
    for (i = 0; i < sequenNumber; i++)
      if (!keepHeader)
        maxLongName = Utils.max(maxLongName, seqsName[i].length());
      else if (seqsInfo != null)
        maxLongName = Utils.max(maxLongName, seqsInfo[i].length());

    if (shortNames && maxLongName > Utils.PHYLIPDISTANCE) {
      maxLongName = Utils.PHYLIPDISTANCE;
      if (keepHeader) {
        try (CppOut err = cerr()) {
          err.endl().p("WARNING: Original sequence header will be cut by charac")
              .p("ter 10").endl();
        }
      }
    }

    for (i = 0; i < sequenNumber; i++) {
      if (!keepHeader)
        file.p(">").p(sub(seqsName[i], 0, maxLongName)).endl();
      else if (seqsInfo != null)
        file.p(">").p(sub(seqsInfo[i], 0, maxLongName)).endl();
      for (j = 0; j < residuesNumber[i]; j += 60)
        file.p(sub(tmpMatrix[i], j, 60)).endl();
    }
  }

  public void alignmentNexusToFile(CppOut file) throws IOException {
    int i, j, k, maxLongName = 0;

    if (!isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences are not aligned. Format (NEXUS) ")
            .p("not compatible with unaligned sequences.").endl().endl();
      }
      return;
    }

    String[] tmpMatrix = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? sequences[i] : Utils.getReverse(sequences[i]);

    for (i = 0; (i < sequenNumber) && (!shortNames); i++)
      maxLongName = Utils.max(maxLongName, seqsName[i].length());

    getTypeAlignment();

    while (aligInfo.indexOf(";") != -1)
      aligInfo = aligInfo.substring(0, aligInfo.indexOf(";"))
          + aligInfo.substring(aligInfo.indexOf(";") + 1);

    file.p("#NEXUS").endl().p("BEGIN DATA;").endl()
        .p(" DIMENSIONS NTAX=").p(sequenNumber).p(" NCHAR=").p(residNumber).p(";").endl();

    if ((dataType == Utils.DNAType) || (dataType == Utils.DNADeg))
      file.p("FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-");
    else if ((dataType == Utils.RNAType) || (dataType == Utils.RNADeg))
      file.p("FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-");
    else if (dataType == Utils.AAType)
      file.p("FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-");

    i = 0;
    while ((j = aligInfo.indexOf(" ", i)) != -1) {
      String w1 = sub(aligInfo, i, j - i);
      String w2 = sub(aligInfo, i, j);
      if (w1.startsWith("MISSING") || w2.startsWith("missing"))
        file.p(" ").p(w1);
      else if (w2.startsWith("MATCHCHAR") || w2.startsWith("matchchar"))
        file.p(" ").p(w1);
      i = j + 1;
    }
    file.p(";").endl();

    for (i = 0; i < sequenNumber; i++)
      file.p("[Name: ").setw(maxLongName + 4).left().p(seqsName[i]).p("Len: ")
          .p(residNumber).p("]").endl();
    file.endl().p("MATRIX").endl();

    for (j = 0; j < residNumber; j += 50) {
      for (i = 0; i < sequenNumber; i++) {
        file.setw(maxLongName + 4).left().p(seqsName[i]);
        for (k = j; k < (j + 50) && k < residNumber; k += 10)
          file.p(" ").p(sub(sequences[i], k, 10));
        file.endl();
      }
      file.endl();
    }
    file.p(";").endl().p("END;").endl();
  }

  public void alignmentMegaToFile(CppOut file) throws IOException {
    int i, j, k;

    if (!isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences are not aligned. Format (MEGA) ")
            .p("not compatible with unaligned sequences.").endl().endl();
      }
      return;
    }

    String[] tmpMatrix = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? sequences[i] : Utils.getReverse(sequences[i]);

    getTypeAlignment();

    file.p("#MEGA").endl().p(filename).endl();

    if ((dataType == Utils.DNAType) || (dataType == Utils.DNADeg))
      file.p("!Format DataType=DNA ");
    else if ((dataType == Utils.RNAType) || (dataType == Utils.RNADeg))
      file.p("!Format DataType=RNA ");
    else if (dataType == Utils.AAType)
      file.p("!Format DataType=protein ");

    file.p("NSeqs=").p(sequenNumber).p(" Nsites=").p(residNumber)
        .p(" indel=- CodeTable=Standard;").endl().endl();

    for (i = 0; i < sequenNumber; i++) {
      file.p("#").p(seqsName[i]).endl();
      for (j = 0; j < residNumber; j += 50) {
        for (k = j; ((k < residNumber) && (k < j + 50)); k += 10)
          file.p(sub(tmpMatrix[i], k, 10)).p(" ");
        file.endl();
      }
      file.endl();
    }
  }

  public void alignmentNBRF_PirToFile(CppOut file) throws IOException {
    int i, j, k;
    String alg_datatype = "";

    String[] tmpMatrix = new String[sequenNumber];
    for (i = 0; i < sequenNumber; i++)
      tmpMatrix[i] = (!reverse) ? sequences[i] : Utils.getReverse(sequences[i]);

    getTypeAlignment();
    if ((dataType == Utils.DNAType) || (dataType == Utils.DNADeg))
      alg_datatype = "DL";
    else if ((dataType == Utils.RNAType) || (dataType == Utils.RNADeg))
      alg_datatype = "RL";
    else if (dataType == Utils.AAType)
      alg_datatype = "P1";

    for (i = 0; i < sequenNumber; i++) {
      if ((seqsInfo != null) && (iformat == oformat))
        file.p(">").p(sub(seqsInfo[i], 0, 2)).p(";").p(seqsName[i])
            .endl().p(sub(seqsInfo[i], 2, Integer.MAX_VALUE)).endl();
      else
        file.p(">").p(alg_datatype).p(";").p(seqsName[i]).endl()
            .p(seqsName[i]).p(" ").p(residuesNumber[i]).p(" bases").endl();

      for (j = 0; j < residuesNumber[i]; j += 50) {
        for (k = j; (k < residuesNumber[i]) && (k < (j + 50)); k += 10)
          file.p(" ").p(sub(tmpMatrix[i], k, 10));

        if (k >= residuesNumber[i]) {
          if ((residuesNumber[i] % 50) == 0)
            file.endl().p(" ");
          else if ((residuesNumber[i] % 10) == 0)
            file.p(" ");
          file.p("*");
        }
        file.endl();
      }
      file.endl();
    }
  }

  public boolean alignmentColourHTML(CppOut file) throws IOException {
    int i, j, kj, upper, k = 0, maxLongName = 0;
    char type;

    if (!isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences are not aligned.").endl().endl();
      }
      return false;
    }

    for (i = 0; i < sequenNumber; i++)
      maxLongName = Utils.max(maxLongName, seqsName[i].length());

    file.p("<!DOCTYPE html>").endl().p("<html><head>").endl().p("    <meta ")
        .p("http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />")
        .endl().p("    <title>readAl v1.5.1</title>").endl()
        .p("    <style type=\"text/css\">").endl()
        .p("    #b  { background-color: #3366ff; }\n")
        .p("    #r  { background-color: #cc0000; }\n")
        .p("    #g  { background-color: #33cc00; }\n")
        .p("    #p  { background-color: #ff6666; }\n")
        .p("    #m  { background-color: #cc33cc; }\n")
        .p("    #o  { background-color: #ff9900; }\n")
        .p("    #c  { background-color: #46C7C7; }\n")
        .p("    #y  { background-color: #FFFF00; }\n")
        .p("    </style>\n  </head>\n\n").p("  <body>\n  <pre>").endl();

    for (j = 0, upper = Utils.HTMLBLOCKS; j < residNumber; j += Utils.HTMLBLOCKS, upper += Utils.HTMLBLOCKS) {
      file.endl();
      file.setw(maxLongName + 19).right().p(j + 10);
      for (i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
        file.setw(10).right().p(i);

      file.endl().setw(maxLongName + 10);
      for (i = j + 1; ((i <= residNumber) && (i <= upper)); i++)
        file.p((i % 10) != 0 ? "=" : "+");

      // precompute each column once per block: the writer used to rebuild the whole
      // column string for every residue, making HTML output O(sequences^2 x columns)
      String[] blockColumns = new String[Math.max(0, Math.min(residNumber, upper) - j)];
      for (k = j; ((k < residNumber) && (k < upper)); k++) {
        StringBuilder column = new StringBuilder(sequenNumber);
        for (kj = 0; kj < sequenNumber; kj++)
          column.append(Utils.charAt(sequences[kj], k));
        blockColumns[k - j] = column.toString();
      }

      for (i = 0; i < sequenNumber; i++) {
        file.endl().setw(maxLongName + 9).left().p(seqsName[i]);

        for (k = j; ((k < residNumber) && (k < upper)); k++) {
          type = Utils.determineColor(Utils.charAt(sequences[i], k), blockColumns[k - j]);
          if (type == 'w')
            file.p(Utils.charAt(sequences[i], k));
          else
            file.p("<span id=").p(type).p(">").p(Utils.charAt(sequences[i], k)).p("</span>");
        }
      }
      file.endl();
    }

    file.p("    </pre>").endl().p("  </body>").endl().p("</html>").endl();
    return true;
  }

  public boolean alignmentSummaryHTML(String destFile, int residues, int seqs,
      int[] selectedRes, int[] selectedSeq, float[] consValues) throws IOException {

    int i, j, k, kj, upper, minHTML, maxLongName;
    float[] simValues;
    boolean[] res, seq;
    int[] gapsValues;
    char type;

    if (!isAligned) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Sequences are not aligned.").endl().endl();
      }
      return false;
    }

    java.io.OutputStream destination;
    try {
      destination = new java.io.FileOutputStream(destFile);
    } catch (IOException e) {
      // ofstream::open reports failure to the caller without throwing.
      return false;
    }
    try (CppOut file = new CppOut(new java.io.BufferedOutputStream(destination, 1 << 16))) {

      maxLongName = 0;
      for (i = 0; i < sequenNumber; i++)
        maxLongName = Utils.max(maxLongName, seqsName[i].length());

      minHTML = Utils.max(25, maxLongName + 10);

      res = new boolean[residNumber];
      seq = new boolean[sequenNumber];

      for (i = 0; i < residues; i++)
        res[selectedRes[i]] = true;
      for (i = 0; i < seqs; i++)
        seq[selectedSeq[i]] = true;

      gapsValues = null;
      if (sgaps != null)
        gapsValues = sgaps.getGapsWindow();
      simValues = null;
      if (scons != null)
        simValues = scons.getMdkwVector();

      file.p("<!DOCTYPE html>").endl().p("<html><head>").endl().p("    <meta ")
          .p("http-equiv=\"Content-Type\" content=\"text/html;charset=ISO-8859-1\" />")
          .endl().p("    <title>trimAl v1.5.1 Summary</title>").endl()
          .p("    <style type=\"text/css\" media=\"all\">").endl()
          .p("    #b  { background-color: #3366ff; }\n")
          .p("    #r  { background-color: #cc0000; }\n")
          .p("    #g  { background-color: #33cc00; }\n")
          .p("    #p  { background-color: #ff6666; }\n")
          .p("    #m  { background-color: #cc33cc; }\n")
          .p("    #o  { background-color: #ff9900; }\n")
          .p("    #c  { background-color: #46C7C7; }\n")
          .p("    #y  { background-color: #FFFF00; }\n")
          .p("    .sel  { background-color: #B9B9B9; }\n")
          .p("    .nsel { background-color: #E9E9E9; }\n")
          .p("    .c1   { background-color: #FFFBF2; }\n")
          .p("    .c2   { background-color: #FFF8CC; }\n")
          .p("    .c3   { background-color: #FAF0BE; }\n")
          .p("    .c4   { background-color: #F0EAD6; }\n")
          .p("    .c5   { background-color: #F3E5AB; }\n")
          .p("    .c6   { background-color: #F4C430; }\n")
          .p("    .c7   { background-color: #C2B280; color: white; }\n")
          .p("    .c8   { background-color: #DAA520; color: white; }\n")
          .p("    .c9   { background-color: #B8860B; color: white; }\n")
          .p("    .c10  { background-color: #918151; color: white; }\n")
          .p("    .c11  { background-color: #967117; color: white; }\n")
          .p("    .c12  { background-color: #6E5411; color: white; }\n")
          .p("    </style>\n  </head>\n\n").p("  <body>\n").p("  <pre>").endl();

      file.p("    <span class=sel>Selected Sequences: ").setw(5).right().p(seqs)
          .p(" /Selected Residues: ").setw(7).right().p(residues).p("</span>")
          .endl().p("    <span class=nsel>Deleted Sequences:  ").setw(5).right()
          .p(sequenNumber - seqs).p(" /Deleted Residues:  ").setw(7).right()
          .p(residNumber - residues).p("</span>").endl();

      if (gapsValues != null)
        file.endl().setw(minHTML).left().p("    Gaps Scores:        ")
            .p("<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>")
            .p("<span  class=c3> <.050 </span><span  class=c4> <.100 </span>")
            .p("<span  class=c5> <.150 </span><span  class=c6> <.200 </span>")
            .p("<span  class=c7> <.250 </span><span  class=c8> <.350 </span>")
            .p("<span  class=c9> <.500 </span><span class=c10> <.750 </span>")
            .p("<span class=c11> <1.00 </span><span class=c12>  =1=  </span>");

      if (simValues != null)
        file.endl().setw(minHTML).left().p("    Similarity Scores:  ")
            .p("<span  class=c1>  =0=  </span><span  class=c2> <1e-6 </span>")
            .p("<span  class=c3> <1e-5 </span><span  class=c4> <1e-4 </span>")
            .p("<span  class=c5> <.001 </span><span  class=c6> <.010 </span>")
            .p("<span  class=c7> <.100 </span><span  class=c8> <.250 </span>")
            .p("<span  class=c9> <.500 </span><span class=c10> <.750 </span>")
            .p("<span class=c11> <1.00 </span><span class=c12>  =1=  </span>");

      if (consValues != null)
        file.endl().setw(minHTML).left().p("    Consistency Scores: ")
            .p("<span  class=c1>  =0=  </span><span  class=c2> <.001 </span>")
            .p("<span  class=c3> <.050 </span><span  class=c4> <.100 </span>")
            .p("<span  class=c5> <.150 </span><span  class=c6> <.200 </span>")
            .p("<span  class=c7> <.250 </span><span  class=c8> <.350 </span>")
            .p("<span  class=c9> <.500 </span><span class=c10> <.750 </span>")
            .p("<span class=c11> <1.00 </span><span class=c12>  =1=  </span>");

      if ((gapsValues != null) || (simValues == null) || (consValues == null))
        file.endl();

      for (j = 0, upper = Utils.HTMLBLOCKS; j < residNumber; j += Utils.HTMLBLOCKS, upper += Utils.HTMLBLOCKS) {

        file.endl().setw(minHTML + 10).right().p(j + 10);
        for (i = j + 20; ((i <= residNumber) && (i <= upper)); i += 10)
          file.setw(10).right().p(i);

        file.endl().setw(minHTML + 1).right();
        for (i = j + 1; ((i <= residNumber) && (i <= upper)); i++)
          file.p((i % 10) != 0 ? "=" : "+");
        file.endl();

        // precompute each column once per block (was rebuilt per residue, O(seq^2 x col))
        String[] blockColumns = new String[Math.max(0, Math.min(residNumber, upper) - j)];
        for (k = j; ((k < residNumber) && (k < upper)); k++) {
          StringBuilder column = new StringBuilder(sequenNumber);
          for (kj = 0; kj < sequenNumber; kj++)
            column.append(Utils.charAt(sequences[kj], k));
          blockColumns[k - j] = column.toString();
        }

        for (i = 0; i < sequenNumber; i++) {
          file.p("    <span class=").p(seq[i] ? "sel>" : "nsel>").p(seqsName[i])
              .p("</span>").setw(minHTML - 4 - seqsName[i].length()).right().p("");

          for (k = j; ((k < residNumber) && (k < upper)); k++) {
            type = Utils.determineColor(Utils.charAt(sequences[i], k), blockColumns[k - j]);
            if (type == 'w')
              file.p(Utils.charAt(sequences[i], k));
            else
              file.p("<span id=").p(type).p(">").p(Utils.charAt(sequences[i], k)).p("</span>");
          }
          file.endl();
        }

        file.endl().setw(minHTML).left().p("    Selected Cols:      ");
        for (k = j; ((k < residNumber) && (k < (j + Utils.HTMLBLOCKS))); k++)
          file.p("<span class=").p(res[k] ? "sel" : "nsel").p("> </span>");
        file.endl();

        if ((gapsValues == null) && (simValues == null) && (consValues == null))
          continue;

        if (gapsValues != null) {
          file.endl().setw(minHTML).left().p("    Gaps Scores:        ");
          for (k = j; ((k < residNumber) && (k < (j + Utils.HTMLBLOCKS))); k++)
            if (gapsValues[k] == 0)
              file.p("<span class=c12> </span>");
            else if (gapsValues[k] == sequenNumber)
              file.p("<span class=c1> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .750)
              file.p("<span class=c11> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .500)
              file.p("<span class=c10> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .350)
              file.p("<span  class=c9> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .250)
              file.p("<span  class=c8> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .200)
              file.p("<span  class=c7> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .150)
              file.p("<span  class=c6> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .100)
              file.p("<span  class=c5> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .050)
              file.p("<span  class=c4> </span>");
            else if (1 - (float) gapsValues[k] / sequenNumber >= .001)
              file.p("<span  class=c3> </span>");
            else
              file.p("<span  class=c2> </span>");
        }
        if (simValues != null) {
          file.endl().setw(minHTML).left().p("    Similarity Scores:  ");
          for (k = j; ((k < residNumber) && (k < (j + Utils.HTMLBLOCKS))); k++)
            if (simValues[k] == 1)
              file.p("<span class=c12> </span>");
            else if (simValues[k] == 0)
              file.p("<span class=c1> </span>");
            else if (simValues[k] >= .750)
              file.p("<span class=c11> </span>");
            else if (simValues[k] >= .500)
              file.p("<span class=c10> </span>");
            else if (simValues[k] >= .250)
              file.p("<span  class=c9> </span>");
            else if (simValues[k] >= .100)
              file.p("<span  class=c8> </span>");
            else if (simValues[k] >= .010)
              file.p("<span  class=c7> </span>");
            else if (simValues[k] >= .001)
              file.p("<span  class=c6> </span>");
            else if (simValues[k] >= 1e-4)
              file.p("<span  class=c5> </span>");
            else if (simValues[k] >= 1e-5)
              file.p("<span  class=c4> </span>");
            else if (simValues[k] >= 1e-6)
              file.p("<span  class=c3> </span>");
            else
              file.p("<span  class=c2> </span>");
        }
        if (consValues != null) {
          file.endl().setw(minHTML).left().p("    Consistency Scores: ");
          for (k = j; ((k < residNumber) && (k < (j + Utils.HTMLBLOCKS))); k++)
            if (consValues[k] == 1)
              file.p("<span class=c12> </span>");
            else if (consValues[k] == 0)
              file.p("<span class=c1> </span>");
            else if (consValues[k] >= .750)
              file.p("<span class=c11> </span>");
            else if (consValues[k] >= .500)
              file.p("<span class=c10> </span>");
            else if (consValues[k] >= .350)
              file.p("<span  class=c9> </span>");
            else if (consValues[k] >= .250)
              file.p("<span  class=c8> </span>");
            else if (consValues[k] >= .200)
              file.p("<span  class=c7> </span>");
            else if (consValues[k] >= .150)
              file.p("<span  class=c6> </span>");
            else if (consValues[k] >= .100)
              file.p("<span  class=c5> </span>");
            else if (consValues[k] >= .050)
              file.p("<span  class=c4> </span>");
            else if (consValues[k] >= .001)
              file.p("<span  class=c3> </span>");
            else
              file.p("<span  class=c2> </span>");
        }
        file.endl();
      }

      file.p("    </pre>").endl().p("  </body>").endl().p("</html>").endl();
    }

    return true;
  }

  // ------------------------------------------------------------------ trimming

  public Alignment cleanGaps(float baseLine, float gapsPct, boolean complementary) throws IOException {
    Alignment ret;
    double cut;

    if (calculateGapStats() != true)
      return null;

    cut = sgaps.calcCutPoint(baseLine, gapsPct);
    ret = cleanByCutValue(cut, baseLine, sgaps.getGapsWindow(), complementary);
    return ret;
  }

  public Alignment cleanConservation(float baseLine, float conservationPct, boolean complementary) throws IOException {
    Alignment ret;
    float cut;

    if (calculateConservationStats() != true)
      return null;

    cut = (float) scons.calcCutPoint(baseLine, conservationPct);
    ret = cleanByCutValue(cut, baseLine, scons.getMdkwVector(), complementary);
    return ret;
  }

  public Alignment clean(float baseLine, float GapsPct, float conservationPct, boolean complementary) throws IOException {
    Alignment ret;
    float cutCons;
    double cutGaps;

    if (calculateGapStats() != true)
      return null;

    if (calculateConservationStats() != true)
      return null;

    cutGaps = sgaps.calcCutPoint(baseLine, GapsPct);
    cutCons = (float) scons.calcCutPoint(baseLine, conservationPct);
    ret = cleanByCutValue(cutGaps, sgaps.getGapsWindow(), baseLine, cutCons, scons.getMdkwVector(), complementary);
    return ret;
  }

  public Alignment cleanCompareFile(float cutpoint, float baseLine, float[] vectValues, boolean complementary) throws IOException {
    Alignment ret;
    float cut;

    float[] vectAux = new float[residNumber];
    Utils.copyVect(vectValues, vectAux, residNumber);
    Utils.quicksort(vectAux, 0, residNumber - 1);
    cut = vectAux[(int) ((float) (residNumber - 1) * (100.0 - baseLine) / 100.0)];
    cut = cutpoint < cut ? cutpoint : cut;
    ret = cleanByCutValue(cut, baseLine, vectValues, complementary);
    return ret;
  }

  public Alignment cleanSpuriousSeq(float overlapColumn, float minimumOverlap, boolean complementary) throws IOException {
    float[] overlapVector = new float[sequenNumber];
    Alignment newAlig;

    if (!calculateSpuriousVector(overlapColumn, overlapVector))
      return null;

    newAlig = cleanOverlapSeq(minimumOverlap, overlapVector, complementary);
    return newAlig;
  }

  public Alignment clean2ndSlope(boolean complementarity) throws IOException {
    Alignment ret;
    int cut;

    if (calculateGapStats() != true)
      return null;

    cut = sgaps.calcCutPoint2ndSlope();
    ret = cleanByCutValue((double) cut, 0, sgaps.getGapsWindow(), complementarity);
    return ret;
  }

  public Alignment cleanCombMethods(boolean complementarity, boolean variable) throws IOException {
    float simCut, first20Point, last80Point;
    float[] simil, vectAux;
    int i, j, acm, gapCut;
    int[] positions, gaps;
    double inic, fin, vlr;

    if (calculateGapStats() != true)
      return null;

    gapCut = sgaps.calcCutPoint2ndSlope();
    gaps = sgaps.getGapsWindow();

    if (calculateConservationStats() != true)
      return null;

    simil = scons.getMdkwVector();

    positions = new int[residNumber];
    Utils.initlVect(positions, residNumber, -1);

    for (i = 0, acm = 0; i < residNumber; i++) {
      if (gaps[i] <= gapCut) {
        positions[i] = i;
        acm++;
      }
    }

    vectAux = new float[acm];
    for (i = 0, j = 0; i < residNumber; i++)
      if (positions[i] != -1)
        vectAux[j++] = simil[i];

    Utils.quicksort(vectAux, 0, acm - 1);

    first20Point = 0;
    last80Point = 0;

    for (i = acm - 1, j = 1; i >= 0; i--, j++) {
      if ((((float) j / acm) * 100.0) <= 20.0)
        first20Point = vectAux[i];
      if ((((float) j / acm) * 100.0) <= 80.0)
        last80Point = vectAux[i];
    }

    // <cmath> selects log10(float), returning float before assignment to double.
    inic = Utils.log10f(first20Point);
    fin = Utils.log10f(last80Point);
    vlr = ((inic - fin) / 10) + fin;
    simCut = (float) Math.pow(10, vlr);

    return cleanStrict(gapCut, sgaps.getGapsWindow(),
        simCut, scons.getMdkwVector(), complementarity, variable);
  }

  public Alignment cleanNoAllGaps(boolean complementarity) throws IOException {
    Alignment ret;

    if (calculateGapStats() != true)
      return null;

    ret = cleanByCutValue((double) (sequenNumber - 1), 0, sgaps.getGapsWindow(), complementarity);
    return ret;
  }

  /** cleanByCutValue, gap version (alignment.cpp:2073-2231). */
  public Alignment cleanByCutValue(double cut, float baseLine,
      int[] gInCol, boolean complementary) throws IOException {

    int i, j, k, jn, oth, pos, block;
    int[] vectAux;
    String[] matrixAux, newSeqsName, newSeqsInfo;
    Alignment newAlig;
    NewValues counter = new NewValues();

    for (i = 0, counter.residues = 0; i < residNumber; i++)
      if (gInCol[i] <= cut) counter.residues++;
      else saveResidues[i] = -1;

    oth = Utils.roundInt(((baseLine / 100.0) - (float) counter.residues / residNumber) * residNumber);

    if (oth > 0) {
      counter.residues += oth;

      vectAux = new int[residNumber];
      Utils.copyVect(gInCol, vectAux, residNumber);
      Utils.quicksort(vectAux, 0, residNumber - 1);
      cut = vectAux[(int) ((float) (residNumber - 1) * (baseLine) / 100.0)];
    }

    for (k = Utils.roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

      for (i = (residNumber / 2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

        for (jn = i; ((jn >= 0) && (saveResidues[jn] != -1) && (oth > 0)); jn--)
          ;

        if ((i - jn) >= k) {
          for (; ((jn >= 0) && (saveResidues[jn] == -1) && (oth > 0)); jn--) {
            if (gInCol[jn] <= cut) {
              saveResidues[jn] = jn;
              oth--;
            } else
              break;
          }
        }
        i = jn;

        for (jn = j; ((jn < residNumber) && (saveResidues[jn] != -1) && (oth > 0)); jn++)
          ;

        if ((jn - j) >= k) {
          for (; ((jn < residNumber) && (saveResidues[jn] == -1) && (oth > 0)); jn++) {
            if (gInCol[jn] <= cut) {
              saveResidues[jn] = jn;
              oth--;
            } else
              break;
          }
        }
        j = jn;
      }
    }

    if (blockSize != 0) {
      for (i = 0, pos = 0, block = 0; i < residNumber; i++) {
        if (saveResidues[i] != -1)
          block++;
        else {
          if (block < blockSize)
            for (j = pos; j <= i; j++)
              saveResidues[j] = -1;
          pos = i + 1;
          block = 0;
        }
      }
      if (block < blockSize)
        for (j = pos; j < i; j++)
          saveResidues[j] = -1;
    }

    if (terminalGapOnly == true)
      if (!removeOnlyTerminal())
        return null;

    if (complementary == true)
      computeComplementaryAlig(true, false);

    counter = removeCols_SeqsAllGaps();

    matrixAux = new String[counter.sequences];
    newSeqsName = new String[counter.sequences];
    newSeqsInfo = new String[counter.sequences];
    for (int z = 0; z < counter.sequences; z++) {
      matrixAux[z] = "";
      newSeqsName[z] = "";
      newSeqsInfo[z] = "";
    }

    fillNewDataStructure(matrixAux, newSeqsName, newSeqsInfo);

    newAlig = new Alignment(filename, aligInfo, matrixAux, newSeqsName, newSeqsInfo,
        counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
        isAligned, reverse, terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber, residNumber,
        residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
        identities, overlaps);

    return newAlig;
  }

  /** cleanByCutValue, conservation version (alignment.cpp:2240-2385). */
  public Alignment cleanByCutValue(float cut, float baseLine,
      float[] ValueVect, boolean complementary) throws IOException {

    int i, j, k, jn, oth, pos, block;
    String[] matrixAux, newSeqsName, newSeqsInfo;
    Alignment newAlig;
    NewValues counter = new NewValues();

    for (i = 0, counter.residues = 0; i < residNumber; i++)
      if (ValueVect[i] > cut) counter.residues++;
      else saveResidues[i] = -1;

    oth = Utils.roundInt(((baseLine / 100.0) - (float) counter.residues / residNumber) * residNumber);
    if (oth > 0) counter.residues += oth;

    for (k = Utils.roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

      for (i = (residNumber / 2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

        for (jn = i; ((jn >= 0) && (saveResidues[jn] != -1) && (oth > 0)); jn--)
          ;

        if ((i - jn) >= k) {
          for (; ((jn >= 0) && (saveResidues[jn] == -1) && (oth > 0)); jn--) {
            if (ValueVect[jn] == cut) {
              saveResidues[jn] = jn;
              oth--;
            } else
              break;
          }
        }
        i = jn;

        for (jn = j; ((jn < residNumber) && (saveResidues[jn] != -1) && (oth > 0)); jn++)
          ;

        if ((jn - j) >= k) {
          for (; ((jn < residNumber) && (saveResidues[jn] == -1) && (oth > 0)); jn++) {
            if (ValueVect[jn] == cut) {
              saveResidues[jn] = jn;
              oth--;
            } else
              break;
          }
        }
        j = jn;
      }
    }

    if (blockSize != 0) {
      for (i = 0, pos = 0, block = 0; i < residNumber; i++) {
        if (saveResidues[i] != -1)
          block++;
        else {
          if (block < blockSize)
            for (j = pos; j <= i; j++)
              saveResidues[j] = -1;
          pos = i + 1;
          block = 0;
        }
      }
      if (block < blockSize)
        for (j = pos; j < i; j++)
          saveResidues[j] = -1;
    }

    if (terminalGapOnly == true)
      if (!removeOnlyTerminal())
        return null;

    if (complementary == true)
      computeComplementaryAlig(true, false);

    counter = removeCols_SeqsAllGaps();

    matrixAux = new String[counter.sequences];
    newSeqsName = new String[counter.sequences];
    newSeqsInfo = new String[counter.sequences];
    for (int z = 0; z < counter.sequences; z++) {
      matrixAux[z] = "";
      newSeqsName[z] = "";
      newSeqsInfo[z] = "";
    }

    fillNewDataStructure(matrixAux, newSeqsName, newSeqsInfo);

    newAlig = new Alignment(filename, aligInfo, matrixAux, newSeqsName, newSeqsInfo,
        counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
        isAligned, reverse, terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber, residNumber,
        residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
        identities, overlaps);

    return newAlig;
  }

  /** cleanByCutValue, combined gaps+conservation version (alignment.cpp:2398-2580). */
  public Alignment cleanByCutValue(double cutGaps, int[] gInCol,
      float baseLine, float cutCons, float[] MDK_Win, boolean complementary) throws IOException {

    int i, j, k, oth, pos, block, jn, blGaps = 0;
    int[] vectAuxGaps;
    float blCons = 0;
    float[] vectAuxCons;
    String[] matrixAux, newSeqsName, newSeqsInfo;
    Alignment newAlig;
    NewValues counter = new NewValues();

    for (i = 0, counter.residues = 0; i < residNumber; i++)
      if ((MDK_Win[i] > cutCons) && (gInCol[i] <= cutGaps)) counter.residues++;
      else saveResidues[i] = -1;

    oth = Utils.roundInt(((baseLine / 100.0) - (float) counter.residues / residNumber) * residNumber);

    if (oth > 0) {
      counter.residues += oth;

      vectAuxCons = new float[residNumber];
      vectAuxGaps = new int[residNumber];

      Utils.copyVect(MDK_Win, vectAuxCons, residNumber);
      Utils.copyVect(gInCol, vectAuxGaps, residNumber);

      Utils.quicksort(vectAuxCons, 0, residNumber - 1);
      Utils.quicksort(vectAuxGaps, 0, residNumber - 1);

      blCons = vectAuxCons[(int) ((float) (residNumber - 1) * (100.0 - baseLine) / 100.0)];
      blGaps = vectAuxGaps[(int) ((float) (residNumber - 1) * (baseLine) / 100.0)];
    }

    for (k = Utils.roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--) {

      for (i = (residNumber / 2), j = (i + 1); (((i > 0) || (j < (residNumber - 1))) && (oth > 0)); i--, j++) {

        for (jn = i; ((jn >= 0) && (saveResidues[jn] != -1) && (oth > 0)); jn--)
          ;

        if ((i - jn) >= k) {
          for (; ((jn >= 0) && (saveResidues[jn] == -1) && (oth > 0)); jn--) {
            if ((MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)) {
              saveResidues[jn] = jn;
              oth--;
            } else
              break;
          }
        }
        i = jn;

        for (jn = j; ((jn < residNumber) && (saveResidues[jn] != -1) && (oth > 0)); jn++)
          ;

        if ((jn - j) >= k) {
          for (; ((jn < residNumber) && (saveResidues[jn] == -1) && (oth > 0)); jn++) {
            if ((MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)) {
              saveResidues[jn] = jn;
              oth--;
            } else
              break;
          }
        }
        j = jn;
      }
    }

    if (blockSize != 0) {
      for (i = 0, pos = 0, block = 0; i < residNumber; i++) {
        if (saveResidues[i] != -1)
          block++;
        else {
          if (block < blockSize)
            for (j = pos; j <= i; j++)
              saveResidues[j] = -1;
          pos = i + 1;
          block = 0;
        }
      }
      if (block < blockSize)
        for (j = pos; j < i; j++)
          saveResidues[j] = -1;
    }

    if (terminalGapOnly == true)
      if (!removeOnlyTerminal())
        return null;

    if (complementary == true)
      computeComplementaryAlig(true, false);

    counter = removeCols_SeqsAllGaps();

    matrixAux = new String[counter.sequences];
    newSeqsName = new String[counter.sequences];
    newSeqsInfo = new String[counter.sequences];
    for (int z = 0; z < counter.sequences; z++) {
      matrixAux[z] = "";
      newSeqsName[z] = "";
      newSeqsInfo[z] = "";
    }

    fillNewDataStructure(matrixAux, newSeqsName, newSeqsInfo);

    newAlig = new Alignment(filename, aligInfo, matrixAux, newSeqsName, newSeqsInfo,
        counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
        isAligned, reverse, terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber, residNumber,
        residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
        identities, overlaps);

    return newAlig;
  }

  public Alignment cleanStrict(int gapCut, int[] gInCol, float simCut,
      float[] MDK_W, boolean complementary, boolean variable) throws IOException {

    int i, num, lenBlock;
    Alignment newAlig;

    IntDeque neighboursBlock = new IntDeque(5);

    // Original unconditionally reads indices 0/1 and length-2: UB below length 2.
    // Guard only those reads; all defined-size neighbor recovery is unchanged.
    for (i = 0; i < residNumber; i++)
      if (gInCol[i] > gapCut)
        saveResidues[i] = -1;

    for (i = 0; i < residNumber; i++)
      if (MDK_W[i] < simCut)
        saveResidues[i] = -1;

    for (i = 0; (i < residNumber) && (i < 5); i++)
      neighboursBlock.pushBack(saveResidues[i]);

    if ((residNumber > 0) && (saveResidues[0] == -1) && (neighboursBlock.size() > 2) &&
        (neighboursBlock.get(1) != -1) &&
        (neighboursBlock.get(2) != -1))
      saveResidues[0] = 0;

    if ((residNumber > 1) && (saveResidues[1] == -1) && (neighboursBlock.size() > 3) &&
        ((neighboursBlock.get(0) != -1) &&
         (neighboursBlock.get(2) != -1) &&
         (neighboursBlock.get(3) != -1)))
      saveResidues[1] = 1;

    if ((residNumber > 1) && (saveResidues[1] == -1) && (neighboursBlock.size() > 2) &&
        ((neighboursBlock.get(0) != -1) &&
         (neighboursBlock.get(2) != -1)))
      saveResidues[1] = 1;

    for (i = 2, num = 0; i < (residNumber - 2); i++, num = 0) {
      if (saveResidues[i] == -1) {
        if (neighboursBlock.get(0) != -1)
          num++;
        if (neighboursBlock.get(1) != -1)
          num++;
        if (neighboursBlock.get(3) != -1)
          num++;
        if (neighboursBlock.get(4) != -1)
          num++;
        saveResidues[i] = (num >= 3) ? i : -1;
      }
      if ((i + 3) < residNumber) {
        neighboursBlock.popFront();
        neighboursBlock.pushBack(saveResidues[i + 3]);
      }
    }

    if ((residNumber > 1) && (saveResidues[(residNumber - 2)] == -1) && (neighboursBlock.size() > 3) &&
        ((neighboursBlock.get(neighboursBlock.size() - 4) != -1) &&
         (neighboursBlock.get(neighboursBlock.size() - 3) != -1) &&
         (neighboursBlock.get(neighboursBlock.size() - 1) != -1)))
      saveResidues[(residNumber - 2)] = (residNumber - 2);

    if ((residNumber > 0) && (saveResidues[(residNumber - 1)] == -1) && (neighboursBlock.size() > 2) &&
        ((neighboursBlock.get(neighboursBlock.size() - 3) != -1) &&
         (neighboursBlock.get(neighboursBlock.size() - 2) != -1)))
      saveResidues[(residNumber - 1)] = (residNumber - 1);

    if (!variable)
      lenBlock = 5;
    else {
      lenBlock = Utils.roundInt(residNumber * 0.01);
      lenBlock = lenBlock > 3 ? (lenBlock < 12 ? lenBlock : 12) : 3;
    }

    blockSize = blockSize > 0 ? blockSize : lenBlock;

    removeSmallerBlocks(blockSize);

    if (terminalGapOnly == true)
      if (!removeOnlyTerminal())
        return null;

    if (complementary == true)
      computeComplementaryAlig(true, false);

    NewValues counter = removeCols_SeqsAllGaps();

    fillNewDataStructure(counter);

    newAlig = new Alignment(filename, aligInfo, counter.matrix, counter.seqsName,
        seqsInfo, counter.sequences, counter.residues, iformat, oformat, shortNames,
        dataType, isAligned, reverse, terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber,
        residNumber, residuesNumber, saveResidues, saveSequences, ghWindow, shWindow,
        blockSize, identities, overlaps);

    return newAlig;
  }

  public Alignment cleanOverlapSeq(float minimumOverlap, float[] overlapSeq,
      boolean complementary) throws IOException {

    Alignment newAlig;
    int i;

    for (i = 0; i < sequenNumber; i++)
      if (overlapSeq[i] < minimumOverlap)
        saveSequences[i] = -1;

    if (complementary == true)
      computeComplementaryAlig(false, true);

    NewValues counter = removeCols_SeqsAllGaps();

    String[] matrixAux = new String[counter.sequences];
    String[] newSeqsName = new String[counter.sequences];
    String[] newSeqsInfo = new String[counter.sequences];
    for (int z = 0; z < counter.sequences; z++) {
      matrixAux[z] = "";
      newSeqsName[z] = "";
      newSeqsInfo[z] = "";
    }

    fillNewDataStructure(matrixAux, newSeqsName, newSeqsInfo);

    newAlig = new Alignment(filename, aligInfo, matrixAux, newSeqsName, newSeqsInfo,
        counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
        isAligned, reverse, terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber, residNumber,
        residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
        identities, overlaps);

    return newAlig;
  }

  public Alignment removeColumns(int[] columns, int init, int size,
      boolean complementary) throws IOException {

    int i, j;

    for (i = init; i < size + init; i += 2)
      for (j = columns[i]; j <= columns[i + 1]; j++)
        saveResidues[j] = -1;

    if (complementary == true)
      computeComplementaryAlig(true, false);

    NewValues counter = removeCols_SeqsAllGaps();

    String[] matrixAux = new String[counter.sequences];
    String[] newSeqsName = new String[counter.sequences];
    String[] newSeqsInfo = new String[counter.sequences];
    for (int z = 0; z < counter.sequences; z++) {
      matrixAux[z] = "";
      newSeqsName[z] = "";
      newSeqsInfo[z] = "";
    }

    fillNewDataStructure(matrixAux, newSeqsName, newSeqsInfo);

    Alignment newAlig = new Alignment(filename, aligInfo, matrixAux, newSeqsName, newSeqsInfo,
        counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
        isAligned, reverse, terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber, residNumber,
        residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
        identities, overlaps);

    return newAlig;
  }

  public Alignment removeSequences(int[] seqs, int init, int size,
      boolean complementary) throws IOException {

    int i, j;

    for (i = init; i < size + init; i += 2)
      for (j = seqs[i]; j <= seqs[i + 1]; j++)
        saveSequences[j] = -1;

    if (complementary == true)
      computeComplementaryAlig(false, true);

    NewValues counter = removeCols_SeqsAllGaps();

    String[] matrixAux = new String[counter.sequences];
    String[] newSeqsName = new String[counter.sequences];
    String[] newSeqsInfo = new String[counter.sequences];
    for (int z = 0; z < counter.sequences; z++) {
      matrixAux[z] = "";
      newSeqsName[z] = "";
      newSeqsInfo[z] = "";
    }

    fillNewDataStructure(matrixAux, newSeqsName, newSeqsInfo);

    Alignment newAlig = new Alignment(filename, aligInfo, matrixAux, newSeqsName, newSeqsInfo,
        counter.sequences, counter.residues, iformat, oformat, shortNames, dataType,
        isAligned, reverse, terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber, residNumber,
        residuesNumber, saveResidues, saveSequences, ghWindow, shWindow, blockSize,
        identities, overlaps);

    return newAlig;
  }

  public void computeComplementaryAlig(boolean residues, boolean sequencesF) {
    int i;

    for (i = 0; i < residNumber && residues; i++)
      saveResidues[i] = (saveResidues[i] == -1) ? i : -1;

    for (i = 0; i < sequenNumber && sequencesF; i++)
      saveSequences[i] = (saveSequences[i] == -1) ? i : -1;
  }

  public boolean removeOnlyTerminal() throws IOException {
    int i;
    int[] gInCol;

    if ((left_boundary == -1) && (right_boundary == -1)) {
      if (calculateGapStats() != true) {
        try (CppOut err = cerr()) {
          err.endl().p("WARNING: Impossible to apply 'terminal-only' method").endl().endl();
        }
        return false;
      }
      gInCol = sgaps.getGapsWindow();

      for (i = 0; i < residNumber && gInCol[i] != 0; i++)
        ;
      left_boundary = i;

      for (i = residNumber - 1; i > -1 && gInCol[i] != 0; i--)
        ;
      right_boundary = i;
    } else if (left_boundary >= right_boundary) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Check your manually set left '").p(left_boundary)
            .p("' and right '").p(right_boundary).p("' boundaries'").endl().endl();
      }
      return false;
    }

    right_boundary += 1;

    for (i = left_boundary; i < right_boundary; i++)
      saveResidues[i] = i;

    return true;
  }

  public void removeSmallerBlocks(int blockSize) {
    int i, j, pos, block;

    if (blockSize == 0)
      return;

    for (i = 0, pos = 0, block = 0; i < residNumber; i++) {
      if (saveResidues[i] != -1)
        block++;
      else {
        if (block < blockSize)
          for (j = pos; j <= i; j++)
            saveResidues[j] = -1;
        pos = i + 1;
        block = 0;
      }
    }

    if (block < blockSize)
      for (j = pos; j < i; j++)
        saveResidues[j] = -1;
  }

  public NewValues removeCols_SeqsAllGaps() throws IOException {
    int i, j, valid, gaps;
    boolean warnings = false;
    NewValues counter = new NewValues();

    for (i = 0, counter.residues = 0; i < residNumber; i++) {
      if (saveResidues[i] == -1)
        continue;

      for (j = 0, valid = 0, gaps = 0; j < sequenNumber; j++) {
        if (saveSequences[j] == -1)
          continue;
        if (Utils.charAt(sequences[j], i) == '-')
          gaps++;
        valid++;
      }
      if (gaps == valid) {
        if (!warnings)
          cerr().endl();
        warnings = true;
        try (CppOut err = cerr()) {
          err.p("WARNING: Removing column '").p(i).p("' composed only by gaps").endl();
        }
        saveResidues[i] = -1;
      } else {
        counter.residues++;
      }
    }

    for (i = 0, counter.sequences = 0; i < sequenNumber; i++) {
      if (saveSequences[i] == -1)
        continue;

      for (j = 0, valid = 0, gaps = 0; j < residNumber; j++) {
        if (saveResidues[j] == -1)
          continue;
        if (Utils.charAt(sequences[i], j) == '-')
          gaps++;
        valid++;
      }
      if (gaps == valid) {
        if (!warnings)
          cerr().endl();
        warnings = true;

        if (keepSequences) {
          try (CppOut err = cerr()) {
            err.p("WARNING: Keeping sequence '").p(seqsName[i])
                .p("' composed only by gaps").endl();
          }
          counter.sequences++;
        } else {
          try (CppOut err = cerr()) {
            err.p("WARNING: Removing sequence '").p(seqsName[i])
                .p("' composed only by gaps").endl();
          }
          saveSequences[i] = -1;
        }
      } else {
        counter.sequences++;
      }
    }
    if (warnings)
      cerr().endl();

    counter.matrix = new String[counter.sequences];
    counter.seqsName = new String[counter.sequences];
    for (int z = 0; z < counter.sequences; z++) {
      counter.matrix[z] = "";
      counter.seqsName[z] = "";
    }

    return counter;
  }

  public void fillNewDataStructure(String[] newMatrix, String[] newNames, String[] newInfos) {
    int i, j, k;

    for (i = 0, j = 0; i < sequenNumber; i++) {
      if (saveSequences[i] == -1)
        continue;

      newNames[j] = seqsName[i];
      // C++ dereferences NULL here for most non-FASTA/PIR trim paths (UB).
      // Use empty metadata only on that undefined path; do not alter strict's
      // separate, observable unfiltered-seqsInfo behavior.
      newInfos[j] = seqsInfo != null ? seqsInfo[i] : "";
      StringBuilder sb = new StringBuilder();
      for (k = 0; k < residNumber; k++) {
        if (saveResidues[k] == -1)
          continue;
        sb.append(Utils.charAt(sequences[i], k));
      }
      newMatrix[j] = sb.toString();
      j++;
    }
  }

  public void fillNewDataStructure(NewValues data) {
    int i, j, k;

    for (i = 0, j = 0; i < sequenNumber; i++) {
      if (saveSequences[i] == -1)
        continue;

      data.seqsName[j] = seqsName[i];
      StringBuilder sb = new StringBuilder();
      for (k = 0; k < residNumber; k++) {
        if (saveResidues[k] == -1)
          continue;
        sb.append(Utils.charAt(sequences[i], k));
      }
      data.matrix[j] = sb.toString();
      j++;
    }
  }

  // ------------------------------------------------------------------ stats

  public boolean calculateGapStats() {
    if (sequences == null)
      return false;

    if (sgaps == null) {
      sgaps = new StatisticsGaps(sequences, sequenNumber, residNumber, dataType);
      sgaps.applyWindow(ghWindow);
    }

    return true;
  }

  public void printStatisticsGapsColumns() throws IOException {
    if (calculateGapStats()) {
      try (CppOut cout = CppOut.toStdout()) {
        sgaps.printGapsColumns(cout);
      }
    }
  }

  public void printStatisticsGapsTotal() throws IOException {
    if (calculateGapStats()) {
      try (CppOut cout = CppOut.toStdout()) {
        sgaps.printGapsAcl(cout);
      }
    }
  }

  public boolean setSimilarityMatrix(SimilarityMatrix sm) throws IOException {
    if (scons == null)
      scons = new StatisticsConservation(sequences, sequenNumber, residNumber, dataType);

    if (!scons.setSimilarityMatrix(sm))
      return false;

    return true;
  }

  public boolean calculateConservationStats() throws IOException {
    if (calculateGapStats() != true)
      return false;

    if (scons == null)
      return false;

    if (scons.isSimMatrixDef() != true)
      return false;

    if (!scons.calculateVectors(sequences, sgaps.getGapsWindow()))
      return false;

    if (scons.isDefinedWindow())
      return true;
    else
      return scons.applyWindow(shWindow);
  }

  public void printStatisticsConservationColumns() throws IOException {
    if (calculateConservationStats()) {
      try (CppOut cout = CppOut.toStdout()) {
        scons.printConservationColumns(cout);
      }
    }
  }

  public void printStatisticsConservationTotal() throws IOException {
    if (calculateConservationStats()) {
      try (CppOut cout = CppOut.toStdout()) {
        scons.printConservationAcl(cout);
      }
    }
  }

  public void printCorrespondence() throws IOException {
    int i;
    try (CppOut cout = CppOut.toStdout()) {
      cout.p("#ColumnsMap\t");
      if (residNumber > 0) {
        for (i = 0; i < residNumber - 1; i++)
          cout.p(saveResidues[i]).p(", ");
        cout.p(saveResidues[i]).endl();
      } else {
        // alignment.cpp:1809 reads saveResidues[0] of a zero-length allocation.
        // Its allocator-dependent integer is UB, not a reproducible column map.
        cout.endl();
      }
    }
  }

  public int getNumSpecies() {
    return sequenNumber;
  }

  public int getNumAminos() {
    return residNumber;
  }

  public void trimTerminalGaps(boolean terminalOnly_, int[] boundaries_) {
    terminalGapOnly = terminalOnly_;
    if (boundaries_ != null) {
      left_boundary = boundaries_[0];
      right_boundary = boundaries_[1];
    }
  }

  public void setWindowsSize(int ghWindow_, int shWindow_) {
    ghWindow = ghWindow_;
    shWindow = shWindow_;
  }

  public void setOutputFormat(int format_, boolean shortNames_) {
    oformat = format_;
    shortNames = shortNames_;
  }

  public void setBlockSize(int blockSize_) {
    blockSize = blockSize_;
  }

  public int getShortNames() {
    return shortNames ? 1 : 0;
  }

  public int getReverse() {
    return reverse ? 1 : 0;
  }

  public void setReverse() {
    reverse = true;
  }

  public void setKeepSequencesFlag(boolean flag) {
    keepSequences = flag;
  }

  public void setKeepSeqsHeaderFlag(boolean flag) {
    keepHeader = flag;
  }

  public int getBlockSize() {
    return blockSize;
  }

  public int getInputFormat() {
    return iformat;
  }

  public int getOutputFormat() {
    return oformat;
  }

  public void getSequencesNames(String[] Names) {
    for (int i = 0; i < sequenNumber; i++)
      Names[i] = seqsName[i];
  }

  public void getSequences(String[] Names, int[] lengths) {
    for (int i = 0; i < sequenNumber; i++) {
      lengths[i] = Utils.removeCharacter('-', sequences[i]).length();
      Names[i] = seqsName[i];
    }
  }

  public void getSequences(String[] Names, String[] Sequences, int[] Lengths) {
    for (int i = 0; i < sequenNumber; i++) {
      Names[i] = seqsName[i];
      Sequences[i] = Utils.removeCharacter('-', sequences[i]);
      Lengths[i] = Sequences[i].length();
    }
  }

  public boolean getSeqNameOrder(String[] names, int[] orderVector) {
    int i, j, numNames;

    for (i = 0, numNames = 0; i < sequenNumber; i++) {
      for (j = 0; j < sequenNumber; j++) {
        if (seqsName[i].equals(names[j])) {
          orderVector[i] = j;
          numNames++;
          break;
        }
      }
    }

    return numNames == sequenNumber;
  }

  public void sequenMatrix() {
    if (seqMatrix == null)
      seqMatrix = new SequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
  }

  public void destroySequenMatrix() {
    seqMatrix = null;
  }

  public void getColumnSeqMatrix(int column, int[] columnSeqMatrix) {
    if (seqMatrix == null)
      seqMatrix = new SequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
    seqMatrix.getColumn(column, columnSeqMatrix);
  }

  public void getColumnSeqMatrix(int value, int row, int[] columnSeqMatrix) {
    if (seqMatrix == null)
      seqMatrix = new SequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
    seqMatrix.getColumn(value, row, columnSeqMatrix);
  }

  public void setSeqMatrixOrder(int[] order) {
    if (seqMatrix == null)
      seqMatrix = new SequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
    seqMatrix.setOrder(order);
  }

  public SequencesMatrix getSeqMatrix() {
    if (seqMatrix == null)
      seqMatrix = new SequencesMatrix(sequences, seqsName, sequenNumber, residNumber);
    return seqMatrix;
  }

  public int getTypeAlignment() {
    if (dataType == 0)
      dataType = Utils.checkTypeAlignment(sequenNumber, residNumber, sequences);
    return dataType;
  }

  public int[] getCorrespResidues() {
    return saveResidues;
  }

  public int[] getCorrespSequences() {
    return saveSequences;
  }

  public boolean isFileAligned() {
    return isAligned;
  }

  // ------------------------------------------------------------------ identities / overlaps

  public void calculateSeqIdentity() throws IOException {
    char indet;

    indet = getTypeAlignment() == Utils.AAType ? 'X' : 'N';

    identities = new float[sequenNumber][];

    // Phase 1: upper triangle. Each unit builds only its own row, so running rows on
    // daemon workers is bit-identical to the serial transcription.
    Parallel.forIndex(0, sequenNumber, row -> {
      float[] values = new float[sequenNumber];
      for (int jj = row + 1; jj < sequenNumber; jj++) {
        int hit = 0, dst = 0;
        for (int k = 0; k < residNumber; k++) {
          char ci = Utils.charAt(sequences[row], k);
          char cj = Utils.charAt(sequences[jj], k);
          if (((ci != indet) && (ci != '-')) ||
              ((cj != indet) && (cj != '-'))) {
            dst++;
            if (ci == cj)
              hit++;
          }
        }
        values[jj] = (float) hit / dst;
      }
      values[row] = 0;
      identities[row] = values;
    });

    // Phase 2: mirror the symmetric half (cheap, restores the original fill order).
    for (int i = 0; i < sequenNumber; i++)
      for (int j = 0; j < i; j++)
        identities[i][j] = identities[j][i];
  }

  public void calculateSeqOverlap() {
    int i, j, k, shared, referenceLength;
    char indet;

    indet = getTypeAlignment() == Utils.AAType ? 'X' : 'N';

    overlaps = new float[sequenNumber][];

    for (i = 0; i < sequenNumber; i++) {
      overlaps[i] = new float[sequenNumber];

      for (j = 0; j < sequenNumber; j++) {
        for (k = 0, shared = 0, referenceLength = 0; k < residNumber; k++) {
          char ci = Utils.charAt(sequences[i], k);
          char cj = Utils.charAt(sequences[j], k);
          if ((ci != indet) && (ci != '-')) {
            referenceLength++;
            if ((cj != indet) && (cj != '-'))
              shared++;
          }
        }
        overlaps[i][j] = (float) shared / referenceLength;
      }
    }
  }

  public int selectMethod() throws IOException {
    float mx, avg, maxSeq = 0, avgSeq = 0;
    int i, j;

    if (identities == null)
      calculateSeqIdentity();

    for (i = 0; i < sequenNumber; i++) {
      for (j = 0, mx = 0, avg = 0; j < sequenNumber; j++) {
        if (i != j) {
          mx = mx < identities[i][j] ? identities[i][j] : mx;
          avg += identities[i][j];
        }
      }
      avgSeq += avg / (sequenNumber - 1);
      maxSeq += mx;
    }

    avgSeq = avgSeq / sequenNumber;
    maxSeq = maxSeq / sequenNumber;

    if (avgSeq >= 0.55)
      return Utils.GAPPYOUT;
    else if (avgSeq <= 0.38)
      return Utils.STRICT;
    else {
      if (sequenNumber <= 20)
        return Utils.GAPPYOUT;
      else {
        if ((maxSeq >= 0.5) && (maxSeq <= 0.65))
          return Utils.GAPPYOUT;
        else
          return Utils.STRICT;
      }
    }
  }

  public void printSeqIdentity() throws IOException {
    int i, j, k, pos, maxLongName;
    float mx, avg, maxAvgSeq = 0, maxSeq = 0, avgSeq = 0;
    float[][] maxs;

    if (identities == null)
      calculateSeqIdentity();

    maxs = new float[sequenNumber][];

    for (i = 0; i < sequenNumber; i++) {
      maxs[i] = new float[2];

      for (k = 0, mx = 0, avg = 0, pos = i; k < sequenNumber; k++) {
        if (i != k) {
          avg += identities[i][k];
          if (mx < identities[i][k]) {
            mx = identities[i][k];
            pos = k;
          }
        }
      }
      avgSeq += avg / (sequenNumber - 1);
      maxAvgSeq += mx;

      maxs[i][0] = mx;
      maxs[i][1] = pos;
    }

    avgSeq = avgSeq / sequenNumber;
    maxAvgSeq = maxAvgSeq / sequenNumber;

    for (i = 0, maxLongName = 0; i < sequenNumber; i++)
      maxLongName = Utils.max(maxLongName, seqsName[i].length());

    try (CppOut cout = CppOut.toStdout()) {
      cout.precision(4);
      cout.fixed();

      for (i = 0, maxSeq = 0; i < sequenNumber; i++)
        if (maxs[i][0] > maxSeq)
          maxSeq = maxs[i][0];

      cout.p("## MaxIdentity\t").p(maxSeq);
      cout.endl().p("#> MaxIdentity\tGet the maximum identity value for any pair ")
          .p("of sequences in the alignment").endl();

      cout.endl().p("## AverageIdentity\t").p(avgSeq);
      cout.endl().p("#> AverageIdentity\tAverage identity between all sequences");

      cout.endl().endl().p("## Identity sequences matrix");
      for (i = 0; i < sequenNumber; i++) {
        cout.endl().setw(maxLongName + 2).left().p(seqsName[i]).p("\t");
        for (j = 0; j < i; j++)
          cout.left().setw(10).p(identities[i][j]).p("\t");
        cout.left().setw(10).p(1.00).p("\t");
        for (j = i + 1; j < sequenNumber; j++)
          cout.left().setw(10).p(identities[i][j]).p("\t");
      }
      cout.endl();

      cout.endl().p("## AverageMostSimilarIdentity\t").p(maxAvgSeq);
      cout.endl().p("#> AverageMostSimilarIdentity\t Average identity between ")
          .p("most similar pair-wise sequences");

      cout.endl().endl().p("## Identity for most similar pair-wise sequences ")
          .p("matrix").endl();
      for (i = 0; i < sequenNumber; i++)
        cout.setw(maxLongName + 2).left().p(seqsName[i])
            .p("\t").left().setw(5).p(maxs[i][0]).p("\t")
            .p(seqsName[(int) maxs[i][1]]).endl();
    }
  }

  public void printSeqOverlap() throws IOException {
    int i, j, k, pos, maxLongName;
    float mx, avg, maxAvgSeq = 0, maxSeq = 0, avgSeq = 0;
    float[][] maxs;

    if (overlaps == null)
      calculateSeqOverlap();

    maxs = new float[sequenNumber][];

    for (i = 0; i < sequenNumber; i++) {
      maxs[i] = new float[2];

      for (k = 0, mx = 0, avg = 0, pos = i; k < sequenNumber; k++) {
        if (i != k) {
          avg += overlaps[i][k];
          if (mx < overlaps[i][k]) {
            mx = overlaps[i][k];
            pos = k;
          }
        }
      }
      avgSeq += avg / (sequenNumber - 1);
      maxAvgSeq += mx;

      maxs[i][0] = mx;
      maxs[i][1] = pos;
    }

    avgSeq = avgSeq / sequenNumber;
    maxAvgSeq = maxAvgSeq / sequenNumber;

    for (i = 0, maxLongName = 0; i < sequenNumber; i++)
      maxLongName = Utils.max(maxLongName, seqsName[i].length());

    try (CppOut cout = CppOut.toStdout()) {
      cout.precision(4);
      cout.fixed();

      for (i = 0, maxSeq = 0; i < sequenNumber; i++)
        if (maxs[i][0] > maxSeq)
          maxSeq = maxs[i][0];

      cout.p("## MaxOverlap\t").p(maxSeq);
      cout.endl().p("#> MaxOverlap\tGet the maximum overlap value for any pair ")
          .p("of sequences in the alignment").endl();

      cout.endl().p("## AverageOverlap\t").p(avgSeq);
      cout.endl().p("#> AverageOverlap\tAverage overlap between all sequences");

      cout.endl().endl().p("## Overlap sequences matrix");
      for (i = 0; i < sequenNumber; i++) {
        cout.endl().setw(maxLongName + 2).left().p(seqsName[i]).p("\t");
        for (j = 0; j < sequenNumber; j++)
          cout.left().setw(10).p(overlaps[i][j]).p("\t");
      }
      cout.endl();
    }
  }

  public boolean calculateSpuriousVector(float overlap, float[] spuriousVector) throws IOException {
    int ovrlap;
    float floatOverlap;
    char indet;

    floatOverlap = overlap * (float) (sequenNumber - 1);
    ovrlap = (int) (overlap * (sequenNumber - 1));

    if (floatOverlap > (float) ovrlap)
      ovrlap++;

    if (spuriousVector == null)
      return false;

    if (getTypeAlignment() == Utils.AAType)
      indet = 'X';
    else
      indet = 'N';

    final int overlapCut = ovrlap;
    final char indeterminate = indet;

    // Each unit owns spuriousVector[i] and only reads the alignment.
    Parallel.forIndex(0, sequenNumber, i -> {
      int seqValue = 0;
      for (int j = 0; j < residNumber; j++) {
        int hit = 0;
        for (int k = 0; k < i; k++) {
          if (sequences[i].charAt(j) == sequences[k].charAt(j))
            hit++;
          else if ((sequences[i].charAt(j) != indeterminate) && (sequences[i].charAt(j) != '-')
              && (sequences[k].charAt(j) != indeterminate) && (sequences[k].charAt(j) != '-'))
            hit++;
        }
        for (int k = (i + 1); k < sequenNumber; k++) {
          if (sequences[i].charAt(j) == sequences[k].charAt(j))
            hit++;
          else if ((sequences[i].charAt(j) != indeterminate) && (sequences[i].charAt(j) != '-')
              && (sequences[k].charAt(j) != indeterminate) && (sequences[k].charAt(j) != '-'))
            hit++;
        }

        if (hit >= overlapCut)
          seqValue++;
      }

      spuriousVector[i] = ((float) seqValue / residNumber);
    });

    return true;
  }

  // ------------------------------------------------------------------ clusters

  int[] represStatic;

  public int[] calculateRepresentativeSeq(float maximumIdent) throws IOException {
    int i, j, pos, clusterNum;
    int[] cluster;
    int[][] seqs;
    float max;

    if (identities == null)
      calculateSeqIdentity();

    seqs = new int[sequenNumber][];
    for (i = 0; i < sequenNumber; i++) {
      seqs[i] = new int[2];
      seqs[i][0] = Utils.removeCharacter('-', sequences[i]).length();
      seqs[i][1] = i;
    }

    Utils.quicksort(seqs, 0, sequenNumber - 1);

    cluster = new int[sequenNumber];
    cluster[0] = seqs[sequenNumber - 1][1];
    clusterNum = 1;

    for (i = sequenNumber - 2; i >= 0; i--) {
      for (j = 0, max = 0, pos = -1; j < clusterNum; j++) {
        if (identities[seqs[i][1]][cluster[j]] > maximumIdent) {
          if (identities[seqs[i][1]][cluster[j]] > max) {
            max = identities[seqs[i][1]][cluster[j]];
            pos = j;
          }
        }
      }

      if (pos == -1) {
        cluster[j] = seqs[i][1];
        clusterNum++;
      }
    }

    represStatic = new int[clusterNum + 1];
    represStatic[0] = clusterNum;
    for (i = 0; i < clusterNum; i++)
      represStatic[i + 1] = cluster[i];

    return represStatic;
  }

  public float getCutPointClusters(int clusterNumber) throws IOException {
    float max, min, avg, gMax, gMin, startingPoint, prevValue = 0, iter = 0;
    int i, j, clusterNum;
    int[] cluster;
    int[][] seqs;

    if (clusterNumber == sequenNumber)
      return 1;
    else if (clusterNumber == 1)
      return 0;

    if (identities == null)
      calculateSeqIdentity();

    for (i = 0, gMax = 0, gMin = 1, startingPoint = 0; i < sequenNumber; i++) {
      for (j = 0, max = 0, avg = 0, min = 1; j < i; j++) {
        if (max < identities[i][j])
          max = identities[i][j];
        if (min > identities[i][j])
          min = identities[i][j];
        avg += identities[i][j];
      }
      for (j = i + 1; j < sequenNumber; j++) {
        if (max < identities[i][j])
          max = identities[i][j];
        if (min > identities[i][j])
          min = identities[i][j];
        avg += identities[i][j];
      }
      startingPoint += avg / (sequenNumber - 1);
      if (max > gMax)
        gMax = max;
      if (min < gMin)
        gMin = min;
    }

    startingPoint /= sequenNumber;

    seqs = new int[sequenNumber][];
    for (i = 0; i < sequenNumber; i++) {
      seqs[i] = new int[2];
      seqs[i][0] = Utils.removeCharacter('-', sequences[i]).length();
      seqs[i][1] = i;
    }
    Utils.quicksort(seqs, 0, sequenNumber - 1);

    cluster = new int[sequenNumber];
    cluster[0] = seqs[sequenNumber - 1][1];

    do {
      clusterNum = 1;
      for (i = sequenNumber - 2; i >= 0; i--) {
        for (j = 0; j < clusterNum; j++)
          if (identities[seqs[i][1]][cluster[j]] > startingPoint)
            break;

        if (j == clusterNum) {
          cluster[j] = seqs[i][1];
          clusterNum++;
        }
      }

      if ((clusterNum == clusterNumber) || (iter > 10))
        break;
      else {
        if (clusterNum > clusterNumber) {
          gMax = startingPoint;
          startingPoint = (gMax + gMin) / 2;
        } else {
          gMin = startingPoint;
          startingPoint = (gMax + gMin) / 2;
        }

        if (prevValue != clusterNum) {
          iter = 0;
          prevValue = clusterNum;
        } else
          iter++;
      }
    } while (true);

    return startingPoint;
  }

  public Alignment getClustering(float identityThreshold) throws IOException {
    int i, j;
    int[] clustering;
    Alignment newAlig;

    clustering = calculateRepresentativeSeq(identityThreshold);

    if (saveSequences == null) {
      // fillMatrices allocates the mapping for aligned input only, but main() lets
      // -clusters/-maxidentity through on unaligned data, where the reference binary
      // segfaults on this exact write. Documented divergence: Java lazily allocates
      // the identity mapping instead of reproducing the crash.
      saveSequences = new int[sequenNumber];
      for (int z = 0; z < sequenNumber; z++)
        saveSequences[z] = z;
    }
    for (i = 0; i < sequenNumber; i++)
      saveSequences[i] = -1;
    for (i = 1; i <= clustering[0]; i++)
      saveSequences[clustering[i]] = clustering[i];

    String[] matrixAux = new String[clustering[0]];
    String[] newSeqsName = new String[clustering[0]];
    String[] newSeqsInfo = new String[clustering[0]];
    for (int z = 0; z < clustering[0]; z++) {
      matrixAux[z] = "";
      newSeqsName[z] = "";
      newSeqsInfo[z] = "";
    }

    for (i = 0, j = 0; i < sequenNumber; i++)
      if (saveSequences[i] != -1) {
        newSeqsName[j] = seqsName[i];
        // alignment.cpp:1515 has the same NULL-seqsInfo UB as manual trimming.
        newSeqsInfo[j] = seqsInfo != null ? seqsInfo[i] : "";
        matrixAux[j] = sequences[i];
        j++;
      }

    newAlig = new Alignment(filename, aligInfo, matrixAux, newSeqsName, newSeqsInfo,
        clustering[0], residNumber, iformat, oformat, shortNames, dataType, isAligned,
        reverse, terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber, residNumber, residuesNumber,
        saveResidues, saveSequences, ghWindow, shWindow, blockSize, identities,
        overlaps);

    return newAlig;
  }

  // ------------------------------------------------------------------ backtranslation

  public Alignment getTranslationCDS(int newResidues, int newSequences, int[] ColumnsToKeep,
      String[] oldSeqsName, SequencesMatrix oldSeqMatrix, Alignment ProtAlig) {

    Alignment newAlig;
    int i, j, k, l, oldResidues;
    int[] mappedSeqs, tmpSequence, selectedRes;

    mappedSeqs = new int[newSequences];
    for (i = 0; i < sequenNumber; i++)
      for (j = 0; j < newSequences; j++)
        if (seqsName[i].equals(oldSeqsName[j])) {
          mappedSeqs[j] = i;
          break;
        }

    oldResidues = oldSeqMatrix.getResidNumber();
    selectedRes = new int[oldResidues];

    for (i = 0; i < oldResidues; i++)
      selectedRes[i] = i;

    if (newResidues == 0) {
      // alignment.cpp:1579/1586 reads ColumnsToKeep[0/-1] for an empty trim (UB).
      // Select nothing safely; successful native cases still emit the same empty
      // output warning, including when -keepseqs retains empty sequence records.
      java.util.Arrays.fill(selectedRes, -1);
    } else {
      for (j = 0; j < ColumnsToKeep[0]; j++)
        selectedRes[j] = -1;

      for (i = 0; i < newResidues - 1; i++)
        for (j = ColumnsToKeep[i] + 1; j < ColumnsToKeep[i + 1]; j++)
          selectedRes[j] = -1;

      for (j = ColumnsToKeep[newResidues - 1] + 1; j < oldResidues; j++)
        selectedRes[j] = -1;
    }

    String[] matrixAux = new String[newSequences];
    StringBuilder[] matAcc = new StringBuilder[newSequences];
    for (int z = 0; z < newSequences; z++) {
      matrixAux[z] = "";
      matAcc[z] = new StringBuilder();
    }
    tmpSequence = new int[oldResidues];

    for (i = 0; i < newSequences; i++)
      if (oldSeqMatrix.getSequence(oldSeqsName[i], tmpSequence)) {
        for (j = 0; j < oldResidues; j++) {
          if ((selectedRes[j] != -1) && (tmpSequence[j] != 0)) {
            for (k = 3 * (tmpSequence[j] - 1), l = 0; l < 3; k++, l++) {
              if (sequences[mappedSeqs[i]].length() > k)
                matAcc[i].append(sequences[mappedSeqs[i]].charAt(k));
              else
                matAcc[i].append('N');
            }
          } else if (selectedRes[j] != -1) {
            matAcc[i].append("---");
          }
        }
        matrixAux[i] = matAcc[i].toString();
      } else {
        return null;
      }

    newAlig = new Alignment(filename, "", matrixAux, oldSeqsName, null, newSequences,
        newResidues * 3, ProtAlig.getInputFormat(), ProtAlig.getOutputFormat(),
        ProtAlig.getShortNames() != 0, Utils.DNAType, true, ProtAlig.getReverse() != 0,
        terminalGapOnly, left_boundary, right_boundary,
        keepSequences, keepHeader, sequenNumber, oldResidues * 3, null, null,
        null, 0, 0, ProtAlig.getBlockSize(), null, null);

    return newAlig;
  }

  public boolean prepareCodingSequence(boolean splitByStopCodon, boolean ignStopCodon,
      Alignment proteinAlig) throws IOException {

    boolean warning = false;
    int i;

    int numbProtSeqs = proteinAlig.getNumSpecies();
    String[] protSeqsNames = new String[numbProtSeqs];
    String[] protSequences = new String[numbProtSeqs];
    int[] protSeqsLengths = new int[numbProtSeqs];

    proteinAlig.getSequences(protSeqsNames, protSequences, protSeqsLengths);

    if (getTypeAlignment() == Utils.AAType) {
      try (CppOut err = cerr()) {
        err.endl().p("ERROR: Check input CDS file. It seems to content protein ")
            .p("residues.").endl().endl();
      }
      return false;
    }

    for (i = 0; i < sequenNumber; i++) {
      int current_prot;
      for (current_prot = 0; current_prot < numbProtSeqs; current_prot++)
        if (protSeqsNames[current_prot].equals(seqsName[i]))
          break;
      if (current_prot == numbProtSeqs)
        continue;

      if (sequences[i].indexOf("-") != -1) {
        if (!warning)
          cerr().endl();
        try (CppOut err = cerr()) {
          err.p("ERROR: Sequence \"").p(seqsName[i]).p("\" has, at least, one gap")
              .endl().endl();
        }
        return false;
      }

      if ((sequences[i].length() % 3) != 0) {
        if (!warning)
          cerr().endl();
        warning = true;
        try (CppOut err = cerr()) {
          err.p("WARNING: Sequence length \"").p(seqsName[i]).p("\" is not ")
              .p("multiple of 3 (length: ").p(sequences[i].length()).p(")").endl();
        }
      }

      if (ignStopCodon)
        continue;

      // stop codon scans: TGA, TAA, TAG (in this order), position multiple of 3
      String[] codons = {"TGA", "TAA", "TAG"};
      for (String codon : codons) {
        int found = -1;
        do {
          found = sequences[i].indexOf(codon, found + 1);

          if ((found != -1) && ((found % 3) == 0)) {
            String protein = protSequences[current_prot];
            int proteinPos = found / 3;
            // C++ string::operator[](size()) yields NUL. Greater indices are UB;
            // safely choose NUL there too, retaining the stop-codon diagnostic.
            char aminoAcid = proteinPos < protein.length() ? toupper(protein.charAt(proteinPos)) : '\0';
            if ((aminoAcid == 'U') || (aminoAcid == 'O') || (aminoAcid == 'X') || (aminoAcid == 'B'))
              continue;

            else if (splitByStopCodon) {
              if (!warning)
                cerr().endl();
              warning = true;
              try (CppOut err = cerr()) {
                err.p("WARNING: Cutting sequence \"").p(seqsName[i]).p("\" at first")
                    .p(" appearance of stop codon \"").p(codon).p("\" (residue \"")
                    .p(aminoAcid).p("\") at position ").p(found + 1).p(" (length: ")
                    .p(sequences[i].length()).p(")").endl();
              }
              sequences[i] = sequences[i].substring(0, found);
            } else {
              if (!warning)
                cerr().endl();
              try (CppOut err = cerr()) {
                err.p("ERROR: Sequence \"").p(seqsName[i]).p("\" has stop codon \"")
                    .p(codon).p("\" (residue \"").p(aminoAcid).p("\") at position ")
                    .p(found + 1).p(" (length: ").p(sequences[i].length()).p(")").endl().endl();
              }
              return false;
            }
          }
        } while (found != -1);
      }
    }

    return true;
  }

  public boolean checkCorrespondence(String[] names, int[] lengths, int totalInputSeqs,
      int multiple) throws IOException {

    int i, j, seqLength, indet;
    boolean warnings = false;
    String tmp;

    for (i = 0; i < sequenNumber; i++) {
      tmp = Utils.removeCharacter('-', sequences[i]);
      seqLength = tmp.length() * multiple;
      indet = (tmp.length() - Utils.min(lastNotOf(tmp, "X"), lastNotOf(tmp, "x"))) - 1;

      for (j = 0; j < totalInputSeqs; j++) {
        if (seqsName[i].equals(names[j])) {
          if (seqLength == lengths[j])
            break;
          else if (seqLength < lengths[j]) {
            if (!warnings)
              cerr().endl();
            warnings = true;
            try (CppOut err = cerr()) {
              err.p("WARNING: Sequence \"").p(seqsName[i]).p("\" will be cut at ")
                  .p("position ").p(seqLength).p(" (length: ").p(lengths[j]).p(")").endl();
            }
            break;
          } else if ((indet > 0) && (indet > (seqLength - lengths[j]) / 3)) {
            if (!warnings)
              cerr().endl();
            warnings = true;
            try (CppOut err = cerr()) {
              err.p("WARNING: Sequence \"").p(seqsName[i]).p("\" has some inde")
                  .p("termination symbols 'X' at the end of sequence. They will be")
                  .p(" included in the final alignment.").endl();
            }
            break;
          } else {
            if (!warnings)
              cerr().endl();
            warnings = true;
            try (CppOut err = cerr()) {
              err.p("WARNING: Sequence \"").p(seqsName[i]).p("\" has less nucleo")
                  .p("tides (").p(lengths[j]).p(") than expected (").p(seqLength)
                  .p("). It will be added N's to complete the sequence").endl();
            }
            break;
          }
        }
      }

      if (j == totalInputSeqs) {
        try (CppOut err = cerr()) {
          err.endl().p("ERROR: Sequence \"").p(seqsName[i]).p("\" is not in ")
              .p("CDS file.").endl().endl();
        }
        return false;
      }
    }

    return true;
  }

  /** C++ string::find_last_not_of(ch): last index whose char != ch, -1 when none. */
  private static int lastNotOf(String s, String ch) {
    for (int i = s.length() - 1; i >= 0; i--)
      if (ch.indexOf(s.charAt(i)) < 0)
        return i;
    return -1;
  }

  // ------------------------------------------------------------------ info

  public void printAlignmentInfo(CppOut file) throws IOException {
    int i, j, valid_res, max, min, max_pos, min_pos, total_res;

    max = 0;
    max_pos = 0;
    min_pos = 0;
    min = residuesNumber[0];

    for (i = 0, total_res = 0; i < sequenNumber; i++) {
      for (j = 0, valid_res = 0; j < residuesNumber[i]; j++)
        valid_res += (sequences[i].charAt(j) != '-' ? 1 : 0);

      total_res += valid_res;

      max_pos = (max > valid_res) ? max_pos : i;
      max = (max > valid_res) ? max : valid_res;
      min_pos = (min < valid_res) ? min_pos : i;
      min = (min < valid_res) ? min : valid_res;
    }

    file.p("## Total sequences\t").p(sequenNumber).endl();
    if (isFileAligned())
      file.p("## Alignment length\t").p(residNumber).endl();
    file.p("## Avg. sequence length\t").p((float) total_res / sequenNumber).endl()
        .p("## Longest seq. name\t'").p(seqsName[max_pos]).p("'").endl()
        .p("## Longest seq. length\t").p(max).endl()
        .p("## Shortest seq. name\t'").p(seqsName[min_pos]).p("'").endl()
        .p("## Shortest seq. length\t").p(min).endl();
  }

  public void calculateColIdentity(float[] ColumnIdentities) {
    int i, j, counter, pos, max, columnLen;
    char letter;
    String column;

    char indet = getTypeAlignment() == Utils.AAType ? 'X' : 'N';
    char gapSymbol = '-';

    for (i = 0, max = 0; i < residNumber; i++, max = 0) {
      column = "";
      for (j = 0; j < sequenNumber; j++) {
        char c = Utils.charAt(sequences[j], i);
        if ((toupper(c) != indet) && (c != gapSymbol))
          column += toupper(c);
      }
      columnLen = column.length();

      StringBuilder col = new StringBuilder(column);
      while (col.length() > 0) {
        letter = col.charAt(0);
        counter = 0;
        pos = 0;

        do {
          counter += 1;
          col.deleteCharAt(pos);
          pos = col.indexOf(String.valueOf(letter), pos);
        } while (pos != -1);

        if (counter > max)
          max = counter;
        if (col.length() < max)
          break;
      }

      if (columnLen != 0)
        ColumnIdentities[i] = (float) max / columnLen;
    }
  }

  public void printColumnsIdentity_DescriptiveStats() throws IOException {
    float[] colIdentities = new float[residNumber];
    float avg, std, max, min;
    int i, positions;

    Utils.initlVect(colIdentities, residNumber, -1);
    calculateColIdentity(colIdentities);

    for (i = 0, max = 0, min = 1, avg = 0, positions = 0; i < residNumber; i++) {
      if (colIdentities[i] != -1) {
        avg += colIdentities[i];
        max = (colIdentities[i] > max) ? colIdentities[i] : max;
        min = (colIdentities[i] < min) ? colIdentities[i] : min;
        positions += 1;
      }
    }
    avg /= positions;

    for (i = 0, std = 0; i < residNumber; i++)
      if (colIdentities[i] != -1)
        // C++ float += pow(float-difference, int) adds in double, then narrows.
        std = (float) (std + Math.pow((colIdentities[i] - avg), 2));
    std = (float) Math.sqrt(std / positions);

    try (CppOut cout = CppOut.toStdout()) {
      cout.p("#maxColIdentity\t").p(max).endl();
      cout.p("#minColIdentity\t").p(min).endl();
      cout.p("#avgColIdentity\t").p(avg).endl();
      cout.p("#stdColIdentity\t").p(std).endl();
    }
  }

  /** newValues struct equivalent. */
  public static final class NewValues {
    public int residues;
    public int sequences;
    public String[] matrix;
    public String[] seqsName;
  }

  /** std::deque<int> restricted to the operations cleanStrict uses (operator[] included). */
  static final class IntDeque {
    private final int[] a;
    private int head = 0;
    private int size = 0;

    IntDeque(int cap) {
      a = new int[Math.max(cap, 1)];
    }

    void pushBack(int v) {
      a[(head + size) % a.length] = v;
      size++;
    }

    void popFront() {
      head = (head + 1) % a.length;
      size--;
    }

    int size() {
      return size;
    }

    int get(int i) {
      return a[(head + i) % a.length];
    }
  }
}
