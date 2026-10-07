# trimAl v1.5.1 / readAl v1.5.1 比对格式 I/O 逐字节考古报告

> 由 Explore 代理于 Phase 0 产出（源码根 `upstream/source/`）。Java 复刻读写器时逐节对照。
> 行号对应当前 upstream 版本。
>
> **实施校正**：这是初期考古报告，并非最终验收证明。实际代码与
> `verify/alignment_audit`、`verify/companions` 的差分结果优先。特别注意
> `ifstream::eof()` 表示失败读触发的 EOF 标志，不是文件指针已到物理末尾；
> Java 加载器已按该差别处理末行和缺最后换行。`utils::quicksort` 的 pivot
> 是末元素（不是中点）。严格保留原版正常输出；空映射、NULL seqsInfo、
> CRLF compareset double-free 等已知 UB 安全分歧必须单列，不计字节一致 PASS。

## 0. 全局常量与格式码（defines.h）

- `DELIMITERS " \t\n"`（空格+tab+换行）、`OTHDELIMITERS " \t\n,:"`、`OTH2DELIMITERS " \n,:;"`（**无 tab** —— NEXUS 解析专用，行内 tab 不切分，要复刻的 quirk）。
- `HTMLBLOCKS 120`、`PHYLIPDISTANCE 10`。
- 格式码：`1`=CLUSTAL，`3`=NBRF/PIR，`8`=FASTA，`11`=PHYLIP3.2（顺序式），`12`=PHYLIP（交错式），`13`=PHYLIP-PAML，`17`=NEXUS，`21`=MEGA 交错，`22`=MEGA 顺序，`99`=-onlyseqs（仅 readAl），`100`=HTML（仅 readAl `-html`）。
- 默认状态（alignment.cpp:66-92 构造函数）：`iformat=0, oformat=0, shortNames=false, forceCaps=false, upperCase=false, lowerCase=false`（三个大小写标志**永远为 false，无 CLI 可改**），`keepSequences=false, keepHeader=false, gapSymbol="-", reverse=false`。

## 1. 输入格式与自动检测

### 1.1 分发
`alignment::loadAlignment`（alignment.cpp:503-534）：先 `iformat = formatInputAlignment(file)`（rwAlignment.cpp:144-274），`oformat = iformat`（509 行），然后 switch 分发到 8 个 loader（1→loadClustal、3→loadNBRF_Pir、8→loadFasta、11→loadPhylip3_2、12→loadPhylip、17→loadNexus、21→loadMegaInterleaved、22→loadMegaNonInterleaved），**default 直接 return false，无任何格式相关报错文案**。

### 1.2 formatInputAlignment 检测算法（rwAlignment.cpp:144-274）
1. 打开文件，`utils::checkFile`（utils.cpp:391-409：存在且 size>0，否则 false）。
2. 反复 `utils::readLine` 直到拿到非 NULL 行；若 EOF 仍无 → `return false`（163-164）。
3. `firstWord = strtok(line, OTHDELIMITERS)`（167）。
4. 按顺序判定（170-206）：
   - `firstWord == "CLUSTAL"` 或 `"clustal"`（strcmp 精确匹配，**"Clustal" 不匹配**）→ 1
   - `firstWord[0]=='>' && firstWord[3]==';'` → 3（NBRF/PIR）。**quirk**：FASTA 名字前 4 字符形如 `>xx;` 时误判为 PIR，例如 `>gi;123` → 误判为 NBRF。
   - `firstWord[0]=='>'` → 8（FASTA）
   - `"#NEXUS"`/`"#nexus"` → 17
   - `"#MEGA"`/`"#mega"` → 逐字节扫描统计以 `#` 开头的行块数（191-202，逻辑怪异：跳到第一个 `#`，然后逐行读、行首读一字符判 `#` 计数）；`blocks>0`→21（交错），否则 22（顺序）（205）。
   - 其它一律进 phylip 分支（209-264）：`sequenNumber=atoi(firstWord)`、下一 token `residNumber=atoi`；`sequenNumber==1 && residNumber!=0` → 12；否则读第 2 行数 token（blocks++），再读第 3 行数 token（blocks--），`blocks==0`→12（交错），否则 11（顺序）；两数任一为 0 → format 保持 0。

### 1.3 检测失败的报错文案（原文）
- readAl（readAl.cpp:82-84）/ trimAl（main.cpp:121-122）：`ERROR: Alignment not loaded: "<infile>" Check the file's content.`（前后各一个空行）

### 1.4 fillMatrices 校验（rwAlignment.cpp:43-113，所有 loader 收尾调用）
- 每条序列每个字符必须 `isalpha || ispunct`，否则（57-61，原文）：`ERROR: The sequence "<seqsName[i]>" has an unknown (<char>) character.`（前后空行）。**注意**：`>` `;` `*` `.` 等标点全部合法；非法仅指控制字符/空格类。
- 长度不齐 → `isAligned=false`；要求对齐但不齐时（72-76）：`ERROR: Sequences should be aligned (all with same length) and there are not. Check your input alignment`。
- aligned 且某条长度 ≠ residNumber（83-90）：`ERROR: The sequence "<name>" (<len>) does not have the same number of residues fixed by the alignment (<residNumber>).`

### 1.5 各 loader 要点（名字/序列切分规则）
- **readLine**（utils.cpp:411-453）：逐字符读到 `\n`，`\r` 被丢弃；只剥**行首**空格与 tab；空行返回 NULL。**行尾空格/tab 保留**。
- **loadFastaAlignment**（673-770）：先数 `>` 行得 sequenNumber（DELIMITERS 切分后首 token 以 `>` 判定）；重读时 `seqsInfo[i+1]` 存**原始 header**（去掉 `>` 的整行，含空格，736）；名字 = OTHDELIMITERS 首 token 再跳过 `>` 前缀（739-750）——名字在空格/tab/`,`/`:` 处截断；序列 token 用 DELIMITERS 拼接。`fillMatrices(false)`（769）——允许未比对输入。
- **loadClustalAlignment**（498-671）：跳过首行（整行存入 `aligInfo`，583）；先扫第一块数序列数（551-561：一行内无任何 `isalpha` 或 `-` 字符 → 认为是保守性得分行，第一块结束）；再 round-robin 读：名字=OTHDELIMITERS 首 token，序列=其余 token 拼接（642-651）。
- **loadPhylipAlignment**（276-381，格式 12，注释写 "Sequential" 但实际按**交错 round-robin** 读）：首行两个 atoi；其后每行 DELIMITERS 首 token 当名字（仅前 sequenNumber 行）、其余 token 全并入当前序列（336-344）；后续文件按 i 循环追加（349-372）。
- **loadPhylip3_2Alignment**（383-496，格式 11，注释写 "Interleaved" 但实际是**顺序式**读：一条序列累计到 `size()==residNumber` 即换下一条，483-486）：名字行 OTHDELIMITERS 切分；首行块数 `blocksFirstLine` 与后续行比较，不符则把误认的新序列并回上一条（471-477）。
- **loadNBRF_PirAlignment**（1297-1419）：`>P1;NAME` 行 → `strtok(line, ">;")` 第一段（如 `P1`）存 seqsInfo，第二段存 seqsName（1372-1377）；下一行整行描述追加到 seqsInfo（1382-1385）；序列行 OTHDELIMITERS 拼接；行尾 `*` 结束一条序列且从序列尾部删除 `*`（1392-1406）。
- **loadNexusAlignment**（772-913）：大写化行首 token 匹配 `BEGIN`/`MATRIX`/`FORMAT`/`DIMENSIONS`（813-837，DIMENSIONS 用 `=;` 再切）；`[..]` 注释剔除（863-874）；矩阵行用 **OTH2DELIMITERS**（881-893，无 tab）round-robin。
- **loadMegaNonInterleavedAlignment**（915-1112）/ **loadMegaInterleavedAlignment**（1114-1295）：`!TITLE` 行可整体替换 `filename`（975-980/1170-1175，若原行不以 `!` 开头则补一个 `!` 前缀），`!FORMAT` 行追加进 aligInfo；`utils::trimLine`（utils.cpp:455-536）剔除 `"..."` 与 `[...]` 注释（不配对时报错：`ERROR: Possible (") mismatch for comments`、`ERROR: Brackets (]) mismatch found`、`ERROR: Brackets ([) mismatch found`）；名字 token 用 `" #\n"`，序列 token 用 `" \n"`。
- 所有 loader 均在开头追加 `filename = "!Title <file>;"`（如 289-291、396-398、512-514、686-688、786-788、930-932、1129-1131、1311-1313）——该字符串**直接成为 MEGA 输出的第二行**。

## 2. 输出 writer 逐字节规格

写出分发：`printAlignment()`（stdout，alignment.cpp:556-606）与 `saveAlignment()`（文件，alignment.cpp:611-690）。两者在 `sequences==NULL` 返回 false；`residNumber==0 || sequenNumber==0` 时打印（564-566 / 620-621，原文）：`WARNING: Output alignment has not been generated. It is empty`（print 版末尾无句号，save 版有句号 "It is empty."），然后返回 **true**。
`saveAlignment` 对格式 1/11/12/13/17/21/22 预检 isAligned（628-636）：`ERROR: Sequences are not aligned. Output format is not compatible with unaligned sequences.`（无 endl 后缀，缺换行是原样 bug）。文件打开失败 → false（640-641）。写完 `file.close()`，无额外 flush。

全部行结束符为 `"\n"`。

### 2.1 FASTA — `alignmentFastaToFile`（rwAlignment.cpp:1620-1663）
- `tmpMatrix[i]` = 原序或 `utils::getReverse`（纯逆序，非互补）。
- `maxLongName`（1638-1650）：`!keepHeader` → max(所有 `seqsName[i].size()`)；`keepHeader && seqsInfo!=NULL` → max(所有 `seqsInfo[i].size()`)；若 `shortNames && maxLongName>10` → `maxLongName=10`，且 keepHeader 时向 **stderr** 输出：`WARNING: Original sequence header will be cut by character 10`（两段拼接："…cut by charac"+"ter 10"）。
- 每条序列：`">" + (keepHeader && seqsInfo!=NULL ? seqsInfo[i].substr(0,maxLongName) : seqsName[i].substr(0,maxLongName)) + "\n"`；随后按 `residuesNumber[i]` 每 60 残基一行 `substr(j,60) + "\n"`（j 从 0，步长 60）。**序列之间与文件末尾均无空行**；空序列只输出 header 行。
- **quirk**：substr(0,maxLongName) 对短名无影响（maxLongName 是全体最大值）；真正截断只发生在 shortNames。无 keepHeader/seqsInfo==NULL 时回落 seqsName（1653-1656）。

### 2.2 FASTA 变体（-onlyseqs，readAl only）— `getSequences(ostream&)`（rwAlignment.cpp:115-142）
- 先 `utils::removeCharacter('-', seq)` 删除**所有** `-`（reverse 时先 `getReverse` 再删，129-130）。
- 每条：`">" + seqsName[i] + "\n"` + 每 60 残基一行 + **每条后额外空行 `"\n"`**（137）。`.` 不删。

### 2.3 CLUSTAL — `alignmentClustalToFile`（1573-1618）
- `maxLongName`（1597-1598）：`!shortNames` 时 = max(name 长度)；**shortNames 时循环被 `&&(!shortNames)` 短路，maxLongName 恒 0**（quirk：`_m10` 不截断名字，只把列宽缩到 5）。
- 头（1601-1604）：`aligInfo` 非空且 `iformat==oformat` → 输出 `aligInfo + "\n\n"`（即把输入 clustal 文件的原 title 行原样回写）；否则 `"CLUSTAL multiple sequence alignment" + "\n\n"`。
- 体（1609-1614）：外层 j 从 0 步长 60；内层每条：`setw(maxLongName+5) << left << seqsName[i]`（**不截断**，短名右侧补空格至宽度 maxLongName+5）+ `tmpMatrix[i].substr(j,60) + "\n"`；每块后 `"\n\n"`（两个空行，包括最后一块之后——文件以 `"\n\n\n"` 结束：最后一条残基行的 endl + 两个 endl）。

### 2.4 NBRF/PIR — `alignmentNBRF_PirToFile`（1805-1859）
- `alg_datatype`（1822-1827）：DNA/DNADeg→`"DL"`，RNA/RNADeg→`"RL"`，AAType→`"P1"`；**dataType==0 时无 else 分支，alg_datatype 为空串**（quirk）。
- 每条序列（1830-1854）：
  - 若 `seqsInfo!=NULL && iformat==oformat`（PIR→PIR 直通，**与 -keepheader 无关**）：`">" + seqsInfo[i].substr(0,2) + ";" + seqsName[i] + "\n" + seqsInfo[i].substr(2) + "\n"`（第 2 行 = 原描述行）。
  - 否则：`">" + alg_datatype + ";" + seqsName[i] + "\n"` + 第二行 `seqsName[i] + " " + residuesNumber[i] + " bases" + "\n"`（如 `seq1 123 bases`）。
  - 序列体：j 步长 50，k 步长 10，每组输出 `" " + substr(k,10)`——**每个 10 残基组前有一个空格**（行首也有）。
  - 末组之后 `if(k >= residuesNumber[i])`：若 `len%50==0` → `"\n" + " "`（`*` 单独成行且前有一空格）；否则若 `len%10==0` → `" "`（`*` 接在行尾空格后）；然后恒输出 `"*"`；再 `"\n"`。len%10!=0 时 `*` 紧贴最后一个残基。
  - 每条序列后空行 `"\n"`（1854）。

### 2.5 MEGA — `alignmentMegaToFile`（1749-1803）（格式 21 与 22 共用同一 writer）
- `#MEGA` + `"\n"` + `filename` + `"\n"`（1776）。filename 即 `!Title <inputfile>;` 或 MEGA 输入 TITLE 行（见 1.5）。**v1.5.1 没有 `#mev` 头，也没有 `EX`/`END-` 块**。
- 格式行（1779-1788）：`"!Format DataType=DNA "` / `"!Format DataType=RNA "` / `"!Format DataType=protein "`（**尾随一个空格**）+ `"NSeqs=" + sequenNumber + " Nsites=" + residNumber + " indel=- CodeTable=Standard;" + "\n\n"`。
- 每条（1791-1798）：`"#" + seqsName[i] + "\n"`；j 步长 50，k 步长 10：`substr(k,10) + " "`（**每组后空格，行尾有尾随空格**），块尾 `"\n"`；每条序列后空行 `"\n"`。文件以最后一条的空行结束。

### 2.6 NEXUS — `alignmentNexusToFile`（1665-1747）
- `maxLongName` 规则同 CLUSTAL（1689-1690，shortNames→0）。
- 先把 `aligInfo` 中所有 `;` 删除（1696-1697）。
- 头（1700-1701）：`"#NEXUS\n"` + `"BEGIN DATA;\n"` + `" DIMENSIONS NTAX=" + sequenNumber + " NCHAR=" + residNumber + ";\n"`（` DIMENSIONS` **有一个前导空格**）。
- FORMAT 行（1704-1725）：DNA/DNADeg→`"FORMAT DATATYPE=DNA INTERLEAVE=yes GAP=-"`，RNA/DNADeg→`"FORMAT DATATYPE=RNA INTERLEAVE=yes GAP=-"`，AA→`"FORMAT DATATYPE=PROTEIN INTERLEAVE=yes GAP=-"`；dataType==0 时什么都不输出（quirk）。随后扫描 `aligInfo`（1713-1724）：以空格分割，词前缀匹配（`compare(0,7,"MISSING")==0` 或 `substr(i,j).compare(0,7,"missing")==0`）→ 追加 `" " + substr(i, j-i)`；`MATCHCHAR`/`matchchar` 前缀 9 同理。**逐字节 quirk**：MISSING 第一比较用 `substr(i, j-i)`、第二比较用 `substr(i, j)`；MATCHCHAR 两处比较都用 `substr(i, j)`；输出统一 `substr(i, j-i)`；且**最后一个词（后面无空格）永远不被扫描**。最后输出 `";" + "\n"`。
- 注释行（1728-1730）：每条 `[Name: ` + `setw(maxLongName+4) left seqsName[i]`（不截断）+ `"Len: " + residNumber + "]" + "\n"`；随后 `"\n" + "MATRIX" + "\n"`（1731）。
- 矩阵（1734-1742）：j 步长 50；每条 `setw(maxLongName+4) left seqsName[i]`；k 步长 10：`" " + substr(k,10)`；行尾 `"\n"`；每块后空行 `"\n"`。
- 尾（1743）：`";" + "\n" + "END;" + "\n"`。

### 2.7 PHYLIP（交错式，outformat 12）— `alignmentPhylipToFile`（1421-1470）
- 未对齐 → stderr：`ERROR: Sequences are not aligned. Format (PHYLIP) not compatible with unaligned sequences.`（前后空行，函数直接 return，**stdout 版不写任何东西**）。
- `maxLongName = 10 (PHYLIPDISTANCE)`；`!shortNames` 时取 `max(10, 所有 name 长度)`（1446-1448）。**与 clustal/nexus 不同：这里 m10 有效截断**。
- 第一行：`" " + sequenNumber + " " + residNumber + "\n"`（**前导空格**）（1452）。
- 第一块（1455-1457）：每条 `setw(maxLongName+3) << left << seqsName[i].substr(0, maxLongName)` + `tmpMatrix[i].substr(0,60)` + `"\n"`（名字**截断到 maxLongName**，右侧补空格至宽度 maxLongName+3，即最短补 3 个空格）。
- 后续块（1461-1465）：i 从 60 步长 60，每条 `tmpMatrix[j].substr(i,60) + "\n"`（无名字、无缩进、残基连续 60 列**无 10 分组**），块后空行 `"\n"`。
- 文件末尾再多一个空行（1466 `file << endl;`）→ 结尾字节 `"\n\n\n"`（最后残基行 endl + 块 endl + 额外 endl）。residNumber≤60 时输出为：头行、名块、空行、空行。

### 2.8 PHYLIP 3.2（顺序式，outformat 11）— `alignmentPhylip3_2ToFile`（1472-1527）
- 头行同 2.7（1503）。maxLongName 规则同 2.7（1497-1499）。
- 每条（1508-1523）：`setw(maxLongName+3) left seqsName[i].substr(0,maxLongName)`；j 步长 50，k 步长 10：`substr(k,10) + " "`（**每组后空格、行尾尾随空格**），行尾 `"\n"`；若 `(j+50) < residNumber` → 下一行先输出 `setw(maxLongName+3) << " "`（**setw 默认右对齐 + 空串 = 恰好 maxLongName+3 个空格的缩进行**，1519）。
- 每条序列后空行 `"\n"`（1522）。文件以最后一条的空行结束。
- **重大 bug（必须复刻）**：1514 用的是 `sequences[i]` 而非 `tmpMatrix[i]` —— `-reverse` 时此格式**不反转**序列。

### 2.9 PHYLIP PAML（outformat 13）— `alignmentPhylip_PamlToFile`（1529-1571）
- 头行同上（1560）。maxLongName 规则同上（1554-1556）。
- 每条**单行**（1564-1566）：`setw(maxLongName+3) left seqsName[i].substr(0,maxLongName)` + **整条序列** + `"\n"`。
- 全部之后一个空行 `"\n"`（1567）。
- **同款 bug**：1566 用 `sequences[i]`，`-reverse` 不生效。

### 2.10 各 writer 一览表

| writer | 行号 | 首行 | 名字宽 setw | 名字截断 | 残基/行 | 组间空格 | 行尾尾随空格 | 末尾 |
|---|---|---|---|---|---|---|---|---|
| fasta | 1620 | `>name` | 无 | substr(0,maxLong) | 60 | 无 | 无 | 末残基行\n |
| onlyseqs | 115 | `>name` | 无 | 无 | 60（去-） | 无 | 无 | 每条后空行 |
| clustal | 1573 | title 或 `CLUSTAL multiple sequence alignment` | maxLong+5 left | 无 | 60 | 无 | 无 | 每块后两空行 |
| nbrf | 1805 | `>XX;name` | 无 | 无 | 50（10 组） | 组**前**空格 | 无（%10==0 时组后补空格接 `*`） | 每条后空行 |
| mega | 1749 | `#MEGA` | 无 | 无 | 50（10 组） | 组**后**空格 | **有** | 每条后空行 |
| nexus | 1665 | `#NEXUS` | maxLong+4 left | 无 | 50（10 组） | 组**前**空格 | 无 | `;\nEND;\n` |
| phylip(12) | 1421 | ` N R` | maxLong+3 left | substr(0,maxLong) | 60（无组） | 无 | 无 | 额外空行（\n\n\n） |
| phylip3.2(11) | 1472 | ` N R` | maxLong+3 left | substr(0,maxLong) | 50（10 组） | 组后空格 | **有** | 每条后空行 |
| paml(13) | 1529 | ` N R` | maxLong+3 left | substr(0,maxLong) | 整条一行 | 无 | 无 | 文件尾一个空行 |

setw 均为 C++ `<iomanip>`：`left` 显式指定；不指定时 setw **默认右对齐**——1519 的 `setw(maxLongName+3) << " "` 就是靠右对齐空串产生纯空格缩进。setw 从不截断，截断只来自 substr。

## 3. 序列名处理汇总

- **无任何非法字符替换、无大小写转换**；名字保留输入原样（除各 loader 的 token 截断）。
- token 截断规则因 loader 而异（见 1.5）：fasta/clustal/phylip11 名字止于空格/tab/`,`/`:`；phylip12 止于空格/tab；nexus 止于空格/换行/`,`/`:`/`;`（**tab 例外地不切**）；nbrf 止于 `>`/`;`；mega 止于空格/`#`/换行。
- **-keepheader**：仅两处生效——(a) FASTA writer 输出 `seqsInfo`（原始 header）而非 `seqsName`（rwAlignment.cpp:1640-1656）；(b) alignment.cpp:1896-1898 只是把 flag 存进对象并在 trimAl 的 alignment 拷贝构造中传递（alignment.cpp:170）。CLI 只在 trimAl（main.cpp:530-532），readAl 没有。注意 NBRF writer 的 header 直通**不看** keepHeader，只看 `seqsInfo!=NULL && iformat==oformat`（1833）。
- **重名**：全代码库无重名检测。重名时 `sequencesMatrix::getSequence`（sequencesMatrix.cpp:181-201）只返回**第一个**匹配行。
- `-fasta_m10/-phylip_m10/-phylip_paml_m10/-phylip3.2_m10`（readAl.cpp:125-174 / main.cpp:164-205）：`outformat` 同基础格式 + `shortNames=true` → fasta 截 10、phylip 三兄弟截 10；clustal/nexus 的 `_m10` 行为见 2.3/2.6（不截断、宽 5/4）。

## 4. 残基字符处理

- 输入不 toupper/tolower，**原样保留大小写**（`forceCaps/upperCase/lowerCase` 成员永远 false，alignment.cpp:70-72）。
- 合法性 = `isalpha || ispunct`（C locale）。`-` `.` `?` `*` 均为 punct 合法。
- **gap 只有 `'-'`**（`gapSymbol="-"`，alignment.cpp:82）；`'.'` 被当作普通残基保留与输出。**不存在 `.`↔`-` 互转**；唯一删 gap 处是 `getSequences`（只删 `-`）与统计函数中的 `removeCharacter('-')`。
- PIR：loader 删除末尾 `*`（rwAlignment.cpp:1392-1406）；writer 每条末尾补 `*`（1850）。
- 反序：`utils::getReverse`（utils.cpp:542-553）纯字符逆序；只有 readAl `-reverse`（readAl.cpp:109-110, 293-294 → `setReverse` alignment.cpp:1883-1885）触发；writer 生效矩阵见第 2 节（11/13 有 bug 不生效）。
- HTML 着色：`utils::determineColor`（utils.cpp:722-832）：`toupper` 后 G→`'o'`、P→`'y'`；`'-'` 直接返回 `'w'`；其余按 CLUSTAL 规则 switch，阈值判断用 `utils::lookForPattern`（utils.cpp:835-851，`toupper` 计数比例 ≥ level，**gap `-` 不计入 dataset 命中但计入 column.size() 分母**）；未命中任何 case → `'w'`。

## 5. HTML 输出

### 5.1 `-html`（readAl，outformat 100）— `alignmentColourHTML`（rwAlignment.cpp:2137-2215）
骨架（2160-2172，逐字节）：
```
<!DOCTYPE html>
<html><head>
    <meta http-equiv="Content-Type" content="text/html;charset=ISO-8859-1" />
    <title>readAl v1.5.1</title>
    <style type="text/css">
    #b  { background-color: #3366ff; }
    #r  { background-color: #cc0000; }
    #g  { background-color: #33cc00; }
    #p  { background-color: #ff6666; }
    #m  { background-color: #cc33cc; }
    #o  { background-color: #ff9900; }
    #c  { background-color: #46C7C7; }
    #y  { background-color: #FFFF00; }
    </style>
  </head>

  <body>
  <pre>
```
（`</style>`、`</head>`、`<body>`、`<pre>` 行是用 `\n` 字面量输出的，且 `  <body>\n  <pre>` 前导两空格；style 无 `media` 属性；title 是 `readAl v1.5.1`。）
- 每块 120 残基（2176-2177）：先 `endl`；列号行 `setw(maxLongName+19) right (j+10)`，然后 i=j+20 步长 10、≤min(residNumber,upper)：`setw(10) right i`（2181-2183）；分隔行 `endl + setw(maxLongName+10)`（右对齐消耗在**第一个** `+`/`=` 上 → 前导 maxLongName+9 个空格），i=j+1 起，`i%10==0`→`'+'` 否则 `'='`（2186-2188）。
- 每条序列：`endl + setw(maxLongName+9) left seqsName[i]`（不截断），随后残基：`'w'` 原字符，否则 `"<span id=" + type + ">" + res + "</span>"`（2194-2205）。
- 尾（2212）：`    </pre>\n  </body>\n</html>\n`。
- 未对齐 → stderr `ERROR: Sequences are not aligned.`（2148-2151），返回 false。

### 5.2 `-htmlout`（trimAl）— `alignmentSummaryHTML`（1861-2135）
调用点：main.cpp:1741-1746（对 **origAlig** 调用，传 singleAlig 的 `getNumAminos()/getNumSpecies()/getCorrespResidues()/getCorrespSequences()` 与 compareVect；失败打印 `ERROR: It's imposible to generate the HTML output file.`）。`getCorrespResidues/Sequences` 即 `saveResidues/saveSequences`（alignment.cpp:2042-2055）。
- `minHTML = max(25, maxLongName+10)`（1895）。
- 头（1923-1955）：同 5.1 骨架但 title=`trimAl v1.5.1 Summary`、style 带 `media="all"`；除 8 个 `#b/r/g/p/m/o/c/y` 外还有：
  - `.sel  { background-color: #B9B9B9; }`、`.nsel { background-color: #E9E9E9; }`
  - **c1..c12 十六进制**（1941-1952）：c1 `#FFFBF2`、c2 `#FFF8CC`、c3 `#FAF0BE`、c4 `#F0EAD6`、c5 `#F3E5AB`、c6 `#F4C430`、c7 `#C2B280`+`color: white`、c8 `#DAA520`+white、c9 `#B8860B`+white、c10 `#918151`+white、c11 `#967117`+white、c12 `#6E5411`+white。**css 里每行以 `\n` 结尾、前导 4 空格、类名后空格数固定**（如 `.c1   {`、`.c10  {`）。
- 选中统计行（1958-1962）：`    <span class=sel>Selected Sequences: ` + `setw(5) right` seqs + ` /Selected Residues: ` + `setw(7) right` residues + `</span>`；第二行 `    <span class=nsel>Deleted Sequences:  ` + `setw(5) right` (sequenNumber-seqs) + ` /Deleted Residues:  ` + `setw(7) right` (residNumber-residues) + `</span>`（注意 "Sequences:" 后 1 空格 vs "Deleted Sequences:" 后 2 空格）。
- 图例行（1965-1990，仅对应向量非 NULL 才打）：`endl + setw(minHTML) left 标签`，标签恰 20 字符（`"    Gaps Scores:        "` / `"    Similarity Scores:  "` / `"    Consistency Scores: "`），后接 12 个 span：c1..c9 形如 `<span  class=cN>  =0=  </span>`（**class 前两个空格**），c10..c12 形如 `<span class=c10> <.750 </span>`（一个空格）；span 文本固定：Gaps=`=0=,<.001,<.050,<.100,<.150,<.200,<.250,<.350,<.500,<.750,<1.00,=1=`；Similarity=`=0=,<1e-6,<1e-5,<1e-4,<.001,<.010,<.100,<.250,<.500,<.750,<1.00,=1=`；Consistency 同 Gaps。
- **怪异条件**（1992-1993）：`if ((gapsValues != NULL) or (simValues == NULL) or (consValues == NULL)) file << endl;` —— 仅当 gapsValues==NULL 且 sim/cons 都非 NULL 时不加换行。
- 块循环（1996-2124，步长 120）：
  - 列号行（2000-2002）：`endl + setw(minHTML+10) right (j+10)`；i=j+20 步长 10：`setw(10) right i`。
  - 分隔行（2005-2008）：`endl + setw(minHTML+1) right` 后循环输出 `+`/`=`（setw 消耗在首字符 → 行首 minHTML 个空格）。
  - 序列行（2011-2027）：`    <span class=` + (`sel>`|`nsel>`) + `seqsName[i]` + `</span>` + `setw(minHTML-4-len) right ""`（=minHTML-4-len 个空格）；残基同 5.1 span 规则；行尾 `endl`。
  - `    Selected Cols:      ` 行（2029-2032）：`endl + setw(minHTML) left` 标签，每列 `"<span class=" + (res[k]?"sel":"nsel") + "> </span>"`。
  - 三个分数行（2039-2122）阈值：
    - Gaps（gapsValues）：`==0`→c12（全 gap）；`==sequenNumber`→c1（无 gap）；`1-(float)v/n` ≥.750→c11、≥.500→c10、≥.350→c9、≥.250→c8、≥.200→c7、≥.150→c6、≥.100→c5、≥.050→c4、≥.001→c3、否则 c2。
    - Similarity（simValues=float* mdkw）：`==1`→c12、`==0`→c1、≥.750→c11、≥.500→c10、≥.250→c9、≥.100→c8、≥.010→c7、≥.001→c6、≥1e-4→c5、≥1e-5→c4、≥1e-6→c3、否则 c2。
    - Consistency（consValues）：同 Gaps 阈值，但 `==1`→c12、`==0`→c1。
    - span 模板同图例（c1..c9 双空格、c10..c12 单空格；行首标签同图例，2040/2068/2096）。块末 `file << endl`（2123）。
- 数据源：`sgaps->getGapsWindow()`（statisticsGaps.cpp:183-186，int*）与 `scons->getMdkwVector()`（statisticsConservation.cpp:269，float*）；仅当统计对象已建才非 NULL。
- 尾（2127）：`    </pre>\n  </body>\n</html>\n`；写 `ofstream`，打开失败返回 false（1885-1887）。

## 6. statisticsFiles 与 sequencesMatrix

### 6.1 statisticsFiles（-sfc/-sft 背后的残基计数矩阵壳）
statisticsFiles.h:37-49 / statisticsFiles.cpp:31-58：极简类，只有 `int **sequencesMatrix; int columns; int columnLength;`。默认构造清零；带参构造按 `species × aminos` 分配 `new int` 矩阵但**从不初始化、从不填充、无任何方法**；析构释放。它实际上是个**半成品/死代码壳**（真正的 -sfc/-sft 统计在 compareFiles.cpp，通过 main.cpp:1600-1605 调 printStatisticsFileColumns/printStatisticsFileAcl）。无格式识别功能。

### 6.2 sequencesMatrix 完整数据结构（sequencesMatrix.h:47-117 / .cpp）
- 字段：`int resNumber; int seqsNumber; int **matrix; string *seqsName;`
- 构造（.cpp:41-72）：复制名字数组；`matrix[i][j]` 按列扫描赋值：非 `'-'` 残基依次编号 1,2,3,...（k 递增），`'-'` 记 0（63-70）。**即"该序列中第几个非 gap 残基"**（用于 -backtrans 密码子对应）。
- `operator=`（74-99）：**bug**：第 94 行 `matrix[i][j] = matrix[i][j]`——自赋值（调用方实际未用该运算符路径）。
- `printMatrix`（120-132）：每 20 列一组，`setw(4) << matrix[j][k] << " "`（每值后一空格），每序列一行 `endl`，每组后额外 `endl`。
- `getColumn(int column, int*)`（134-145）：越界列填 0。`getColumn(int value, int row, int*)`（147-160）：在 row 中找首个==value 的列，取该列全员；找不到填 -1。
- `setOrder(int *order)`（162-179）：按 order 重排行。
- `getSequence(string, int*)`（181-201）：按名查行（第一个匹配），拷贝整行。
- 头文件声明了 `removeColumns(int,int,int*,int*)`（.h:109）但 **.cpp 中没有实现**——死声明。

## 7. 排序调用点与浮点/整数输出精度

### 7.1 排序（全部为 `utils::quicksort`）
| 调用点 | 所在函数 | 排序对象 |
|---|---|---|
| alignment.cpp:865 | cleanCompareFile | float vectAux（一致性值副本），取 baseline 分位 |
| alignment.cpp:1048 | cleanCombMethods | float vectAux（保守值），找 20%/80% 点 |
| alignment.cpp:1272 | calculateRepresentativeSeq | int**（[去gap长度, 序号]），按长度排 |
| alignment.cpp:1390 | getCutPointClusters | int** 同上 |
| alignment.cpp:2107 | cleanByCutValue(double cut,...) | int* gInCol 副本 |
| alignment.cpp:2442-2443 | cleanByCutValue(双阈值) | float MDK 副本 + int gInCol 副本 |
| statisticsConservation.cpp:324 | calcCutPoint | float MDK_Window 副本 |
| statisticsConservation.cpp:396 | printConservationAcl | float 副本（输出前排序） |
| compareFiles.cpp:467 | printStatisticsFileAcl | float 副本 |

### 7.2 浮点/整数输出点与精度
- **`cout` 默认精度 6（defaultfloat，%g 风格）**：`printAlignmentInfo` 的 `## Avg. sequence length\t`（rwAlignment.cpp:2250，`(float)total/sequenNumber`）。
- **`cout.precision(4) + cout << fixed`**：`printSeqIdentity`（autAlignment.cpp:255-256）与 `printSeqOverlap`（335-336）：`## MaxIdentity/## AverageIdentity/## MaxOverlap/## AverageOverlap` 为 4 位定点小数；矩阵值 `setiosflags(ios::left) << setw(10) << 值 << "\t"`（identities 下三角含对角 1.00；overlap 全矩阵），most-similar 矩阵值 `setw(5)`（271-289, 351-353）。
- **`cout.precision(10)`（无 fixed → 10 位有效数字、可能科学计数法）**：gaps 每列/累计、similarity 每列/累计、compareFiles 统计（见各 writer 章节）。
- **注意 setw 不截断**：值超过宽度时整段输出。
- 其它：similarityMatrix.cpp:425 `cerr << setw(8) << setprecision(4) << right`；sequencesMatrix.cpp:126 `setw(4)`。
- `printCorrespondence`（alignment.cpp:1802-1811）：`#ColumnsMap\t` + `saveResidues[i] << ", "`（除最后一个，末元素无逗号）。

## 8. 必须复刻的怪异行为/边界 bug 清单

1. **phylip3.2(11) 与 paml(13) 的 `-reverse` 失效**：writer 体用 `sequences[i]` 而非已反转的 `tmpMatrix[i]`（rwAlignment.cpp:1514、1566）。outformat 12 正常用 tmpMatrix（1463）。
2. **`_m10` 对 clustal/nexus 不截断名字**：`maxLongName` 循环被 `&&(!shortNames)` 短路为 0（1597-1598、1689-1690），列宽塌缩为 5/4，长名破坏对齐——原样复刻。
3. **格式检测 quirk**：`>xx;` 前缀的 FASTA 误判为 PIR（174）；只有全大写 `CLUSTAL`/全小写 `clustal` 可识别（170）；任意非空首词都进 phylip 分支，靠 atoi 失败兜底返回 0。
4. **PIR 直通条件**：`seqsInfo!=NULL && iformat==oformat` 时自动保留原 `>XX;name`+描述行，**与 -keepheader 无关**（1833-1835）；dataType 未识别时 `alg_datatype` 为空串 → 输出 `">;name"`。
5. **NEXUS FORMAT 行组装**：MISSING/MATCHCHAR 比较 substr 上下界不对称（1715-1721）；无尾空格的最后一个词被忽略；dataType==0 时输出裸 `";"` 行；aligInfo 中 `;` 全删。
6. **MEGA 检测的分块计数逻辑**（191-202）本身古怪（先读掉 `#` 行首字符、行内循环条件为 `c != '\n'`），按代码逐行移植即可。
7. **readLine 只剥行首空白、保留行尾空白、丢弃 `\r`**（utils.cpp:427-442）→ 行尾空格会进入 aligInfo/名字等（clustal title 行）。
8. **setw 一次性语义**：`alignmentSummaryHTML` 分隔行 `setw(minHTML+1) << right` 只填充**第一个** `+`/`=`（2005-2007）；2186 `setw(maxLongName+10)` 未写对齐——继承 2181 的 `right`。
9. **PIR 末尾 `*` 的三种排布**（1845-1851）：`len%50==0` → `*` 独占一行且行首一空格；`len%10==0` → 行尾空格后接 `*`；其余紧贴。len==0 时循环不执行、**无 `*`**。
10. **空序列**：fasta/onlyseqs 只输出 header 行（for 循环空转）；onlyseqs 删 `-` 后可能为空。
11. **saveAlignment 未对齐预检报错行缺 endl**（alignment.cpp:632-634，cerr 文本后无 endl）。
12. **空比对输出告警后返回 true**（"It is empty"/"It is empty." 句号差异，alignment.cpp:564-566 vs 620-622）。
13. **main.cpp `if(stats < -1) cout << endl` 死代码**（1559、1566 等）：stats 从 0 递增，条件永假，无输出影响。
14. **-htmlout 作用于 origAlig** 且在 `-out` 写文件之前（main.cpp:1741-1746 vs 1763）；`getCorrespResidues` 在未 trim 时是 0..n-1 恒等表。
15. **sequencesMatrix::operator= 自赋值 bug**（sequencesMatrix.cpp:94）。
16. **OTH2DELIMITERS 无 tab**：NEXUS 行内 tab 会残留在 token 中进序列。
17. `loadPhylip3_2Alignment:451` `firstLine = 1`（恒真重赋值，无害但移植时别"优化"）。
18. **readAl 只在 `argc==1`/`-h`/`--version` 位于 argv[1] 时生效**（readAl.cpp:51-68）；trimAl 同（main.cpp:76-95）。

## 9. readAl（readAl.cpp）与 trimAl（main.cpp）功能差异清单

共同点：两者共用同一 `alignment` 类与全部 loader/writer；`-in/-out` 语义一致；`loadAlignment` 失败文案相同；格式码一致。

**readAl 独有（trimAl 没有）**：
- `-format`（输入格式 + 是否对齐，readAl.cpp:96-97, 239-264；输出 `## Input filename\t'..'`、`## Input file format\t..`、`## Input file aligned\tYES|NO`）
- `-type`（readAl.cpp:100-101, 266-282；`## Input file datatype\t` + `nucleotides:dna|nucleotides:dna_degenerate_codes|nucleotides:rna|nucleotides:rna_degenerate_codes|amino-acids|unknown`）
- `-info`（readAl.cpp:105-106, 284-285 → `printAlignmentInfo`）
- `-reverse`（readAl.cpp:109-110, 293-294）
- `-onlyseqs`（outformat 99，readAl.cpp:183-184）
- `-html`（outformat 100，readAl.cpp:179-180；trimAl 的 HTML 是 `-htmlout` 走 alignmentSummaryHTML，readAl 无 -htmlout）
- readAl 的 stdout 直通：`-out` 缺省时 `inAlig.printAlignment()`（readAl.cpp:304-306）
- readAl 参数互斥检查（202-230）：必须 -in；必须有动作；无 -out 时信息选项与输出格式二选一；仅有信息选项时不得带 -out。

**trimAl 独有（readAl 没有）**：
- `-htmlout <file>`（main.cpp:145-149, 1741-1746）
- `-keepheader`（main.cpp:530-532, 1497）、`-keepseqs`（525-527, 1496）
- 全部修剪方法：`-nogaps/-noallgaps/-gappyout/-strict/-strictplus/-automated1`、手动阈值 `-gt/-gapthreshold、-st/-simthreshold、-ct/-conthreshold、-cons`、`-selectcols/-selectseqs {..}`、`-block`、`-terminalonly`、`--set_boundaries`、`-resoverlap/-seqoverlap`、`-clusters/-maxidentity`、`-complementary`
- 统计输出：`-sgc/-sgt/-ssc/-sst/-sfc/-sft/-sident/-soverlap`（main.cpp:1556-1596）
- 比较集：`-compareset/-forceselect/-ct` 体系（232-272）、`-matrix`、`--alternative_matrix`
- CDS：`-backtrans/-ignorestopcodon/-splitbystopcodon`
- `-w/-gw/-sw/-cw` 窗口、`-colnumbering`（1126, 1774-1775 → printCorresponding）
- `-h` 会额外打印 examples（main.cpp:86-89, 1923-1974）
- trimAl **没有** `-html/-onlyseqs/-reverse/-format/-type/-info`。

## 附：Java 逐字节复刻的三大技术要点

1. **C++ ostream 格式化语义**：`setw` 不截断、一次性；`left/right` 操纵符粘滞（每处 setw 几乎都紧跟显式对齐符，唯一例外 2186 `setw(maxLongName+10)` 继承前文 `right`）。浮点：`precision(10)` 无 fixed = 10 位有效数字 %g；`fixed+precision(4)` = 4 位定点；默认 6 位 %g。
2. **`\n` vs `endl`**：所有写出均为 `endl`（= `\n` + flush），字节层面就是 `\n`。
3. **`string::substr` 越界**：`substr(pos,len)` len 超长则取到末尾；`tmpMatrix[i].substr(0,60)` 对短序列安全截短——Java `substring` 需自行 clamp 防异常。
