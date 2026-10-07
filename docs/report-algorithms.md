# trimAl v1.5.1 算法考古报告 —— 统计与剪裁

> 由 Explore 代理于 Phase 0 产出（源码根 `upstream/source/`，行号已核对）。
> Java 复刻时逐节对照；怪异行为清单（§9）是字节级一致的关键。
>
> **数值 oracle 校正（实施后）**：本报告是初期源码梳理，不是最终验收依据。
> `exp(-Q[i])` 和 `log10(first20Point/last80Point)` 在参考编译环境选择的是
> **float 重载 `expf/log10f`**，不是 double 计算后 cast；Java 已按 glibc 2.35
> 参考实现复刻，并由 `verify/numeric/math.results.json` 做位级校准。
> `statisticsGaps::calcCutPoint` 的 `(double) columnLength * gapThreshold` 首先
> 提升 columnLength 到 double，因此乘法是 double；下文旧注释不能替代实际代码。
> `vectAux[-1]`、空列 slope、AA 替代矩阵越界等 UB 的观测与限制详见
> `verify/numeric/ub.results.json` 及项目 README 的差异清单；非可移植 C++ 语义。
> 正式验证以各套件保存的原始输出与最终报告为准。

## 0. 全局约定与关键事实（先读）

- 数据类型常量：`DNAType=1, RNAType=2, AAType=3, DNADeg=4, RNADeg=5`（defines.h:39-43）；`GAPPYOUT=1, STRICT=2`（defines.h:48-49）；`HTMLBLOCKS=120`、`PHYLIPDISTANCE=10`（defines.h:55-56）。
- 不确定字符：AA 用 `'X'`，核酸用 `'N'`（statisticsGaps.cpp:50-53、statisticsConservation.cpp:129、autAlignment.cpp:41 等）。
- **本版本不存在 `-colscore`、`-automated2`、`-mw` 选项**（全库 grep 无命中）。`statAl.cpp` 另有 `-scolidentt`（列一致性描述统计，autAlignment.cpp:411-446）。`calcCutPointMixSlope`（statisticsGaps.cpp:242-321）是**死代码**，无任何调用点。
- 所有排序都是 `utils::quicksort`（自写快排，**pivot 取末元素**），无 std::sort/qsort（utils.cpp:230-383）。
- 浮点精度：全部 `float`（simMat/distMat/MDK/identities），求值顺序处多为 float 与 double 混算，下面逐一标注。

## 1. 列统计计算链

### 1.1 statisticsGaps（gap 统计）

**构造函数**（statisticsGaps.cpp:39-82）：
- `columnLength = species`（序列数），`columns = aminos`（列数），`maxGaps=0, halfWindow=0`。
- 逐列统计：`alignmentMatrix[j][i] == '-'` → `gapsInColumn[i]++`；`== indet`（X/N）→ `aminosXInColumn[i]++`（statisticsGaps.cpp:69-75）。
- 直方图：`numColumnsWithGaps[gapsInColumn[i]]++`（长度 `species+1`）；`gapsWindow[i] = gapsInColumn[i]`（**窗口=0 时 gapsWindow 就是原始计数**）；更新 maxGaps（statisticsGaps.cpp:78-80）。

**applyWindow(_halfWindow)**（statisticsGaps.cpp:132-174）：
- 拒绝条件：`(halfWindow == _halfWindow) || (_halfWindow > columns/4)` → 返回 false（**静默忽略**，无报错；`columns/4` 是 int 整除）。构造后 halfWindow=0，所以 `applyWindow(0)` 直接返回 false，保持原始计数（statisticsGaps.cpp:139-140）。
- 窗口求和（**环形反射边界**）：

```cpp
for(j = i - halfWindow, gapsWindow[i] = 0; j <= i + halfWindow; j++) {
  if(j < 0)                gapsWindow[i] += gapsInColumn[-j];
  else if(j >= columns)    gapsWindow[i] += gapsInColumn[((2 * columns - j) - 2)];
  else                     gapsWindow[i] += gapsInColumn[j];
}
```

（statisticsGaps.cpp:155-162。注意左端映射是 `-j`（不是 `-j-1`），右端是 `2*columns-j-2`。例如 halfWindow=2、i=0：j=-2→列2, j=-1→列1, j=0→列0, j=1→列1, j=2→列2。）
- **"中值"实为均值四舍五入**：`gapsWindow[i] = utils::roundInt(((double) gapsWindow[i]/window))`，其中 `window = 2*halfWindow+1`（statisticsGaps.cpp:150,165）。`utils::roundInt(x) = (int)((double)x + 0.5)`（utils.cpp:99-102）——正数即"half-up"；无偶数窗口取中值规则（本版没有真 median；`-mw` 也不存在）。
- 然后 `numColumnsWithGaps[gapsWindow[i]]++`，更新 maxGaps（statisticsGaps.cpp:167-171）。
- 调用点：alignment.cpp:1668（`calculateGapStats` 内，仅首次创建 sgaps 时 `sgaps->applyWindow(ghWindow)`）。

**calcCutPoint(minInputAlignment, gapThreshold)**（statisticsGaps.cpp:188-230，返回 double）：

```cpp
cuttingPoint_gapThreshold = (double) columnLength * gapThreshold;   // columnLength(int) * gapThreshold(float) → float 乘，再转 double
cuttingPoint_MinimumConserv = utils::roundInt(((double)(columns * minInputAlignment) / 100.0));
// 注意: columns(int) * minInputAlignment(float) 先得 float，再 cast double，再 /100.0 (double)
if(cuttingPoint_MinimumConserv > columns) cuttingPoint_MinimumConserv = columns;
for(i = 0, acum = 0; i < columnLength; i++) {
  acum += numColumnsWithGaps[i];
  if(acum >= cuttingPoint_MinimumConserv) break;   // acum(int) >= double 比较
}
if(numColumnsWithGaps[i])
  cuttingPoint_MinimumConserv = (double)(i - ((float)(acum - cuttingPoint_MinimumConserv)/numColumnsWithGaps[i]));
  // i 是 int; (acum - cp) int 减 double → double, cast float; float/int → float; i - float → float; 最后转 double
else cuttingPoint_MinimumConserv = 0;
return utils::max(cuttingPoint_MinimumConserv, cuttingPoint_gapThreshold);   // double max
```

- 语义：取两者较大值（保留更多列的一侧，"-cons 保底优先"）。
- 边界坑：若直方图累计从未达标（如所有列为满 gap 且 baseline=100%），循环退出时 `i == columnLength`，访问 `numColumnsWithGaps[columnLength]`（数组大小 columnLength+1，合法但通常是 0→走 else 分支置 0；若该桶非 0 则产生 `columnLength - 小数` 的怪值）。

**calcCutPoint2ndSlope()**（statisticsGaps.cpp:332-397，gappyout 用）：

```cpp
float maxSlope = -1, *secondSlopeVector;  int prev, pprev, maxIter, act = 0, max = 0;
secondSlopeVector = new float[maxGaps+1];
utils::initlVect(secondSlopeVector, maxGaps, -1.0);   // 只初始化 0..maxGaps-1; 下标 maxGaps 未初始化!
maxIter = maxGaps + 1;
// 前置扫描: 找第一个非空桶之前的索引
for(act = 0, max = 0; numColumnsWithGaps[act] == 0; act++) max = act + 1;
act = 0;
while(act < maxIter) {
  while((numColumnsWithGaps[act]) == 0) act++;      // 无界推进(见bug列表#9.4)
  pprev = act;  if((act+1) >= maxIter) break;
  do { act++; } while((numColumnsWithGaps[act]) == 0);
  prev = act;   if((act+1) >= maxIter) break;
  do { act++; } while((numColumnsWithGaps[act]) == 0);
  if(act >= maxIter) break;
  secondSlopeVector[act] = ((float)(act - pprev) / columnLength);
  secondSlopeVector[act] /= ((float)(numColumnsWithGaps[act] + numColumnsWithGaps[prev]) / columns);
  if(secondSlopeVector[pprev] != -1.0) {
    if((secondSlopeVector[act]/secondSlopeVector[pprev]) > maxSlope) {
      maxSlope = (secondSlopeVector[act]/secondSlopeVector[pprev]);
      max = pprev;
    }
  } else if(secondSlopeVector[prev] != -1.0) {
    if((secondSlopeVector[act]/secondSlopeVector[prev]) > maxSlope) {
      maxSlope = (secondSlopeVector[act]/secondSlopeVector[prev]);
      max = pprev;     // 注意: 这里 max 也赋 pprev (不是 prev) —— 原样复刻
    }
  }
  act = prev;
}
return max;
```

（statisticsGaps.cpp:374-389 为精确转录。`act - pprev` 是 int 相减再 cast float；除法顺序：先 `(float)Δ/columnLength`，再 `/((float)(n[act]+n[prev])/columns)`。）

**calcCutPointMixSlope()**（statisticsGaps.cpp:242-321）：死代码但完整存在（一阶斜率比 + 二阶斜率比取 delta 最大者）。如需完整复刻语义可参考，但 v1.5.1 行为不包含它。

**打印**（见 §6）。

### 1.2 statisticsConservation（保守性统计）

**matrixIdentity（序列两两相异度矩阵）**（statisticsConservation.cpp:123-163）：
- 对每对 `i<j`、每列 k：若 `seq[i][k]` 非 `-` 非 indet：若 `seq[j][k]` 同样有效且相等 → `sum++`；且 `length++`（只要 i 有效 length 就 +1，无论 j 是否有效）。若 i 无效但 j 有效 → `length++`。
- `matrixIdentity[j][i] = (100.0 - ((float) sum/ length) * 100.0); matrixIdentity[i][j] = matrixIdentity[j][i];`
  - **若 length==0（两序列全 gap/X）→ 0/0 = NaN**（bug，见 §9）。表达式顺序：`(float)sum / length`（float/int→float），`* 100.0`（float*double→double），`100.0 - double` → double，再赋给 float（二次舍入）。

**calculateVectors(alignmentMatrix, gaps)**（statisticsConservation.cpp:172-213）——每列 Q/MDK：

```cpp
for(i = 0; i < columns; i++) {
  for(j = 0, num = 0, den = 0; j < sequences; j++) {
    if((alignmentMatrix[j][i] != '-') && (alignmentMatrix[j][i] != indet))
      for(k = j + 1; k < sequences; k++)
        if((alignmentMatrix[k][i] != '-') && (alignmentMatrix[k][i] != indet)) {
          num += matrixIdentity[j][k] * simMatrix->getDistance(alignmentMatrix[j][i], alignmentMatrix[k][i]);
          den += matrixIdentity[j][k];
        }
  }
  Q[i] = (den == 0) ? 0 : num / den;          // float/float
  MDK[i] = (float) exp(-Q[i]);                // exp 以 double 计算(-Q 提升为 double), 再截断为 float
  if(gaps != NULL)
    if(((float) gaps[i] / sequences) >= 0.8) MDK[i] = 0;   // gaps 是【窗口化后】的 getGapsWindow(); 0.8 为 double 字面量
  if(MDK[i] > 1) MDK[i] = 1;
}
```

- 关键：getDistance 返回的是**欧氏距离矩阵 distMat**（不是相似度分数！见 §3）；Q = 以 identity 为权的平均距离；MDK = exp(-Q)。80% gap 判定用窗口化后的 gap 数（int）cast float 除以 sequences。
- `gaps` 实参来自 alignment.cpp:1757（`sgaps->getGapsWindow()`），即已含 -gw 窗口效果。
- 累加顺序：j 外层升序、k 从 j+1 升序——`num += identity*dist` 的浮点累加顺序必须逐对一致才能位级一致。

**applyWindow（保守性）**（statisticsConservation.cpp:223-248）：
- 拒绝条件同 gaps：`(halfWindow == _halfWindow) || (_halfWindow > columns/4)`；构造函数把 halfWindow 初始化为 **-1**（statisticsConservation.cpp:45），因此 shWindow=0 也会真正执行（window=1，`MDK_Window[i] = MDK[i] / (float)window` = 除以 1.0f，位级不变）。
- 环形反射求和后 `MDK_Window[i] = MDK_Window[i] / (float) window;`（statisticsConservation.cpp:239-245）。**注意 MDK_Window 不清零地复用**：构造时 initlVect 清 0（statisticsConservation.cpp:55），重复调用 applyWindow 会累加（实际主流程每对象只调用一次）。
- `isDefinedWindow()` 返回 `halfWindow != -1`（statisticsConservation.cpp:257-260）；alignment.cpp:1763-1766 的逻辑：若已定义则跳过，否则 applyWindow(shWindow)。

**calcCutPoint(minInputAlignment, similThreshold)**（statisticsConservation.cpp:306-346，返回 double）：

```cpp
utils::copyVect(MDK_Window, vectAux, columns);
utils::quicksort(vectAux, 0, columns-1);
for(i = columns - 1; i >= 0; i--) if(vectAux[i] < similThreshold) break;
cuttingPoint_SimilThreshold = vectAux[i];      // BUG: 若所有值 >= 阈值, i == -1 → vectAux[-1] 越界读(见§9.1)
highestPos = (int) ((double)(columns - 1) * (100.0 - minInputAlignment)/100.0);
highestPos = highestPos < (columns - 1) ? highestPos : columns - 1;
cuttingPoint_MinimumConserv = vectAux[highestPos];
return (cuttingPoint_MinimumConserv < cuttingPoint_SimilThreshold ? cuttingPoint_MinimumConserv : cuttingPoint_SimilThreshold);
```

（`(100.0 - minInputAlignment)`：double - float → double；`(int)` 截断。）

**printConservationColumns / printConservationAcl**：见 §6。

### 1.3 序列级 identity/overlap（autAlignment.cpp）

- `calculateSeqIdentity`（autAlignment.cpp:35-74）：dst = "至少一方有效（非 `-`、非 indet）"的列数；hit = dst 内两字符相同的列数（**两方都无效的列不计入 dst，因此不会 hit**）；`identities[i][j] = (float) hit/dst`（dst 可为 0 → NaN/inf）。对称矩阵，`identities[i][i]=0`。
- `calculateSeqOverlap`（autAlignment.cpp:76-109）：以 i 为参考，`overlaps[i][j] = (float) shared/referenceLength`，**非对称**；referenceLength=0 → NaN。
- `calculateRelaxedSeqIdentity`（autAlignment.cpp:111-141）：hit 直接比全长（`/residNumber`），当前无调用点（grep 确认）。
- `calculateColIdentity`（autAlignment.cpp:363-409，statAl `-scolidentt` 用）：每列取非 gap 非 indet 字符（toupper），在可变 string 上反复 `erase+find` 计数最高频残基 `max`，`ColumnIdentities[i] = float(max)/columnLen`；全 gap/X 列保持初值 -1。提前终止条件 `if((int) column.size() < max) break;`。
- `printColumnsIdentity_DescriptiveStats`（autAlignment.cpp:411-446）：avg、max、min 在单循环内累计；`avg /= positions`；`std = sqrt(std_sum/positions)`（`pow((colIdentities[i] - avg), 2)` 是 double 运算）；`positions==0` 时除零。打印无精度设置（默认 6 位有效数字）。

## 2. 剪裁方法全清单

### 2.0 公共底层

**saveResidues/saveSequences 语义**：载入时 `saveResidues[i]=i, saveSequences[i]=i`（rwAlignment.cpp:99-108）；删除标记为 `-1`。

**cleanByCutValue（gap 版，double cut）**（alignment.cpp:2073-2231）：
1. 保留 `gInCol[i] <= cut` 的列，否则 `saveResidues[i] = -1`；计数 counter.residues（:2084-2086）。
2. `oth = utils::roundInt((((baseLine/100.0) - (float) counter.residues/residNumber)) * residNumber)`（:2093）——求值顺序：`(baseLine/100.0)`（float/double→double）；`((float)counter.residues/residNumber)`（float/int→float）；double - float → double；`* residNumber`（double*int→double）；roundInt。
3. 若 oth>0：`counter.residues += oth`；对 gInCol 副本快排，**重设 cut**：`cut = vectAux[(int)((float)(residNumber - 1) * (baseLine)/100.0)]`（:2107-2108；注意这里乘的是 `baseLine` 不是 `100-baseLine`，因为 gap 越小越好，取升序数组的前 baseline% 处）。
4. 放宽回收循环（:2118-2166）：`for(k = utils::roundInt(0.005 * residNumber); (k >= 0) && (oth > 0); k--)`——块大小从 0.5% 长度递减到 0；内层从 `i = residNumber/2, j = i+1` 同时向两侧扫（`i--, j++`），先左后右；每侧先跳过已保留列找到"连续 -1 块头"，若块长 `>= k` 则在块内回收 `gInCol[jn] <= cut` 的列（`saveResidues[jn] = jn; oth--;`），遇不满足立即 break。**i = jn / j = jn 回写**（:2142,2162）。
5. `-block`（blockSize != 0）：把小于 blockSize 的连续保留块整体置 -1，末块单独检查（:2169-2192）。
6. `terminalGapOnly` → `removeOnlyTerminal()`（§2.11）；`complementary` → `computeComplementaryAlig(true,false)`（翻转列选择，alignment.cpp:2889-2897）；`removeCols_SeqsAllGaps()`（§7）；`fillNewDataStructure`；构造新 alignment（把 `saveResidues/saveSequences` 压缩数组传给拷贝构造，alignment.cpp:209-226）。

**cleanByCutValue（保守性版，float cut）**（alignment.cpp:2240-2385）：
1. 保留 `ValueVect[i] > cut`（**严格大于**，注意与 gap 版 `<=` 相反）（:2251-2253）。
2. oth 同式（:2260）。
3. **不重算 cut**（没有快排分支）——放宽时只回收 `ValueVect[jn] == cut`（**严格等于 cut 的列**，:2286,2309）。这是与 gap 版最大的差异，务必原样复刻。
4. 余下与 gap 版相同（blockSize/terminal/complementary/removeCols/fill/construct）。

**cleanByCutValue（组合版）**（alignment.cpp:2398-2580）：
1. 保留 `(MDK_Win[i] > cutCons) && (gInCol[i] <= cutGaps)`（:2411-2413）。
2. oth>0 时：双副本快排，`blCons = vectAuxCons[(int)((float)(residNumber-1)*(100.0-baseLine)/100.0)]`（高位），`blGaps = vectAuxGaps[(int)((float)(residNumber-1)*(baseLine)/100.0)]`（低位）（:2439-2446）。
3. 放宽回收条件是 **OR**：`(MDK_Win[jn] >= blCons) || (gInCol[jn] <= blGaps)`（:2480,2504）。
4. 余下同上。

**removeSmallerBlocks(blockSize)**（alignment.cpp:2946-2975）：与 cleanByCutValue 内联版逻辑一致（`block < blockSize` 即删 `[pos, i]`；末尾 `for(j = pos; j < i; j++)`，此时 i == residNumber）。

**removeOnlyTerminal()**（alignment.cpp:2903-2942）：
- 若 `left_boundary == -1 && right_boundary == -1`（自动模式）：用 `sgaps->getGapsWindow()` 找第一个/最后一个 `gInCol == 0` 的列作为边界（:2918-2922）。
- 手动边界（--set_boundaries）时若 `left >= right` 报错（:2925-2929）。
- **`right_boundary += 1`**（:2933）——注意这会修改成员变量（多次调用会持续 +1，但每个 alignment 对象通常只走一遍）。然后把 `[left, right)` 全部 `saveResidues[i] = i` 恢复（:2938-2939）。

**removeCols_SeqsAllGaps()**（alignment.cpp:2980-3049）：见 §7。

### 2.1 手动 -gt（cleanGaps，main.cpp:1661-1662 / alignment.cpp:706-734）
- `-gt 0.9` → 内部 `gapThreshold = 1 - 0.9 = 0.1`（main.cpp:307，浮点减法！0.9 的 float 表示被 atof 成 double 再 1-x，gapThreshold 是 float，有表示误差——Java 需用同样 float 存储复刻）。
- `cleanGaps(baseLine=conserve, gapsPct=gapThreshold, complementary)`：先 calculateGapStats，`cut = sgaps->calcCutPoint(baseLine, gapsPct)`（double），`cleanByCutValue(cut, baseLine, getGapsWindow(), complementary)`。
- `-nogaps` = `cleanGaps(0, 0, complementary)`（main.cpp:1615）→ cut=0 → 只留 gap 数为 0 的列。

### 2.2 -st（cleanConservation，main.cpp:1652-1657 / alignment.cpp:750-778）
- `cut = (float) scons->calcCutPoint(baseLine, conservationPct)`（**double→float 显式转换**，alignment.cpp:764）→ cleanByCutValue(float 版)。

### 2.3 -gt + -st 组合（clean，alignment.cpp:799-836）
- `cutGaps = sgaps->calcCutPoint(...)`（double）；`cutCons = scons->calcCutPoint(...)`（赋给 float，隐式转换，:822）→ cleanByCutValue(组合版)。

### 2.4 -cons（保底百分比）
- 只是传入 baseLine；`conserve == -1` 时归 0（main.cpp:1491-1492）；范围校验 0..100（main.cpp:408）。

### 2.5 -nogap/-noallgaps
- `-noallgaps` = `cleanNoAllGaps(complementary)`（alignment.cpp:1102-1123）：`cleanByCutValue((sequenNumber - 1), 0, getGapsWindow(), complementarity)`——cut=序列数-1，baseLine=0。

### 2.6 -gappyout（clean2ndSlope，alignment.cpp:950-978）
- `cut = sgaps->calcCutPoint2ndSlope()`（int）→ `cleanByCutValue(cut, 0, getGapsWindow(), complementarity)`（baseLine=0，无放宽——oth=roundInt((0 - ratio)*n)，可能为负，循环不执行）。

### 2.7 -strict / -strictplus（cleanCombMethods，alignment.cpp:986-1093 + cleanStrict:2592-2739）

cleanCombMethods：
1. calculateGapStats；`gapCut = sgaps->calcCutPoint2ndSlope()`；`gaps = getGapsWindow()`。
2. calculateConservationStats；`simil = scons->getMdkwVector()`。
3. `positions[i] = i` 当 `gaps[i] <= gapCut`，统计 `acm`（:1029-1034）。
4. `vectAux[0..acm-1]` = 被选列的 simil 值（顺序原样）；`utils::quicksort(vectAux, 0, acm-1)`（升序，:1048）。
5. 20%/80% 点（:1057-1062）：

```cpp
for(i = acm - 1, j = 1; i >= 0; i--, j++) {
  if((((float) j/acm) * 100.0) <= 20.0) first20Point = vectAux[i];
  if((((float) j/acm) * 100.0) <= 80.0) last80Point  = vectAux[i];
}
```

由于 j 递增、条件在尾部成立，最终 first20Point = 满足 `j <= 0.2*acm` 的最大 j 对应值（高端第 floor(0.2*acm) 个值），last80Point 同理（floor(0.8*acm)）。
6. simCut（:1069-1072）：

```cpp
inic = log10(first20Point);      // double (C log10)
fin  = log10(last80Point);
vlr  = ((inic - fin) / 10) + fin;  // 除以 int 10 → double
simCut = (float) pow(10, vlr);
```

若 first20Point 或 last80Point 为 0（列全被 80%gap 规则清零）→ log10(0) = -inf → simCut = 0（pow(10,-inf)）。NaN 场景也需复刻语义。
7. `cleanStrict(gapCut, getGapsWindow(), simCut, getMdkwVector(), complementarity, variable)`；`-strict` → variable=false，`-strictplus` → variable=true（main.cpp:1623-1627）。

cleanStrict（alignment.cpp:2592-2739）：
1. `if(gInCol[i] > gapCut) saveResidues[i] = -1;`（:2601-2603）。
2. `if(MDK_W[i] < simCut) saveResidues[i] = -1;`（:2606-2608）。
3. deque 恢复孤立删除列（3/4 邻居规则）：载入前 `min(5, residNumber)` 个值（:2614-2615）；特殊情形（列 0；列 1 两种宽度；:2619-2638）；主循环 `for(i = 2; i < residNumber - 2; i++)`：deque 布局 `[i-2, i-1, i, i+1, i+2]`，检查槽 0,1,3,4 四邻居 `num >= 3` 则 `saveResidues[i] = i`，随后 `if((i+3) < residNumber) { pop_front(); push_back(saveResidues[i+3]); }`（:2641-2659）；末两列特殊情形用 deque 尾部值（`size()-4, size()-3, size()-1` 与 `size()-3, size()-2`，:2665-2678）。
4. 块大小：`!variable → lenBlock = 5`；`variable → lenBlock = utils::roundInt(residNumber * 0.01)` 且夹在 [3,12]：`lenBlock = lenBlock > 3 ? (lenBlock < 12 ? lenBlock : 12) : 3;`（:2682-2687）。
5. `blockSize = blockSize > 0 ? blockSize : lenBlock;`（用户 -block 优先；:2690）→ `removeSmallerBlocks(blockSize)`。
6. terminal / complementary / removeCols_SeqsAllGaps / `fillNewDataStructure(&counter)`（注意 strict 走 newValues 版本，分配在 removeCols_SeqsAllGaps 内，:2720）→ 构造。

### 2.8 -automated1（autAlignment.cpp:150-204 + main.cpp:1629-1634）
- **没有组合枚举，也没有一致性比较**：只是启发式选路。
- `selectMethod()`：确保 identities；对每行 i 求 `mx = max_{j≠i} identities[i][j]`、`avg = Σ_{j≠i} identities[i][j]`；`avgSeq += avg/(sequenNumber - 1); maxSeq += mx;` 循环后 `avgSeq /= sequenNumber; maxSeq /= sequenNumber;`（:180-181）。
- 判定：`avgSeq >= 0.55 → GAPPYOUT`；`avgSeq <= 0.38 → STRICT`；否则 `sequenNumber <= 20 → GAPPYOUT`；否则 `(maxSeq >= 0.5 && maxSeq <= 0.65) → GAPPYOUT`，else STRICT。
- main 中：GAPPYOUT → `clean2ndSlope(complementary)`；STRICT → `cleanCombMethods(complementary, false)`（main.cpp:1630-1633）。

### 2.9 -ct/-compareset/-forceselect（cleanCompareFile，alignment.cpp:851-890）

```cpp
utils::copyVect(vectValues, vectAux, residNumber);
utils::quicksort(vectAux, 0, residNumber-1);
cut = vectAux[(int)((float)(residNumber - 1) * (100.0 - baseLine)/100.0)];
cut = cutpoint < cut ? cutpoint : cut;    // 取较小者(更宽松)
cleanByCutValue(cut, baseLine, vectValues, complementary);   // float 版: 保留 > cut
```

main 中 `-ct` 时 `conserve` 作 baseLine（main.cpp:1638-1639）。

### 2.10 -resoverlap/-seqoverlap（cleanSpuriousSeq，alignment.cpp:906-938 + 1132-1234 + 2743-2788）
1. `floatOverlap = overlap * float(sequenNumber-1); ovrlap = int(overlap * (sequenNumber-1)); if(floatOverlap > float(ovrlap)) ovrlap++;`（即 ceil，:1140-1144）。
2. 每序列 i、每列 j：对 k≠i 统计 hit：`seq[i][j]==seq[k][j]`（**含同为 '-' 或同为 indet 的情形，也算 hit**）或（两者都有效且不同也 hit，即"双方均有效"就 hit）（:1179-1209）。
3. `if(hit >= ovrlap) seqValue++;`；`spuriousVector[i] = ((float) seqValue / residNumber);`（:1217-1226）。
4. `cleanOverlapSeq(minimumOverlap, overlapVector, complementary)`：`overlapSeq[i] < minimumOverlap → saveSequences[i] = -1`（:2753-2755）；随后（complementary 只翻序列，`computeComplementaryAlig(false, true)`）；**main 中之后还会串一个 `cleanNoAllGaps(false)`**（main.cpp:1643-1648）；`seqOverlap/100` 转换在 main.cpp:1644。

### 2.11 -terminalonly / --set_boundaries
- 解析：`utils::readNumbers_StartEnd`（utils.cpp:640-659，只允许一个逗号两个数）；校验 `boundaries[1] < numAminos`（main.cpp:1230-1235）。
- `trimTerminalGaps(terminal, boundaries)`（alignment.cpp:1831-1841）取 `[0]`、`[1]`。各 cleanByCutValue/cleanStrict 在生成选择后调用 `removeOnlyTerminal()`（§2.0）。**若 -terminalonly 单独使用且无方法 → main.cpp:1238-1244 报错**。

### 2.12 -selectcols / -selectseqs
- 解析 `utils::readNumbers`（utils.cpp:662-719）：`numbers[0]=nElems`，成对展开 `a-b` 或单点 `a,a`；校验 `numbers[i-2] >= 0 && numbers[i-1] >= numbers[i-2]` 否则 NULL。main 校验上界（main.cpp:1674-1709）。
- `removeColumns(columns, 1, delColumns[0], complementary)`（alignment.cpp:2792-2836）：`for(i = init; i < size + init; i += 2) for(j = columns[i]; j <= columns[i+1]; j++) saveResidues[j] = -1;`。
- `removeSequences` 同构（:2841-2885）；随后 main 里对 selectseqs 结果再跑 `cleanNoAllGaps(false)`（main.cpp:1703-1707）。
- `-complementary` 与 `-colnumbering`（printCorrespondence，alignment.cpp:1802-1811，格式 `#ColumnsMap\tv0, v1, ...`）。

### 2.13 -clusters / -maxidentity（alignment.cpp:1249-1476）
- `calculateRepresentativeSeq(maximumIdent)`：按"去 gap 长度"对 `seqs[i] = {len, idx}` 快排（int** 版）；`cluster[0] = seqs[sequenNumber-1][1]`（最长者）；从次长往下，若与某簇代表 `identities[seqs[i][1]][cluster[j]] > maximumIdent` 取最大者入簇，否则 `cluster[clusterNum++] = seqs[i][1]`。返回 `repres = [clusterNum, cluster...]`（**static 局部指针**，:1253）。
- `getCutPointClusters(clusterNumber)`：clusterNumber==n → 返回 1；==1 → 返回 0；否则以平均 identity 为起点二分（上下界 gMax/gMin 更新为 startingPoint 再取均值；`iter > 10` 防死循环；:1410-1466）。簇判定条件是 `> startingPoint`（:1417）。
- `getClustering`（:1482-1542）生成选序列的新比对；main 里随后 `cleanNoAllGaps(false)`（main.cpp:1714-1729）。

## 3. similarityMatrix.cpp

### 3.1 默认矩阵（values.h，`#ifdef SIMMatrix` 守卫，仅被 similarityMatrix.cpp 包含）
- `listNTSym[6]="ACGTU"`（values.h:32）；`listAASym[21]="ARNDCQEGHILKMFPSTWYV"`（values.h:34）；`listNTDegenerateSym[16]="ACGTURYKMSWBDHV"`（values.h:36）。
- **AA 默认 = BLOSUM62**（`defaultAAMatrix[20][20]`，values.h:76-97），20×20 整数值（4/-1/-2/…/11），顺序按 `ARNDCQEGHILKMFPSTWYV`。存为 float。
- 核酸默认 = 5×5 单位阵（values.h:46-52）。
- **degenerated_nt_identity**：`defaultNTDegeneratedMatrix[15][15]`（values.h:54-73）——A/C/G/T/U 对角 1.0；两符号并集交集比：R/Y/K/M/S/W 行用 `1/4.`；B/D/H/V 行用 `1/6.`（浮点字面量，如 `1/4.` = 0.25f 精确，`1/6.` = float(0.16666667…)）。精确值需按行转录（values.h:63-72 给出完整 15 行）。
- 备选阵 `alternative_1_NTDegeneratedMatrix[15][15]` = 全单位阵（values.h:103-119），由 `--alternative_matrix degenerated_nt_identity`（main.cpp:218-226 → `alternativeSimilarityMatrices(1, datatype)`，similarityMatrix.cpp:351-414）。

### 3.2 矩阵装载/初始化
- `loadSimMatrix`（similarityMatrix.cpp:142-256）：首行去空格得符号表（必须 A-Z，`toupper`）；支持带/不带首列两种格式；**非对称阵自动对称化**：`if(simMat[i][j] != simMat[j][i]) { float value = (simMat[i][j] + simMat[j][i]) / 2.0; simMat[i][j]=value; simMat[j][i]=value; }`（:229-237，float+float 除 double→double→float）。
- **distMat = 基于相似度列向量的欧氏距离**（四处相同实现，:242-252 / :277-287 / :308-318 / :338-348）：

```cpp
for(j = 0; j < numPositions; j++)
  for(i = 0; i < numPositions; i++)
    if((i != j) && (distMat[i][j] == 0.0)) {
      for(k = 0, sum = 0; k < numPositions; k++)
        sum += ((simMat[k][j] - simMat[k][i]) * (simMat[k][j] - simMat[k][i]));
      sum = (float) sqrt(sum);          // sqrt(double) 后 cast float
      distMat[i][j] = sum; distMat[j][i] = sum;
    }
```

**求值顺序固定：外层 j、内层 i、k 从 0 递增**——k 的累加顺序影响 float 和的位级结果；`simMat[k][j] - simMat[k][i]` 顺序固定（j 列减 i 列）。例如核酸单位阵：不同核苷酸距离 = sqrt(2.0)→float(1.4142135f)，相同 = 0。
- `defaultAASimMatrix / defaultNTSimMatrix / defaultNTDegeneratedSimMatrix`（:258-349）只是免对称化直接装入再算距离。vhash 映射 `listXSym[i]-'A'`；`getDistance(a,b)` 大写化后查 `distMat[vhash[..]][vhash[..]]`，查不到打印错误并返回 -1（:438-465）。
- main.cpp:1526-1552 决定用哪个矩阵（strict/strictplus/automated1/st/-ssc/-sst 才建；`-matrix` 文件加载返回值**未检查**——main.cpp:1530）。
- `printMatrix()`：`cerr << setw(8) << setprecision(4) << right << simMat[i][j]`（:421-428）。

## 4. 滑动窗口 -w/-gw/-sw/-cw 与"median"

- 选项解析（main.cpp:868-999）：`-w/-gw/-sw/-cw` 均要求正整数（`atoi`，`<=0` 报错）；互斥规则：`-w` 与 gw/sw/cw 冲突；与自动方法、select 方法冲突。
- 应用（main.cpp:1500-1510）：`windowSize != -1 → gapWindow = simWindow = windowSize`；否则未设者 = 0。`setWindowsSize(gapWindow, simWindow)`。
- **均为"半窗口"**：window = 2*half+1。gap：statisticsGaps.applyWindow（§1.1，结果 `roundInt(mean)`，即本版所谓"median"实为四舍五入均值；**偶数/奇数窗口无取中值规则，无半值处理**）。similarity：statisticsConservation.applyWindow（纯平均 `sum/(float)window`）。consistency：compareFiles::applyWindow（compareFiles.cpp:379-429，平均，`columnsValue[i] /= window`；拒绝条件只有 `halfWindow > columns/4`）。
- 环形反射边界三处一致：`j<0 → vect[-j]`；`j>=n → vect[2n-j-2]`。
- compare 窗口应用时机（main.cpp:1385-1400）：`windowSize` 优先于 `conWindow`；forceComparison 分支同。
- `-cw` 在单比对模式下被静默忽略（除非同时 -matrix 触发报错，main.cpp:1202-1205）。
- 窗口超过 `columns/4` 时静默不生效（三个 applyWindow 均返回 false，无警告）。

## 5. 所有排序调用点（全部为自写 `utils::quicksort`，pivot=末元素，Lomuto 变体 utils.cpp:230-262(float)/287-318(int)/338-374(int**)）

| 调用点 | 对象 | 对输出的影响 |
|---|---|---|
| alignment.cpp:865 | cleanCompareFile：consistency 副本 | 选 percentile 值 → cut |
| alignment.cpp:1048 | cleanCombMethods：被选列 simil 副本 | first20/last80 点 → simCut |
| alignment.cpp:1272 | calculateRepresentativeSeq：{len,idx} 按 len | **不稳定！等长序列的相对顺序取决于交换过程 → 影响簇代表选择 → 影响输出序列集合**（Java 必须逐语句复刻排序含 swap） |
| alignment.cpp:1390 | getCutPointClusters：同上 | 同上 |
| alignment.cpp:2107 | cleanByCutValue(gap)：gInCol 副本（oth>0 时） | 重设 cut |
| alignment.cpp:2442-2443 | cleanByCutValue(组合)：MDK 与 gaps 副本 | blCons/blGaps |
| compareFiles.cpp:467 | printStatisticsFileAcl 副本 | 仅打印分布 |
| statisticsConservation.cpp:324 | calcCutPoint：MDK_Window 副本 | cut 值 |
| statisticsConservation.cpp:396 | printConservationAcl 副本 | 仅打印 |

quicksort 细节（utils.cpp:230-262）：`elem_div = vect[fin]`；`i=ini-1, j=fin`；`while(1){ while(vect[++i] < elem_div) if(i==fin) break; while(vect[--j] > elem_div) if(j==0) break; if(i<j) swap else break; }` 然后 `swap(vect[i], vect[fin])`，递归 `[ini, i-1]` 与 `[i+1, fin]`。入口守卫 `if((ini >= fin) || (fin < 0)) return;`。**对纯值数组结果与排序算法无关（多重集相同），但对 int**（带索引）结果依赖交换历史，必须逐位复刻。**

## 6. 浮点输出/比较点（Java 复刻格式）

- statisticsGaps.cpp:421-429（-sgc）：表头三行固定文本；`cout.precision(10)`；每行 `"  " << setw(5) << i << "\t\t" << setw(10) << (vectAux[i] * 100.0)/columnLength << "\t" << setw(7) << 1 -((vectAux[i] * 1.0)/columnLength)`。表达式：`int * 100.0`→double；`int * 1.0`→double；`1 - double`→double。**precision(10) 会污染后续所有 cout**。
- statisticsGaps.cpp:447-465（-sgt）：`cout.precision(10)`；`(numColumnsWithGaps[i] * 100.0)/columns`、`(acm * 100.0)/columns`、`(i * 1.0)/columnLength`、`1 - ((i * 1.0)/columnLength)`；首列带 `setiosflags(ios::left)`。
- statisticsConservation.cpp:360-374（-ssc/-sst 单列）：`cout.precision(10)`；`setw(7) << MDK_Window[i]`（float 值以默认格式、10 位有效数字打印）。
- statisticsConservation.cpp:399-427（-sst 累积）：`cout.precision(10)`；`((float) num/columns * 100.0)`（float/int→float 再 *double）、`setw(15) << refer`（float）。
- compareFiles.cpp:187-189（-compareset 进度输出，仅 verbosity）：`"Pond. Hits: " << setw(8) << value << "\t%Consistency: " << value/numResiduesAlig[i]`——**未设 precision（6 位有效数字）**；`205: "File Selected:\t" << fileNames[alig] << ... << "Value:\t\t" << max`。
- compareFiles.cpp:434-451（-sfc）：`cout.precision(10)`；`setw(5) << i + 1`（**从 1 开始计数**），值无 setw。
- compareFiles.cpp:457-520（-sft）：`cout.precision(10)`；`((float) num/numAminos * 100.0)` 等；`setw(15) << refer`。
- autAlignment.cpp:255-256 与 335-336（-sident/-soverlap）：`cout.precision(4); cout << fixed;`（4 位小数定点）；矩阵行 `setw(10) << identities[i][j] << "\t"`，对角线打印 `setw(10) << 1.00`；most-similar 表 `setw(5) << maxs[i][0]`。**fixed 状态同样污染后续输出**（之后若再有默认格式输出会受影响——注意 printSeqIdentity 后可能继续打印比对，Java 需复刻此状态泄漏或确认输出顺序）。
- autAlignment.cpp:442-445（-scolidentt）：`#maxColIdentity/#minColIdentity/#avgColIdentity/#stdColIdentity` 无精度设置（6 位有效）。
- similarityMatrix.cpp:421-428：`cerr << setw(8) << setprecision(4)`。
- 无任何 `printf` 格式串（全 iostream）；HTML 摘要所有数字为整数（rwAlignment.cpp:1958-1962 等）。
- main.cpp:91-95：版本行 `trimAl v1.5.rev1 build[2025-11-25]`（defines.h:34-37：BUILD "2025-11-25"、VERSION 1.5、REVISION 1——注意 VERSION 是 double 1.5，`"trimAl v" << VERSION` 打 "1.5"）。

## 7. 剪裁后同步机制

- **sequencesMatrix**（sequencesMatrix.cpp:41-72）：构造时 `matrix[i][j] = k++`（每行从 1 开始对非 gap 字符编号，gap 为 0）。剪裁时**不就地修改**；新 alignment 对象的 seqMatrix=NULL，懒构建（alignment.cpp:1968-1971）。用于：比对比较（getColumn（:134-160）`getColumn(value,row,..)` 是"在 row 行找等于 value 的列，返回该列全序列向量，找不到返回 -1 向量"——注意 `for(i=0;...) if(matrix[row][i]==value) break;` 若 value=0 会匹配第一个 gap！但调用方先判断 `columnSeqMatrix[k] != 0` 才调用）、`setOrder`（按名字序重排行，:162-179）、backtranslation 的 `getSequence`（:181-201）。
- **backtranslation**（main.cpp:1609-1610, 1750-1759；alignment.cpp:1551-1650）：`seqMatrix = origAlig->getSeqMatrix()`（**原始比对**的矩阵，在剪裁前抓取）；`getTranslationCDS(newResidues, newSequences, singleAlig->getCorrespResidues(), seqNames, seqMatrix, singleAlig)`：`selectedRes[]` = 全 i，然后被 `ColumnsToKeep` 之间的区间置 -1（:1576-1587，注意 `for(j = 0; j < ColumnsToKeep[0]; j++) selectedRes[j] = -1;`——把第一个保留列之前的全删，之后对相邻保留列之间置 -1，最后尾部置 -1）；对每个保留残基输出对应 3 个密码子（`sequences[mappedSeqs[i]][k]` 越界补 'N'，:1603-1611）；gap 列输出 `---`（:1612-1614）。新比对宽度 `newResidues*3`。
- **HTML 摘要**（rwAlignment.cpp:1861-2135，`alignmentSummaryHTML`）：main.cpp:1741-1746 用 `origAlig`（原比对）+ `singleAlig->getCorrespResidues()/getCorrespSequences()`（压缩后的映射数组，长度 = 新列/序列数）重建 `res[]/seq[]` 布尔表；`gapsValues = sgaps->getGapsWindow()`、`simValues = scons->getMdkwVector()`（**原比对的统计**，含窗口）；consValues = compareVect（一致性，可为 NULL）。色阶阈值：gap 用 `1 - (float)gapsValues[k]/sequenNumber` 与 .750/.500/.350/.250/.200/.150/.100/.050/.001 比较，==0 → c12、==sequenNumber → c1（:2042-2065）；similarity 用 MDK 值与 1/0/.750/.500/.250/.100/.010/.001/1e-4/1e-5/1e-6（:2070-2093）；consistency 与 1/0/.750/.500/.350/.250/.200/.150/.100/.050/.001（:2098-2121）。块宽 HTMLBLOCKS=120，`minHTML = max(25, maxLongName+10)`。残基着色（CLUSTAL 方案）在 `utils::determineColor`（utils.cpp:722-832）+ `lookForPattern`（utils.cpp:835-851，`count/column.size() >= level`）。
- **complementary 对齐**：`computeComplementaryAlig(residues, sequences)`（alignment.cpp:2889-2897）在最终 `removeCols_SeqsAllGaps` **之前**翻转选择，因此互补比对同样会剔除"选择范围内的全 gap 列/序列"（并打印相应 WARNING）。
- `removeCols_SeqsAllGaps()`（alignment.cpp:2980-3049）：两遍扫描（列先、序列后），只统计有效选择内的全 gap；每次删除打 `WARNING: Removing column 'i' composed only by gaps` / 序列同理；`keepSequences` 时保留序列但计数（:3029-3037）；返回 `counter` 并分配 `counter.matrix/seqsName`。`fillNewDataStructure`（:3053-3091）逐字符 `resize(size+1, ch)` 复制。

## 8. compareFiles.cpp（-compareset/-forceselect/-ct/-cs…）

**algorithm()**（compareFiles.cpp:40-238）——选"最一致"文件：
1. `numSeqs = vectAlignments[0]->getNumSpecies()`；校验其余文件序列数相同、且名字可与 [0] 建立映射（getSeqNameOrder，alignment.cpp:1937-1960）；然后**把每个比对的 seqMatrix 按第 0 个比对的序列顺序重排**（setSeqMatrixOrder，:92-96）。
2. 对每个比对 i：`vectHits[i][j]` = 列 j 的一致性。计算：取第 i 比对第 j 列的残基位置向量 `columnSeqMatrix[k]`（残基在该序列内的序号，gap=0）；对每个非 0 的 k，向其余**所有**比对 l（先 l<i 后 l>i）取 `getColumnSeqMatrix(columnSeqMatrix[k], k, aux)`；对 m>k 的非 0 残基：`columnSeqMatrix[m] == columnSeqMatrixAux[m]` → hits++；pairRes++（:134-170）。
3. `if(pairRes != 0) { vectHits[i][j] += ((1.0 * hits)/pairRes); value += vectHits[i][j]; }`（:176-179；`1.0*hits` → double 除法再赋 float）。
4. verbosity 输出（§6）。`if((value/numResiduesAlig[i]) > max) { alig = i; max = value/numResiduesAlig[i]; }`（:194-197）——**严格大于**，初值 max=0/alig=0：全 0 时选第 0 个文件；并列时选先出现者。
5. 打印 `"File Selected:\t..."`（:203-206）；`columnsValue` 初始化为 -1 再填 `vectHits[alig]`（:212-216，**长度按选中文件的列数**）。返回 alig 或 -1。
- main.cpp:1331-1410 流程：读文件列表（逐行，最多 256 字节行缓冲计数 numfiles），逐个 load+`sequenMatrix()`+数据类型一致性检查+跟踪 maxAminos；`referFile = algorithm(...)`（verbosity = `(stats >= 0) && (outfile != NULL)`，main.cpp:1380-1383）；然后 `origAlig->loadAlignment(filesToCompare[referFile])` 重新加载选中文件。

**forceComparison()**（compareFiles.cpp:244-374）：以 selected（-forceselect 载入的 origAlig）为参考；同样先检查+重排所有比对矩阵（:277-306，注意这里循环从 i=0 开始，包含 selected 自身？不——vectAlignments 是 compAlig 数组，不含 selected）；每列：`columnsValue[i] += ((1.0 * hit)/pairRes)`（:355），对每个非 0 残基 j 遍历全部 k 个比对（**包含与自身比较**，k ∈ [0,numAlignments)，:329-348——与 algorithm 不同，algorithm 排除了 i 自身！这是两函数的差异点）。

**applyWindow**（compareFiles.cpp:379-429）：`if(halfWindow > columns/4) return false;`（无"已算过"检查）；in place 平均。

## 9. 怪异行为 / bug / 魔法数字（必须复刻或知晓）

1. **statisticsConservation.cpp:326-329 越界读**：`for(i = columns-1; i >= 0; i--) if(vectAux[i] < similThreshold) break; cuttingPoint_SimilThreshold = vectAux[i];` —— 所有列 >= 阈值时 i=-1，`vectAux[-1]`。C++ 实际值取决于堆内容；Java 复刻需决策（常见实际后果：读到栈/堆邻近 float，通常表现为极大值 → cut 取 MinumumConserv 分支。建议复刻"语义近似"：如 simThreshold<=所有值时用 vectAux[0]?—— 无法确定，需对照二进制/实测输出校准）。
2. **cleanStrict 末两列特殊情形**（alignment.cpp:2665-2678）：residNumber<2 时 `saveResidues[residNumber-2]` 负下标写（极短比对崩溃）；正常长度时逻辑如上转录。
3. **cleanByCutValue(conservation) 放宽只回收 `== cut` 的列**（:2286,2309）而 gap 版回收 `<= cut` 且会重算 cut——两版不对称是有意行为。
4. **calcCutPoint2ndSlope 未初始化与越界**：`secondSlopeVector = new float[maxGaps+1]` 只初始化前 maxGaps 项（statisticsGaps.cpp:338-339）；内层 `while(numColumnsWithGaps[act]==0) act++` 无上界（:352-353 等），理论上可越过 numColumnsWithGaps（大小 columnLength+1）末元素——当 maxGaps==columnLength 且尾桶为 0 时越界读。复刻时建议按"循环至 maxIter 封顶"的可见行为处理并实测校准。
5. **calculateMatrixIdentity length==0 → NaN**（statisticsConservation.cpp:159）；calculateSeqIdentity dst==0 → inf（autAlignment.cpp:71）；calculateSeqOverlap referenceLength==0 → inf（autAlignment.cpp:106）。NaN/inf 会沿 `num/den` 传播——Java 用 Float.NaN/Infinity 与 C float 语义一致（IEEE754），比较行为（NaN 参与的 `>`/`<` 恒 false）也一致。
6. **static 局部指针**：`utils::readNumbers/readNumbers_StartEnd`（utils.cpp:643,666）与 `calculateRepresentativeSeq` 的 `repres`（alignment.cpp:1253）——多次调用共享缓冲；单次运行通常无影响，但要知道。
7. **`operator=` 与拷贝构造的双重 bug**：alignment.cpp:334-339 `delete [] residuesNumber; if(old.residuesNumber){...} residuesNumber = NULL;`（拷贝后立即置 NULL——分配被丢弃）；:390-397 `delete seqMatrix; seqMatrix = old.seqMatrix;`（**浅共享指针**，双方析构会 double-free——单线程主流程中因旧对象先被 delete 的顺序而侥幸存活）；构造函数 :230-273 先建 identities/overlaps 再立刻置 NULL（:272）——**拷贝构造里 identities 被构建后又被覆写为 NULL（内存泄漏 + 丢失）**：新对象 identities 总是 NULL，之后用到会重算。复刻时按"identities=NULL"的最终行为处理。
8. **compareFiles::forceComparison 包含与自身比较**（:329 `for(k = 0; k < numAlignments; k++)`）而 algorithm 排除自身（:141,156 两个区间绕开 i）——两函数 hit/pairRes 基数不同。
9. **main.cpp:1530** `similMatrix->loadSimMatrix(matrix)` 返回值未检查——坏矩阵文件会以半初始化状态继续（memoryAllocation 已跑，vhash/simMat 部分填充）。
10. **main.cpp:1734-1746 NULL 解引用路径**：若唯一指定的修剪方法返回 NULL，`singleAlig = origAlig; origAlig = NULL;` 后 `origAlig->alignmentSummaryHTML(...)` 对 NULL 调用（仅当 -htmlout 与失败方法同时出现）。
11. **applyWindow 拒绝条件用 `columns/4`（列数而非序列数）**——对 gap 和 conservation 一致；-w 过大时静默失效。
12. **precision 状态泄漏**：`cout.precision(10)`（statisticsGaps/Conservation/compareFiles 打印函数）与 `cout << fixed; cout.precision(4)`（autAlignment）设置后不恢复；后续所有 cout（包括 `#ColumnsMap`、序列输出?——序列输出走 `file`（ofstream/cout 引用），printAlignment 用 cout 时也会受影响——但 printAlignment 在 stats 打印之后执行，main.cpp:1769-1770，**如果同时给了 -sgc 之类的统计与 stdout 输出，比对输出的浮点不受影响（无浮点），但若有任何浮点打印将延续 precision(10)**）。
13. **-gt 的 `1 - atof(x)`**（main.cpp:307）float 存储后再参与 `(double)columnLength * gapThreshold`——0.1 之类值有表示误差，cut 的 `.5` 进位边界可能与"精确十进制"不同。
14. 魔法数字集合：0.8（MDK 清零 gap 比例，statisticsConservation.cpp:206）、3/4 邻居（cleanStrict）、5（strict 块）、[3,12]（strictplus 块，1% 长度）、0.005（放宽块初始 0.5%）、20/80（strict 分位点）、/10 与 pow(10)（对数插值）、0.55/0.38/0.5/0.65/20（selectMethod 阈值）、10 次迭代上限（getCutPointClusters）、25/maxLongName+10（HTML 空白）、120（HTML 块）、60（FASTA/Clustal 行宽）、50/10（Phylip/MEGA/NEXUS 块）、10（PHYLIPDISTANCE 名字宽）、100（checkTypeAlignment 采样上限常量已失效——utils.cpp:589 被注释，现扫描全序列）、0.95（核酸比例阈值，utils.cpp:611-614）。
15. **-sgc/-sgt/-ssc/-sst 的 stats 计数怪癖**（main.cpp:1044-1113, 1556-1596）：每 flag `stats--`；`if(stats<0) stats--;`（再次递减）；打印时 `stats++` 且 `if(stats < -1) cout << endl;`——控制多统计输出之间的空行。
16. **statAl `-scolidentt`** 是唯一使用 `calculateColIdentity` 的入口；trimAl 主程序没有 -colscore。
17. **printCorrespondence**（alignment.cpp:1802-1811）：`for(i = 0; i < residNumber - 1; i++) cout << saveResidues[i] << ", "; cout << saveResidues[i] << endl;`——循环结束后 i==residNumber-1，恰好打印最后一个；**residNumber==0 时 i 未初始化即用**（理论 UB）。
18. `alignment::printAlignment`/`saveAlignment` 空比对仅告警并返回 true（alignment.cpp:563-567, 619-623）。
19. `sequencesMatrix::operator=` 自拷贝 bug（sequencesMatrix.cpp:90-95 `matrix[i][j] = matrix[i][j]`，复制的是自身未初始化内存）——该运算符当前无调用点。
20. **`numColumnsWithGaps` 大小是 `species+1`（序列数），而直方图索引上限是 maxGaps（≤species）**——calcCutPoint 循环 `i < columnLength`（species）会遍历到最后一个合法桶。

## 10. -w 系列在三种数据流中的精确接线（补充表）

| 选项 | 解析 | 作用于 |
|---|---|---|
| -w | windowSize | gapWindow=simWindow=windowSize（main.cpp:1500-1503）；compare 窗口（1386,1398） |
| -gw | gapWindow | 仅 statisticsGaps.applyWindow（alignment.cpp:1668） |
| -sw | simWindow | 仅 statisticsConservation.applyWindow（alignment.cpp:1766） |
| -cw | conWindow | 仅 compareFiles::applyWindow（main.cpp:1388,1400）；单比对时无效 |
| 默认 | 未设 → 0 | gap: 原始计数；cons: MDK_Window=MDK（window=1） |

## 11. main 的方法分派顺序（互斥并非全互斥，按代码顺序覆盖，main.cpp:1614-1730）

`nogaps → noallgaps → gappyout → strict → strictplus → automated1 → (comThreshold) → (resOverlap&seqOverlap) → (simThreshold[+gapThreshold]) → else gapThreshold → (selectCols/selectSeqs) → (maxIdentity / clusters)`。任何方法都未触发时 `singleAlig = origAlig; origAlig = NULL`（:1734-1737）。`-ct` 与 `-in` 互斥（-ct 必须配 -compareset/-forceselect，main.cpp:352-383）；`-sfc/-sft` 亦与 -in 互斥（:1087,1103）。

---

**校准建议**：Java 复刻的位级风险点按优先级为：(a) calculateVectors 的浮点累加顺序与 float/double 混合舍入（§1.2）；(b) distMat 欧氏距离的 k 循环顺序与 (float)sqrt 截断（§3.2）；(c) roundInt 的 `(int)(x+0.5)`（负数行为不同，本库输入均非负）；(d) strict 的 20/80 分点与 log10/pow 链（§2.7）；(e) int** 快排的不稳定性对 clusters/-maxidentity 的影响（§5）；(f) statisticsConservation::calcCutPoint 的 i==-1 越界（需实测对照 upstream/source/trimal 二进制输出校准，upstream/source/ 内已有编译好的 trimal、statal、readal 可用于差分测试）。
