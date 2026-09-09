# Semiconductor Wafer Process Variation and Quality Analysis

## 半導體晶圓製程變異與品質分析

本專案以 **R** 進行半導體晶圓製程資料之統計品質管制分析，使用 245 片 wafer 的量測資料與 2,740 個二元製程變數，探討整體製程分布、LOT 間差異、統計管制圖、高相關製程變數、變數篩選，以及平均與變異數模型。

本 repository 為原始統計品質管制課程專案的重新整理版本，將原先分散於 R Markdown、課堂程式與分析紀錄中的內容整合為單一可執行的 R script，保留課程專案所採用的主要分析方法，並移除重複、測試性與不必要的程式碼，以提升可讀性與可重現性。

---

## 專案文件

- [完整專案報告](report/半導體晶圓製程變異與品質分析_專案報告.pdf)
- [R 分析程式](wafer_quality_analysis.R)
- [資料集](data/new_wafer_data.csv)

---

## 專案概述

資料包含：

- `LOT`：批次編號，共 10 個 LOT
- `WAFER`：各 LOT 內的 wafer 編號
- `OBSERVATION`：品質量測值
- `T0001` ~ `T2740`：2,740 個二元製程變數

資料規模：

| 項目 | 數量 |
|---|---:|
| Wafer samples | 245 |
| LOTs | 10 |
| Process variables | 2,740 |
| Constant process variables | 15 |
| Missing values | 0 |
| Duplicate rows | 0 |

---

## 分析目標

本專案主要探討以下問題：

1. `OBSERVATION` 的整體分布與 LOT 間差異為何？
2. 製程是否存在明顯的批次間與晶圓間變異？
3. X-bar、R、S 與 Individuals Chart 顯示哪些製程異常？
4. 2,740 個高維製程變數中，是否存在高度相關的「孿生變數」？
5. OGA 與 LASSO 會篩選出哪些重要製程變數？
6. 是否可以同時找出影響製程平均水準與變異程度的製程變數？

---

## 分析流程

```text
Data Loading & Cleaning
        ↓
Exploratory Data Analysis
        ↓
Process Capability
        ↓
Control Charts
        ↓
Twin-variable Correlation Analysis
        ↓
OGA Variable Selection
        ↓
LASSO Variable Selection
        ↓
Two-stage OGA + HDIC + Trim
```

---

## 分析方法

### 1. 探索性資料分析

針對 `OBSERVATION` 進行整體分布與 LOT 間差異分析，包含：

- Summary statistics
- Histogram
- Normal Q-Q plot
- Shapiro-Wilk normality test
- LOT-level boxplots
- Wafer sequence visualization
- Kruskal-Wallis test

Shapiro-Wilk 常態性檢定結果：

```text
W = 0.96979
p-value = 4.574e-05
```

結果顯示 `OBSERVATION` 並不完全符合常態分布。

不同 LOT 的分布位置亦存在顯著差異：

```text
Kruskal-Wallis chi-squared = 178.92
df = 9
p-value < 2.2e-16
```

---

### 2. 製程能力分析

專案保留課堂分析中使用的不同規格定義方式。

#### 課堂假設規格

```text
LSL = 12
USL = 20
```

分析結果：

| 指標 | 數值 |
|---|---:|
| Mean | 16.2925 |
| SD | 1.8345 |
| Cp | 0.7268 |
| Cpk | 0.6737 |
| Outside assumed limits | 2.8571% |

#### 1% / 99% 分位數規格代理

以資料的 1% 與 99% 分位數作為規格代理：

```text
Cp ≈ 0.741
```

#### Mean ± 3 SD 自然界限

以：

```text
Mean ± 3 SD
```

作為自然界限時，樣本超出範圍的比例為：

```text
0.4082%
```

> 分位數界限與 Mean ± 3 SD 並非正式工程規格，因此此部分主要作為課堂中製程能力與資料分布分析的實作。

---

### 3. 統計製程管制

建立以下 Statistical Process Control（SPC）圖表：

- Individuals Chart
- X-bar Chart
- R Chart
- S Chart
- 各 LOT 每 5 片 wafer 的平均變化

Global Individuals Chart 的主要參數：

```text
Center = 16.29249
StdDev = 0.9653863
LCL = 13.39633
UCL = 19.18865
```

X-bar Chart 顯示多個 LOT 平均值超出共同管制界限，反映資料中存在明顯的 LOT-to-LOT variation。

---

### 4. 孿生製程變數分析

由於資料包含大量二元製程變數，因此先移除 15 個 constant variables，再計算 Pearson correlation。

課堂專案將：

```text
0.8 < |r| < 1
```

定義為高度相關的孿生製程變數。

分析結果：

```text
Highly correlated records = 35,026
Connected clusters = 38
Independent variables = 253
Cluster representatives = 38
Variables after cluster reduction = 291
```

每個 correlation cluster 中，選擇與 `OBSERVATION` 絕對相關程度最高的變數作為該群集代表。

---

### 5. Orthogonal Greedy Algorithm（OGA）

OGA 依序選擇與目前 residual 關聯程度最大的變數，用於高維製程資料的變數篩選。

#### Raw OGA

直接使用全部 2,740 個製程變數：

```text
Selected variables = 27
```

前幾個被選出的變數：

```text
T0404, T1370, T1241, T2395, T1202, T0447, ...
```

#### Twin-aware OGA

先進行孿生變數群集整理，再對 291 個代表變數與獨立變數執行 OGA：

```text
Selected variables = 32
```

前幾個被選出的變數：

```text
T0404, T0050, T2045, T0450, T0447, T0872, ...
```

---

### 6. LASSO

使用 `glmnet` 建立 LASSO regression，並利用 10-fold cross-validation 選擇懲罰參數 λ。

```r
set.seed(42)

cv.glmnet(
  x = X,
  y = y,
  alpha = 1,
  nfolds = 10
)
```

分析結果：

| Criterion | Lambda | Selected Variables |
|---|---:|---:|
| `lambda.min` | 0.05726786 | 54 |
| `lambda.1se` | 0.2311913 | 18 |

在兩種 λ 設定中皆被選取的變數包含：

```text
T0404
T0447
T0510
T0851
T1026
T1203
T1490
T2352
T2398
```

---

### 7. Two-stage OGA + HDIC + Trim

最後使用 Two-stage OGA + HDIC + Trim，同時進行：

- Mean model variable selection
- Dispersion model variable selection

#### Mean model

選取變數：

```text
T0404
T1370
```

模型結果：

```text
E(OBSERVATION)
= 14.992340
+ 2.811701 × T0404
+ 1.698212 × T1370
```

#### Dispersion model

選取變數：

```text
T0456
```

模型結果：

```text
log(σ²)
= -0.2905104
+ 0.5958349 × T0456
```

此方法提供同時分析製程平均水準與變異程度的方式，可用於辨識可能影響製程中心位置與穩定性的製程變數。

---

## 主要結果

| 分析項目 | 主要結果 |
|---|---|
| Shapiro-Wilk | `p = 4.574e-05` |
| LOT comparison | Kruskal-Wallis `p < 2.2e-16` |
| Quantile-based Cp | `0.741` |
| Mean ± 3 SD outside rate | `0.4082%` |
| Twin-variable correlation records | `35,026` |
| Twin-variable clusters | `38` |
| Variables after cluster reduction | `291` |
| Raw OGA | `27` variables |
| Twin-aware OGA | `32` variables |
| LASSO `lambda.min` | `54` variables |
| LASSO `lambda.1se` | `18` variables |
| Two-stage mean variables | `T0404`, `T1370` |
| Two-stage dispersion variable | `T0456` |

---

## 圖表

R script 會產生以下主要圖表：

- OBSERVATION histogram
- Normal Q-Q plot
- LOT-level boxplots
- LOT × wafer sequence plots
- Individuals Chart
- LOT grouped mean plot
- X-bar Chart
- R Chart
- S Chart
- LASSO cross-validation curve
- Two-stage mean / dispersion model plot

完整圖表與分析說明可參考：

[半導體晶圓製程變異與品質分析－專案報告](report/半導體晶圓製程變異與品質分析_專案報告.pdf)

---

## 專案結構

```text
semiconductor-wafer-quality-analysis/
│
├── README.md
├── wafer_quality_analysis.R
│
├── data/
│   └── new_wafer_data.csv
│
└── report/
    └── 半導體晶圓製程變異與品質分析_專案報告.pdf
```

---

## 執行方式

### 1. 安裝所需套件

```r
install.packages(c(
  "dplyr",
  "ggplot2",
  "qcc",
  "igraph",
  "glmnet"
))
```

### 2. 資料路徑

資料檔預設位於：

```text
data/new_wafer_data.csv
```

R script 使用相對路徑讀取資料：

```r
file_path <- file.path(
  "data",
  "new_wafer_data.csv"
)
```

因此執行程式時，working directory 應位於 repository 根目錄。

### 3. 執行分析

在 RStudio 中開啟：

```text
wafer_quality_analysis.R
```

由專案根目錄依序執行程式，即可重現主要分析結果與圖表。

---

## 使用工具

- R
- RStudio
- `dplyr`
- `ggplot2`
- `qcc`
- `igraph`
- `glmnet`

---

## 專案背景與性質

本專案原為統計品質管制課程之團隊專案，主要針對半導體晶圓製程資料進行統計分析與品質管制方法實作。

原始分析為團隊合作成果。

本 repository 為後續重新整理與整合後的版本，將原先分散於 R Markdown、課堂程式與分析紀錄中的內容重新整合為單一 R script，並重新整理專案文件與分析報告，以提升程式碼可讀性、分析流程一致性與結果可重現性。

---

## 分析限制

- `T0001`–`T2740` 為匿名二元製程變數，缺乏實際設備與製程 metadata，因此結果主要解讀為統計關聯。
- 部分製程能力分析使用假設規格或資料分位數作為 proxy，並非正式工程 specification。
- 高相關變數代表統計上的相似行為，不應直接視為實際相同設備或製程步驟。
- 本專案重點為統計品質管制方法的實作與比較，而非建立正式量產製程監控系統。
