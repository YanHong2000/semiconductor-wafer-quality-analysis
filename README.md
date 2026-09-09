# Semiconductor Wafer Process Variation and Quality Analysis

半導體晶圓製程變異與品質分析

這是一個以 **R** 完成的統計品質管制課程專案，使用 245 片 wafer 的量測資料與 2,740 個二元製程變數，分析製程分布、LOT 間差異、管制圖、高相關製程變數、變數篩選，以及平均與變異數模型。

本 repository 將原始課堂程式重新整理成單一可執行的 R script，保留課程中使用並與老師討論過的主要分析方法，同時移除重複、測試性與不必要的程式碼。

---

## Project Overview

資料中包含：

- `LOT`：批次編號，共 10 個 LOT
- `WAFER`：各 LOT 內的 wafer 編號
- `OBSERVATION`：品質量測值
- `T0001` ~ `T2740`：2,740 個二元製程變數

資料規模：

| Item | Value |
|---|---:|
| Wafer samples | 245 |
| LOTs | 10 |
| Process variables | 2,740 |
| Constant process variables | 15 |
| Missing values | 0 |
| Duplicate rows | 0 |

---

## Objectives

本專案主要回答以下問題：

1. `OBSERVATION` 的整體分布與 LOT 間差異如何？
2. 製程是否存在明顯的批次間與晶圓間變異？
3. X-bar、R、S 與 Individuals Chart 顯示哪些製程異常？
4. 2,740 個高維製程變數中，是否存在高度相關的「孿生變數」？
5. OGA 與 LASSO 會選出哪些重要製程變數？
6. 是否可以同時找出影響平均水準與變異程度的製程變數？

---

## Analysis Workflow

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

## Methods

### 1. Exploratory Data Analysis

分析 `OBSERVATION` 的整體分布，包含：

- Summary statistics
- Histogram
- Normal Q-Q plot
- Shapiro-Wilk normality test
- LOT-level boxplots
- Wafer sequence visualization
- Kruskal-Wallis test

Shapiro-Wilk 檢定結果：

```text
W = 0.96979
p-value = 4.574e-05
```

不同 LOT 的分布位置也存在顯著差異：

```text
Kruskal-Wallis chi-squared = 178.92
df = 9
p-value < 2.2e-16
```

---

### 2. Process Capability

專案保留課堂中使用的不同規格定義。

#### Assumed engineering limits

```text
LSL = 12
USL = 20
```

結果：

| Metric | Value |
|---|---:|
| Mean | 16.2925 |
| SD | 1.8345 |
| Cp | 0.7268 |
| Cpk | 0.6737 |
| Outside assumed limits | 2.8571% |

#### 1% / 99% quantile proxy

以資料的 1% 與 99% 分位數作為規格代理：

```text
Cp ≈ 0.741
```

#### Natural limits

以：

```text
Mean ± 3 SD
```

作為自然界限時，樣本超出比例為：

```text
0.4082%
```

> Quantile limits 與 natural limits 並非正式工程規格，因此相關結果主要作為課堂中的製程能力與資料分布練習。

---

### 3. Statistical Process Control

建立以下管制圖：

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

X-bar Chart 顯示多個 LOT 平均值超出共同管制界限，反映明顯的 LOT-to-LOT variation。

---

### 4. Twin-variable Correlation Analysis

製程變數為大量二元欄位，因此先移除 15 個 constant variables，再計算 Pearson correlation。

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

每個 correlation cluster 選擇與 `OBSERVATION` 絕對相關最高的變數作為代表。

---

### 5. Orthogonal Greedy Algorithm (OGA)

OGA 依序選擇與目前 residual 關聯最大的變數，逐步建立高維變數篩選結果。

#### Raw OGA

直接使用全部 2,740 個製程變數：

```text
Selected variables = 27
```

前幾個變數：

```text
T0404, T1370, T1241, T2395, T1202, T0447, ...
```

#### Twin-aware OGA

先進行孿生變數群集整理，再對 291 個代表 / 獨立變數執行 OGA：

```text
Selected variables = 32
```

前幾個變數：

```text
T0404, T0050, T2045, T0450, T0447, T0872, ...
```

---

### 6. LASSO

使用 `glmnet` 建立 LASSO regression，並使用 10-fold cross-validation 選擇懲罰參數。

```r
set.seed(42)

cv.glmnet(
  x = X,
  y = y,
  alpha = 1,
  nfolds = 10
)
```

結果：

| Criterion | Lambda | Selected Variables |
|---|---:|---:|
| `lambda.min` | 0.05726786 | 54 |
| `lambda.1se` | 0.2311913 | 18 |

在兩種 lambda 設定中反覆出現的變數包含：

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

最後使用 Two-stage OGA + HDIC + Trim，同時選擇：

- Mean model variables
- Dispersion model variables

#### Mean model

選出：

```text
T0404
T1370
```

模型：

```text
E(OBSERVATION)
= 14.992340
+ 2.811701 × T0404
+ 1.698212 × T1370
```

#### Dispersion model

選出：

```text
T0456
```

模型：

```text
log(σ²)
= -0.2905104
+ 0.5958349 × T0456
```

此模型提供一個同時檢查製程平均水準與變異程度的分析方式。

---

## Key Results

| Analysis | Main Result |
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

## Visualizations

The R script generates:

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

---

## Repository Structure

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

> 若資料檔不適合公開，可將 `data/` 加入 `.gitignore`，並在 README 中保留資料欄位與格式說明即可。

---

## How to Run

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

### 2. Modify the data path

在 `wafer_quality_analysis.R` 中修改：

```r
file_path <- file.path(
  "data",
  "new_wafer_data.csv"
)
```

為自己的資料路徑。

### 3. Run the analysis

在 RStudio 中開啟：

```text
wafer_quality_analysis.R
```

依序執行即可重現主要分析結果與圖表。

---

## Tools

- R
- RStudio
- `dplyr`
- `ggplot2`
- `qcc`
- `igraph`
- `glmnet`

---

## Project Context

本專案原為統計品質管制課程之團隊專案，主要針對半導體晶圓製程資料進行統計分析與品質管制方法實作。

本 repository 為原始課堂專案的整理版本，將原先分散於 R Markdown、課堂程式與分析紀錄中的內容重新整合為單一 R script，包含資料前處理、探索性資料分析、製程能力分析、管制圖、OGA、LASSO，以及 Two-stage OGA + HDIC + Trim 等方法。

整理後的版本主要著重於提升程式碼的可讀性、分析流程的一致性與結果的可重現性，方便後續閱讀、維護與作品集展示。

本專案原始分析為團隊合作成果，本 repository 則為後續重新整理與整合後的版本。

---

## Limitations

- `T0001`–`T2740` 為匿名二元製程變數，缺乏實際設備 / 製程 metadata，因此結果主要解讀為統計關聯。
- 部分製程能力分析使用假設規格或資料分位數作為 proxy，並非正式工程 specification。
- 高相關變數代表的是統計上的相似行為，不應直接視為實際相同設備。
- 本專案的重點為統計品質管制課程方法的實作與比較，而非建立正式量產製程監控系統。

---

## Author

**Lin Yan Hong / 林彥宏**
