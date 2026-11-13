## Normalization and Pseudo-control–based Batch Effect Correction

RNA-seq and ATAC-seq datasets were normalized and corrected for technical variation and production-center effects using the workflow below. In general, data generated within the same production center were initially processed together.

### Step 1. Within-center normalization and removal of unwanted variation

RNA-seq gene counts were first normalized using the relative log expression (RLE) method. Residual variation across exposure groups was then estimated using `RUVr` from the `RUVSeq` package.

Three unwanted factors (`k = 3`) were incorporated into a general linear modeling framework to account for technical variation unrelated to exposure, including variation associated with tissue dissection and library preparation.

The same normalization framework was applied to ATAC-seq peak-count data.

**Script**

- `1. RLE_normalization_and_RUVr_function.R`

### Step 2. Preparation of exposure and control datasets

For each exposure, exposed samples were matched with the corresponding control samples according to relevant biological factors, including tissue, age, and sex. RUVr-processed data and the corresponding raw count matrices were then combined across exposure datasets for construction of the cross-center pseudo-control reference.

**Scripts**

- `2.1 combine_Exposure_RUVr_output.R`
- `2.2 combine_Exposure_raw_table.R`

### Step 3. Construction of pseudo-control references

Pseudo-controls were generated separately for female and male samples using matched controls across production centers.

For each biological condition, the **pseudo-control (CP)** represents the mean signal of the corresponding control samples across production centers. A **center control (CC)** was also calculated as the mean signal of control samples within each production center.

The difference between each center-specific control profile and the corresponding pseudo-control reference was used to characterize systematic production-center effects.

**Script**

- `3. build_pseudo_control.R`

### Step 4. Pseudo-control–based production-center correction

For each production center, a center-specific normalization factor was derived from the difference between the center control (CC) and the corresponding pseudo-control (CP).

Control samples were first aligned to the pseudo-control reference. The same center-specific adjustment was then applied to exposed samples generated at that center, thereby reducing systematic production-center effects while maintaining the relative differences between control and exposed samples.

**Script**

- `4. Run_correction_based_on_pseudo_ctrl.R`

### Step 5. Differential analysis

The corrected datasets were used for downstream differential analyses.

For RNA-seq, differentially expressed genes (DEGs) were identified using:

- [`Run_DEG_analysis.R`](https://github.com/Zhang-lab/TaRGET-II-Integrative-Analysis/blob/main/Method/Pseudo_Ctrl_correction/RNA/5.%20Run_DEG_analysis.R)

For ATAC-seq, differentially accessible regions (DARs) were identified using:

- [`Run_DAR_analysis.R`](https://github.com/Zhang-lab/TaRGET-II-Integrative-Analysis/blob/main/Method/Pseudo_Ctrl_correction/ATAC/5.%20Run_DAR_analysis.R)
