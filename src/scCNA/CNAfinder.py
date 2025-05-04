import numpy as np
from utils import sort_genes_by_location
import scanpy as sc
from typing import Sequence, Optional
import logging
from anndata import AnnData
import os
os.environ["R_HOME"] = "/usr/lib/R"
import infercnvpy as cnv
import pandas as pd
import scipy.sparse
import re
def find_cnas(adata, min_cells=50, threshold=0.5, window_size=250,
            exclude_chromosomes: Optional[Sequence[str]] = None,
            reference_key: str | None = None,
            reference_cat: None | str | Sequence[str] = None,
            reference: np.ndarray | None = None,
            lfc_clip: float = 3,
            step: int = 10,
            dynamic_threshold: float | None = 1.5,

            ):
    """
    Detect CNAs using sliding window and annotate adata.obs['detect_CNA'].

    Parameters:
    adata
        annotated data matrix
    reference_key
        Column name in adata.obs that contains tumor/normal annotations.
        If this is set to None, the average of all cells is used as reference.
    reference_cat
        One or multiple values in `adata.obs[reference_key]` that annotate
        normal cells.
    reference
        Directly supply an array of average normal gene expression. Overrides
        `reference_key` and `reference_cat`.
    min_cells: 
        Minimum number of cells to define a CNA
    threshold: 
        Fold change threshold
    window_size: 
        Number of genes per window
    Returns:
    - cna_list: List of CNA region dictionaries
    """
    if not adata.var_names.is_unique:
        raise ValueError("Ensure your var_names are unique!")
    if {"chromosome", "start", "end"} - set(adata.var.columns) != set():
        raise ValueError(
            "Genomic positions not found. There need to be `chromosome`, `start`, and `end` columns in `adata.var`. "
        )
    var_mask = adata.var["chromosome"].isnull()
    if np.sum(var_mask):
        logging.warning(f"Skipped {np.sum(var_mask)} genes because they don't have a genomic position annotated. ")  # type: ignore
    if exclude_chromosomes is not None:
        var_mask = var_mask | adata.var["chromosome"].isin(exclude_chromosomes)

    tmp_adata = adata[:, ~var_mask]
    if reference is None:
        reference = _get_reference(adata, reference_key, reference_cat, reference)
        var_mask = np.array(var_mask)  # 如果还不是 array
        reference = reference.loc[:, ~var_mask]
        

    adata = tmp_adata.copy()
    adata = sort_genes_by_location(adata)
    genes = adata.var
    # 如果基因名是 index
    gene_names = adata.var.index.values  # numpy array of gene names

# 如果基因名是某一列，比如 "gene_name"
    cells = adata.obs_names
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    expr_norm = adata.X
    if not isinstance(expr_norm, np.ndarray):
        expr_norm = expr_norm.toarray()
    
    cna_list = []
    cell_cna_dict = {cell: [] for cell in cells}
    var = tmp_adata.var.loc[:, ["chromosome", "start", "end"]]
    for celltype in reference.index:
        selected_cells = adata.obs_names[adata.obs['cell_type'] == celltype]
        for chrom in genes['chromosome'].unique():
            idx = genes['chromosome'] == chrom
            gene_subset = genes[idx]
            # 获取 selected_cells 对应的整数索引
            selected_cell_indices = np.array([np.where(adata.obs_names == cell)[0][0] for cell in selected_cells])
            # 选择与当前染色体相关的基因表达数据
            expr_subset = expr_norm[selected_cell_indices,:]
            expr_subset=expr_subset[:, idx.values]

            num_genes = expr_subset.shape[1]
            for start in range(0, num_genes, window_size):
                end = min(start + window_size, num_genes)
                if end - start < 50:
                    continue
                ref_profile = reference.loc[celltype]  # shape: (n_genes,) or (n_cells, n_genes)
                ref_window = ref_profile.iloc[:, start:end] if ref_profile.ndim == 2 else ref_profile[start:end]
                genes_list= gene_names[start:end]
                smothed_windoe_expr=_infercnv_chunk(expr_subset[:, start:end],var,ref_window,lfc_clip,window_size,step,dynamic_threshold,genes_list)
                smothed_windoe_expr=smothed_windoe_expr.toarray()
                # print(smothed_windoe_expr)
                # print(np.all(smothed_windoe_expr == 0))
                window_expr = np.mean(smothed_windoe_expr,axis=1)
                # print(window_expr.shape)
                # 获取指定细胞类型的列
                # ref_profile = reference.loc[celltype]  # shape: (n_genes,) or (n_cells, n_genes)
                # ref_window = ref_profile.iloc[:, start:end] if ref_profile.ndim == 2 else ref_profile[start:end]
                reference_window_expr = np.mean(ref_window)

                # Step 2: Compute fold-change relative to reference (this will identify relative expression levels)
    
                # window_expr = np.mean(expr_subset[:, start:end], axis=1)
                # fold_change = window_expr / reference_window_expr  # Calculate fold change against reference
                
                # gain_mask = fold_change > threshold  # Gain if fold-change is higher than threshold
                # loss_mask = fold_change < 1 / threshold
                gain_mask=window_expr>threshold
                loss_mask=window_expr<-threshold
                gain_mask = gain_mask.flatten()
                loss_mask=loss_mask.flatten()
                # print(selected_cells.shape)
                # print(gain_mask.shape)
                gain_cells = selected_cells[gain_mask]
                loss_cells = selected_cells[loss_mask]
                window_start_pos = gene_subset.iloc[start]['start']
                window_end_pos = gene_subset.iloc[end - 1]['start']
                window_label_gain = f"{chrom}:{window_start_pos}-{window_end_pos}_gain"
                window_label_loss = f"{chrom}:{window_start_pos}-{window_end_pos}_loss"

                if len(gain_cells) >= min_cells:
                    cna_list.append({'chromosome': chrom, 'start': window_start_pos, 'stop': window_end_pos, 'type': 'gain'})
                    for cell in gain_cells:
                        cell_cna_dict[cell].append(window_label_gain)

                if len(loss_cells) >= min_cells:
                    cna_list.append({'chromosome': chrom, 'start': window_start_pos, 'stop': window_end_pos, 'type': 'loss'})
                    for cell in loss_cells:
                        cell_cna_dict[cell].append(window_label_loss)
    # 写入 adata.obs['detect_CNA']
    adata.obs['detect_CNA'] = adata.obs_names.map(lambda x: ';'.join(cell_cna_dict[x]) if cell_cna_dict[x] else 'none')
    return adata

def _get_reference(
    adata: AnnData,
    reference_key: str | None,
    reference_cat: None | str | Sequence[str],
    reference: np.ndarray | None,
) -> np.ndarray:
    """Parameter validation extraction of reference gene expression.

    If multiple reference categories are given, compute the mean per
    category.

    Returns a 2D array with reference categories in rows, cells in columns.
    If there's just one category, it's still a 2D array.
    """
    if reference is None:
        if reference_key is None or reference_cat is None:
            logging.warning(
                "Using mean of all cells as reference. For better results, "
                "provide either `reference`, or both `reference_key` and `reference_cat`. "
            )  # type: ignore
            reference = np.mean(adata.X, axis=0)

        else:
            obs_col = adata.obs[reference_key]
            if isinstance(reference_cat, str):
                reference_cat = [reference_cat]
            reference_cat = np.array(reference_cat)
            reference_cat_in_obs = np.isin(reference_cat, obs_col)
            if not np.all(reference_cat_in_obs):
                raise ValueError(
                    "The following reference categories were not found in "
                    "adata.obs[reference_key]: "
                    f"{reference_cat[~reference_cat_in_obs]}"
                )
            reference = np.vstack([np.mean(adata.X[obs_col.values == cat, :], axis=0) for cat in reference_cat])
            reference = pd.DataFrame(reference,
                    index=reference_cat,            # 设置行为 cell type 名
                    columns=adata.var_names)
    if reference.ndim == 1:
        reference = reference[np.newaxis, :]

    if reference.shape[1] != adata.shape[1]:
        raise ValueError("Reference must match the number of genes in AnnData. ")

    return reference

def _infercnv_chunk(tmp_x, var, reference, lfc_cap, window_size, step, dynamic_threshold,genes):
    """The actual infercnv work is happening here.

    Process chunks of serveral thousand genes independently since this
    leads to (temporary) densification of the matrix.

    Parameters see `infercnv`.
    tmp_x:
    partial expr gene list: [expr[i : i + chunksize, :] for i in range(0, adata.shape[0], chunksize)]
    """
    # Step 1 - compute log fold change. This densifies the matrix.
    # Per default, use "bounded" difference calculation, if multiple references
    # are available. This mitigates cell-type specific biases (HLA, IGHG, ...)
    reference = pd.DataFrame([reference]) 
    # if reference.shape[0] == 1:
    x_centered = tmp_x - reference.to_numpy()[0, :]
    # else:
    #     ref_min = np.min(reference, axis=0)
    #     ref_max = np.max(reference, axis=0)
    #     # entries that are between the two "bounds" are considered having a logFC of 0.
    #     x_centered = np.zeros(tmp_x.shape, dtype=tmp_x.dtype)
    #     above_max = tmp_x > ref_max
    #     below_min = tmp_x < ref_min
    #     x_centered[above_max] = _ensure_array(tmp_x - ref_max)[above_max]
    #     x_centered[below_min] = _ensure_array(tmp_x - ref_min)[below_min]

    x_centered = _ensure_array(x_centered)
    # Step 2 - clip log fold changes
    x_clipped = np.clip(x_centered, -lfc_cap, lfc_cap)
    # Step 3 - smooth by genomic position
    chr_pos, x_smoothed = _running_mean_by_chromosome(
        x_clipped, var, window_size=window_size, step=step,genes=genes)

    x_res = x_clipped - np.median(x_smoothed, axis=1)[:, np.newaxis]
    # print("x_smoothed[:5]", x_smoothed[:5])
    # print("median:", np.median(x_smoothed, axis=1)[:5])
    # step 5 - standard deviation based noise filtering
    if dynamic_threshold is not None:
        noise_thres = dynamic_threshold * np.std(x_res)
        x_res[np.abs(x_res) < noise_thres] = 0
    x_res = scipy.sparse.csr_matrix(x_res)
    return x_res

def _ensure_array(a):
    """If a is a matrix, turn it into an array."""
    if isinstance(a, np.matrix):
        return a.A
    else:
        return a

def _running_mean(
    x: np.ndarray | scipy.sparse.spmatrix,
    n: int = 50,
    step: int = 10,
    gene_list: list = None
) -> tuple[np.ndarray, pd.DataFrame | None]:
    """
    Compute a pyramidially weighted running mean.

    Densifies the matrix. Use `step` and `chunksize` to save memory.

    Parameters
    ----------
    x
        matrix to work on
    n
        Length of the running window
    step
        only compute running windows every `step` columns, e.g. if step is 10
        0:99, 10:109, 20:119 etc. Saves memory.
    gene_list
        List of gene names to be used in the convolution
    """
    if n < x.shape[1]:  # regular convolution: the filter is smaller than the #genes
        r = np.arange(1, n + 1)
        pyramid = np.minimum(r, r[::-1])
        smoothed_x = np.apply_along_axis(
            lambda row: np.convolve(row, pyramid, mode="valid"),
            axis=1,
            arr=x,) / np.sum(pyramid)

        ## get the indices of the genes used in the convolution
        convolution_indices = get_convolution_indices(x, n)[np.arange(0, smoothed_x.shape[1], step)]
        ## Pull out the genes used in the convolution
        convolved_gene_names = gene_list[convolution_indices]
        smoothed_x = smoothed_x[:, np.arange(0, smoothed_x.shape[1], step)]

        return smoothed_x

    else:  # If there is less genes than the window size, set the window size to the number of genes and perform a single convolution
        n = x.shape[1]  # set the filter size to the number of genes
        r = np.arange(1, n + 1)
        ## As we are only doing one convolution the values should be equal
        pyramid = np.array([1] * n)
        smoothed_x = np.apply_along_axis(
            lambda row: np.convolve(row, pyramid, mode="valid"),
            axis=1,
            arr=x,
        ) / np.sum(pyramid)
        return smoothed_x

def get_convolution_indices(x, n):
    indices = []
    for i in range(x.shape[1] - n + 1):
        indices.append(np.arange(i, i + n))
    return np.array(indices)

def _running_mean_by_chromosome(
    expr, var, window_size, step, genes
) -> tuple[dict, np.ndarray, pd.DataFrame | None]:
    """Compute the running mean for each chromosome independently. Stack the resulting arrays ordered by chromosome.

    Parameters
    ----------
    expr
        A gene expression matrix, appropriately preprocessed
    var
        The var data frame of the associated AnnData object
    window_size
        size of the running window (number of genes in to include in the window)
    step
        only compute every nth running window where n = `step`. Set to 1 to compute
        all windows.

    Returns
    -------
    chr_start_pos
        A Dictionary mapping each chromosome to the index of running_mean where
        this chromosome begins.
    running_mean
        A numpy array with the smoothed gene expression, ordered by chromosome
        and genomic position
    """
    chromosomes = _natural_sort([x for x in var["chromosome"].unique()])
    # print([x for x in var["chromosome"].unique() if x != "chrM"])
    running_means = [ _running_mean_for_chromosome(chr, expr, var, window_size, step,genes) for chr in chromosomes
    ]
    chr_start_pos = {}
    for chr, i in zip(chromosomes, np.cumsum([0] + [x.shape[1] for x in running_means]), strict=False):
        chr_start_pos[chr] = i

    ## Concatenate the gene dfs

    return chr_start_pos, np.hstack(running_means)

def _natural_sort(l: Sequence):
    """Natural sort without third party libraries.
    Adapted from: https://stackoverflow.com/a/4836734/2340703
    """
    def convert(text):
        return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    return sorted(l, key=alphanum_key)

def _running_mean_for_chromosome(chr, expr, var, window_size, step,genes):
    # genes = var.loc[var["chromosome"] == chr].sort_values("start").index.values
    # print('gene',genes.shape)
    # print(expr.shape)
    # tmp_x = expr[:, var.index.get_indexer(genes)]
    x_conv= _running_mean(
        expr, n=window_size, step=step, gene_list=genes
    )
    return x_conv

if __name__ == '__main__':
    import multiprocessing
    multiprocessing.set_start_method("spawn")
    adata=sc.read_h5ad("../../PBMC_simulated_cnas_041025.h5ad")
    # import numpy as np
    # chrom_counts = adata.var["chromosome"].value_counts()
    # single_gene_chroms = chrom_counts[chrom_counts == 1].index.tolist()
    # print("只有一个基因的染色体有：", single_gene_chroms)
    # cnv.tl.infercnv(
    # adata,
    # reference_key="cell_type",
    # exclude_chromosomes= ['HSCHR20_1_CTG3', 'HSCHR19_4_CTG2', 'HSCHR19_2_CTG3_1', 'GL000009.2', 'HSCHR21_3_CTG1_1', 'HSCHR21_6_CTG1_1', 'HSCHR22_1_CTG3', 'HSCHR18_3_CTG2_1', 'HSCHR22_1_CTG4', 'HSCHR22_3_CTG1', 'HSCHRX_2_CTG3', 'KI270734.1', 'HSCHR18_5_CTG1_1', 'GL000218.1', 'HSCHR18_2_CTG2_1', 'HSCHR18_2_CTG2', 'HG926_PATCH', 'HG1343_HG173_HG459_PATCH', 'HG2213_PATCH', 'HG2521_PATCH', 'HSCHR2_1_CTG5', 'HSCHR2_4_CTG1', 'HSCHR2_6_CTG1', 'HSCHR3_9_CTG2_1', 'GL000195.1', 'HSCHR7_1_CTG4_4', 'HSCHR7_2_CTG1', 'HSCHR7_3_CTG6', 'HSCHR10_1_CTG2', 'HSCHR13_1_CTG1', 'GL000219.1', 'HSCHR17_3_CTG1', 'HSCHR17_3_CTG2', 'HSCHR17_7_CTG4', 'HSCHR18_1_CTG1_1', 'HSCHR16_2_CTG3_1'],
    # reference_cat=[
    #             'CD4 T cell',
    #             'CD14 monocyte',
    #             'B cell',          
    #             'CD8 T cell',
    #             'NK cell', 
    #             'FCGR3A monocyte',   
    #             'Dendritic',
    #             'Megakaryocyte'],
    # window_size=250,
    # )
    # print(adata.obs.columns)
# no_cnv_cells = adata.obs[
#     (adata.obs['simulated_cnvs'].isna()) | 
#     (adata.obs['simulated_cnvs'] == '') | 
#     (adata.obs['simulated_cnvs'] == 'none')
# ].index
# var_mask = adata.var["chromosome"].isnull()
# tmp_adata = adata[:, ~var_mask]
# adata=tmp_adata.copy()
# adata_no_cnv = adata[no_cnv_cells]
# # 获取表达矩阵 (cells × genes)
# expr = pd.DataFrame(adata_no_cnv.X.toarray() if hasattr(adata_no_cnv.X, 'toarray') else adata_no_cnv.X,
#                     index=adata_no_cnv.obs_names,
#                     columns=adata_no_cnv.var_names)  # 检查重复的列
# # 加入细胞类型信息
# expr['cell_type'] = adata_no_cnv.obs['cell_type'].values

# # 分组计算每个基因的均值
# referencs = expr.groupby('cell_type').mean().drop(columns='cell_type', errors='ignore')
adata=find_cnas(adata,reference_key="cell_type",
                reference_cat=[
                'CD4 T cell',
                'CD14 monocyte',
                'B cell',          
                'CD8 T cell',
                'NK cell', 
                'FCGR3A monocyte',   
                'Dendritic',
                'Megakaryocyte'],window_size=250)
# print(adata.obs.columns)
print(adata.obs['detect_CNA'].value_counts())
# adata.var.loc[:, ["ensg", "chromosome", "start", "end"]].head()
