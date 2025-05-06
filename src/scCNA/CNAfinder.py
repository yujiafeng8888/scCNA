import numpy as np
from .utils import sort_genes_by_location
import scanpy as sc
from typing import Sequence
import logging
from anndata import AnnData
import os
os.environ["R_HOME"] = "/usr/lib/R"
from scipy.sparse import csr_matrix

def find_cnas(adata, min_cells=25, threshold=10, window_size=250,exclude_chromosomes: Sequence[str] | None = ("chrX", "chrY"),
            reference_key: str | None = None,
            reference_cat: None | str | Sequence[str] = None,
            reference: np.ndarray | None = None,):
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
    import numpy as np
    import pandas as pd
    import scanpy as sc
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

    tmp_adata = adata[:, ~var_mask].copy()
    if isinstance(tmp_adata.X, np.ndarray):
        tmp_adata.X = csr_matrix(tmp_adata.X)
        print('converting to sparse matrix...')
    adata = tmp_adata.copy()
    reference = _get_reference(adata, reference_key, reference_cat, reference)
    adata = sort_genes_by_location(adata)
    genes = adata.var
    cells = adata.obs_names
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)

    expr_norm = adata.X

    cna_list = []
    cell_cna_dict = {cell: [] for cell in cells}
    for i in range(0,len(reference_cat)):
        cell_type = reference_cat[i] # select celltype
        selected_cells = adata.obs_names[adata.obs['cell_type'] == cell_type]
        for chrom in genes['chromosome'].unique():
            idx = genes['chromosome'] == chrom
            gene_subset = genes[idx]
            # get selected_cells interger index
            selected_cell_indices = np.array([np.where(adata.obs_names == cell)[0][0] for cell in selected_cells])

            # Select gene expression data related to the current chromosome.
            expr_subset = expr_norm[selected_cell_indices,:]
            expr_subset=expr_subset[:, idx.values]

            num_genes = expr_subset.shape[1]
            for start in range(0, num_genes, window_size):
                end = min(start + window_size, num_genes)
                if end - start < 50:
                    continue

                window_expr = np.mean(expr_subset[:, start:end], axis=1)
                
                part_reference=reference[i, start:end]
                print(part_reference.shape)
                reference_window_expr = np.mean(part_reference,axis=1)
                
                # Step 2: Compute fold-change relative to reference (this will identify relative expression levels)
                fold_change = window_expr / reference_window_expr  # Calculate fold change against reference

                # Step 3: Detect gain and loss based on fold-change
                gain_mask = fold_change > threshold  # Gain if fold-change is higher than threshold
                loss_mask = fold_change < 1 / threshold
                gain_mask = gain_mask.A.flatten()
                loss_mask=loss_mask.A.flatten()
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
    # write adata.obs['detect_CNA']
    # print(cell_cna_dict)
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

    if reference.ndim == 1:
        reference = reference[np.newaxis, :]

    if reference.shape[1] != adata.shape[1]:
        raise ValueError("Reference must match the number of genes in AnnData. ")

    return reference

# adata=sc.read_h5ad("../../adata.sim.h5ad")
# adata=find_cnas(adata,reference_key="cell_type",
#                 reference_cat=[
#                 'CD4 T cell',
#                 'CD14 monocyte',
#                 'B cell',          
#                 'CD8 T cell',
#                 'NK cell', 
#                 'FCGR3A monocyte',   
#                 'Dendritic',
#                 'Megakaryocyte'])

# print(adata.obs.columns)
# print(adata.obs['cell_type'].value_counts())
# adata.var.loc[:, ["ensg", "chromosome", "start", "end"]].head()
# print(adata.X.shape)
