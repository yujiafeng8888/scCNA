import numpy as np
import random

def simulate_cnas_windowed(
    adata,
    window_size=100,
    cna_types=("gain", "loss", "homo_del"),
    frequencies=(0.2, 0.5, 0.9),
    chroms=("1", "2", "3"),
    cell_key="cell_type"
):
    """
    Simulate CNAs in expression matrix by modifying expression of genomic windows.
    Writes `simulated_cnvs` to adata.obs (overwrites if exists).

    Parameters:
    - adata: AnnData object
    - window_size: number of genes per simulated CNA window
    - cna_types: list of CNA types to simulate
    - frequencies: list of cell-level frequencies for each simulated CNA
    - chroms: chromosomes to simulate CNAs on
    """
    if "simulated_cnvs" not in adata.obs.columns:
        adata.obs["simulated_cnvs"] = ""
    else:
        adata.obs["simulated_cnvs"] = ""

    if not isinstance(adata.X, np.ndarray):
        X = adata.X.toarray()
    else:
        X = adata.X

    var = adata.var
    simulated_log = []

    for chrom in chroms:
        genes_chr = var[var["chromosome"] == chrom]
        gene_indices = genes_chr.index.to_numpy()

        if len(gene_indices) < window_size:
            continue

        start_idx = random.randint(0, len(gene_indices) - window_size)
        window_genes = gene_indices[start_idx:start_idx + window_size]
        window_starts = var.loc[window_genes, "start"]
        region_start = window_starts.min()
        region_end = window_starts.max()

        # Choose CNA type and frequency
        cna_type = random.choice(cna_types)
        freq = random.choice(frequencies)
        n_cells = int(len(adata) * freq)
        selected_cells = np.random.choice(adata.obs_names, n_cells, replace=False)

        # Apply effect
        for cell in selected_cells:
            idx = adata.obs_names.get_loc(cell)
            if cna_type == "gain":
                X[idx, adata.var_names.isin(window_genes)] *= 1.5
            elif cna_type == "loss":
                X[idx, adata.var_names.isin(window_genes)] *= 0.5
            elif cna_type == "homo_del":
                X[idx, adata.var_names.isin(window_genes)] = 0

            region_str = f"{chrom}:{region_start}-{region_end} (CN {0 if cna_type == 'homo_del' else ('+' if cna_type == 'gain' else '-')})"
            if adata.obs.at[cell, "simulated_cnvs"] == "":
                adata.obs.at[cell, "simulated_cnvs"] = region_str
            else:
                adata.obs.at[cell, "simulated_cnvs"] += f";{region_str}"

        simulated_log.append((chrom, region_start, region_end, cna_type, freq))

    adata.X = X
    return adata, simulated_log
