import pandas as pd

def extract_gene_position(gtf_file):
    gtf = pd.read_csv(gtf_file,sep="\t",comment="#", header=None,
                      names=["chromosome","source","feature",
                      "start","end","score","strand","frame","attribute"])
    gtf["gene_name"] = gtf["attribute"].str.extract(r'gene_name "([^"]+)"')
    gene_position = gtf[gtf["feature"] == 'gene']
    gene_position=gene_position[['gene_name','chromosome','start','end']]
    print(gene_position.columns)
    return gene_position


def annotate_gene_position(adata,gtf_filepath):
    gene_position = extract_gene_position(gtf_filepath)
    gene_position_filtered = gene_position.set_index('gene_name')
    gene_position_filtered = gene_position_filtered[~gene_position_filtered.index.duplicated(keep='first')]
    gene_position_filtered = gene_position_filtered.reindex(adata.var.index)
    
    adata.var = adata.var.join(gene_position_filtered[['chromosome', 'start', 'end']])

    return adata

def sort_genes_by_location(adata):
    """
    Assumes .var has 'chromosome' and 'start' columns.
    Sorts genes by (chromosome, start).
    """
    if 'chromosome' not in adata.var.columns or 'start' not in adata.var.columns:
        raise ValueError("adata.var must contain 'chromosome' and 'start' columns")
    
    return adata[:, adata.var.sort_values(['chromosome', 'start']).index]

