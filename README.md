# scCNA
`scCNA` is a Python package designed for detecting and annotating Copy Number Aberrations (CNA) in single-cell RNA-sequencing data. The package includes functions for identifying CNAs and annotating genes with their chromosomal positions.

## Installation

To install the `scCNA` package, you can either clone the repository or install it directly from the source. Make sure you have `pip` installed and an environment with Python 3.10+.

### Installing from Source

Clone the repository:

```bash
git clone https://github.com/yujiafeng8888/scCNA
cd scCNA
pip install .
```
## Usage

The input data of scCNA can be raw count, scCNA will normalize the count
### Annotate_gene_position

For adata does not have gene position annotation, you can use function `annotate_gene_position`.
The parameters of `find_cnas` are listed as follows:
-`adata`: input adata

-`gtf_filepath`: filepath of gtf annotation file
For adata do not have Please download a GTF annotation file. Make sure the file includes gene-level entries in the feature column. We recommend using the GENCODE annotation. For more information of GENCODE, see https://www.gencodegenes.org.

### Find_cnas
scCNA accepts both raw and normalized data, but requires cell annotations.
To identify CNAs in scRNA-seq data, you can use the function `find_cnas`.
The parameters of `find_cnas` are listed as follows:

-`adata`:inpupt adata

-`reference_key`:Column name in adata.obs that contains tumor/normal annotations.If this is set to None, the average of all cells is used as reference.

-`reference_cat`: One or multiple values in `adata.obs[reference_key]` that annotate normal cells.

-`reference`: Directly supply an array of average normal gene expression. Overrides `reference_key` and `reference_cat`.

-`min_cells`: Minimum number of cells to define a CNA.

-`threshold`: Fold change threshold.

-`window_size`: Number of genes per window.

#### Example
 ```python
import scCNA as cna
ad_def = cna.find_cnas(
    adata.copy(),
    reference_key='cell_type',
    reference_cat=[
        'CD4 T cell','CD14 monocyte','B cell','CD8 T cell',
        'NK cell','FCGR3A monocyte','Dendritic','Megakaryocyte'
    ],
    threshold=2,
    min_cells=5,
    window_size=100
)
 ```
## Referenced Code
The python module of calculate reference is adapted from icbi-lab:
Author: Gregor Sturm
Repository: https://github.com/icbi-lab/infercnvpy
License: BSD 3-Clause

## Contributing
Contributions and PRs are welcome!
## License
This project is licensed under the MIT License. See the LICENSE file for details.
## Contact
For questions or suggestions, please open an issue.






