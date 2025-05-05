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
## usage

The input data of scCNA can be raw count, scCNA will normalize the count
## annotate_gene_position
Please download a GTF annotation file. Make sure the file includes gene-level entries in the feature column. We recommend using the GENCODE annotation. 
## find_cnas
## Referenced Code
The python module of calculate reference and window mean expression is adapted from icbi-lab:
Author: Gregor Sturm
Repository: https://github.com/icbi-lab/infercnvpy
License: BSD 3-Clause

## Contributing
Contributions and PRs are welcome!
## License
This project is licensed under the MIT License. See the LICENSE file for details.
## Contact
For questions or suggestions, please open an issue.



