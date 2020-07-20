CCI v1.0 (Counterpart Composite Index) was designed for projecting carcinoma cells into reference cells, to identify their healthy counterparts that could reveal the subtypes and lineage stage of carcinoma cells based on single cell RNA sequencing data (scRNA-seq).
Details of CCI were described in paper Integrated decoding hematopoiesis and leukemogenesis at single-cell resolution and its clinical implication. (Qin et al., doi: https://doi.org/10.1101/2020.02.21.960401)

###Installation###
Please compile cpp codes of CCI by GNU c++, to do this just run shells

   bash makefile.sh

Compilation should be finished in seconds.CCI also need R >3.4 and several R packages available in your computer including Seurat, ggplot2 and RColorBrewer.


###Running instructions###
Please implement CCI calculation following steps below:

1.Construct normal cell atlas as reference
Reference cells were stored in folder data_ref/ including 4 necessary files as 1)expression matrix (cell as rows and gene as columns), 2)gene list, 3)cell types and 4)identity of each cell as in expression matrix. Please format your files as the toy data. We provide a comprehensive cell atlas of normal BMMCs here, which could be used for study of hematopoiesis and leukemogenesis.

2.Cluster your target cells
CCI projects cell clusters to reference instead of single cell due to prevalent droplets in scRNA-seq. Please properly cluster your target cells before running CCI.
Target cells are stored in folder data_obj/ with expression matrix file of each sub-cluster. 

3.Run CCI
To calculate CCI score, please simplely run 

   bash run_cci.sh AML1_cluster1 (toy data in data_obj/) 

For the toy data which contains 100 cells and 4K genes, computating will be finished in several minutes by one normal CPU spread.
CCI score (BF probability) could be found in folder result/.

4.Visualization.
Users should visualize CCI results by their own. We provide example plots as PDFs in result/ folder.

Any question please feel free to contact qinpf@sustech.edu.cn




