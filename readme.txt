CCI v1.0 (Counterpart Composite Index) was designed for projecting carcinoma cells into reference cells, to identify their healthy counterparts that could reveal the subtypes and lineage stage of carcinoma cells based on single cell RNA sequencing data (scRNA-seq).


###Installation###
Please compile cpp codes of CCI by GNU c++, to do this just run shells

   bash makefile.sh

Compilation should be finished in seconds.CCI also need R >3.4 and several R packages available in your computer including Seurat, ggplot2 and RColorBrewer.


###Running instructions###
Please implement CCI calculation following below steps:

1.Construct normal cell atlas as reference.
Reference cells were stored in folder data_ref/ including 4 files as expression matrix (cell as rows and gene as columns), gene list, cell types and identity of each cells as in expression matrix. We provide a comprehensive cell atlas of normal BMMCs here, which could be used for study of hematopoiesis and leukemogenesis.

2.Cluster your target cells
CCI projects cell clusters to reference instead of single cell due to prevalent droplets in scRNA-seq. Please properly cluster your target cells before running CCI.
Target cells are stored in folder data_obj/ with expression matrix files of each sub-cluster. 

3.Run CCI.
To calculation CCI score, please simplely run 

   bash run_cci.sh AML1_cluster1 (toy data in data_obj/) 

For the toy data which contains 100 cells and 4K genes, computating will be finished in several minutes by one normal CPU spread.
CCI score (BF probability) could be found in folder result/.

4.Visualization.
Users should visualize the CCI results by their own. We provide example plots as PDFs in result/ folder.

Any questions please feel free to contact qinpf@sustech.edu.cn




