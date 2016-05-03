April 4, 2016

Bioconductors:

We are pleased to announce Bioconductor 3.3, consisting of 1211
software packages, 293 experiment data packages, and 916
up-to-date annotation packages.

There are 108 new software packages, and many updates and improvements
to existing packages; Bioconductor 3.3 is compatible with R 3.3,
and is supported on Linux, 32- and 64-bit Windows, and Mac OS X.  This
release includes an updated Bioconductor [Amazon Machine Image][1]
and [Docker containers][2].

Visit [http://bioconductor.org][3]
for details and downloads.

[1]: http://bioconductor.org/help/bioconductor-cloud-ami/
[2]: http://bioconductor.org/help/docker/
[3]: http://bioconductor.org

Contents
--------

* Getting Started with Bioconductor 3.3
* New Software Packages
* NEWS from new and existing packages
* Packages removed from Bioconductor since the last release

Getting Started with Bioconductor 3.3
======================================

To update to or install Bioconductor 3.3:

1. Install R 3.3.  Bioconductor 3.3 has been designed expressly
for this version of R.

2. Follow the instructions at
[http://bioconductor.org/install/](http://bioconductor.org/install/) .

New Software Packages
=====================

There are 108 new packages in this release of Bioconductor.

AneuFinder - This package implements functions for CNV calling, plotting, export and analysis from whole-genome single cell sequencing data.

bacon - Bacon can be used to remove inflation and bias often observed in epigenome- and transcriptome-wide association studies. To this end bacon constructs an empirical null distribution using a Gibbs Sampling algorithm by fitting a three-component normal mixture on z-scores.

BadRegionFinder - BadRegionFinder is a package for identifying regions with a bad, acceptable and good coverage in sequence alignment data available as bam files. The whole genome may be considered as well as a set of target regions. Various visual and textual types of output are available.

BasicSTARRseq - Basic peak calling on STARR-seq data based on a method introduced in "Genome-Wide Quantitative Enhancer Activity Maps Identified by STARR-seq" Arnold et al. Science. 2013 Mar 1;339(6123):1074-7. doi: 10.1126/science. 1232542. Epub 2013 Jan 17.

BatchQC - Sequencing and microarray samples often are collected or processed in multiple batches or at different times. This often produces technical biases that can lead to incorrect results in the downstream analysis. BatchQC is a software tool that streamlines batch preprocessing and evaluation by providing interactive diagnostics, visualizations, and statistical analyses to explore the extent to which batch variation impacts the data. BatchQC diagnostics help determine whether batch adjustment needs to be done, and how correction should be applied before proceeding with a downstream analysis. Moreover, BatchQC interactively applies multiple common batch effect approaches to the data, and the user can quickly see the benefits of each method. BatchQC is developed as a Shiny App. The output is organized into multiple tabs, and each tab features an important part of the batch effect analysis and visualization of the data. The BatchQC interface has the following analysis groups: Summary, Differential Expression, Median Correlations, Heatmaps, Circular Dendrogram, PCA Analysis, Shape, ComBat and SVA.

BgeeDB - A package for the annotation and gene expression data download from Bgee database, and TopAnat analysis: GO-like enrichment of anatomical terms, mapped to genes by expression patterns.

biomformat - This is an R package for interfacing with the BIOM format. This package includes basic tools for reading biom-format files, accessing and subsetting data tables from a biom object (which is more complex than a single table), as well as limited support for writing a biom-object back to a biom-format file. The design of this API is intended to match the python API and other tools included with the biom-format project, but with a decidedly "R flavor" that should be familiar to R users. This includes S4 classes and methods, as well as extensions of common core functions/methods.

BioQC - BioQC performs quality control of high-throughput expression data based on tissue gene signatures

biosigner - Feature selection is critical in omics data analysis to extract restricted and meaningful molecular signatures from complex and high-dimension data, and to build robust classifiers. This package implements a new method to assess the relevance of the variables for the prediction performances of the classifier. The approach can be run in parallel with the PLS-DA, Random Forest, and SVM binary classifiers. The signatures and the corresponding 'restricted' models are returned, enabling future predictions on new datasets. A Galaxy implementation of the package is available within the Workflow4metabolomics.org online infrastructure for computational metabolomics.

cellity - A support vector machine approach to identifying and filtering low quality cells from single-cell RNA-seq datasets.

cellTree - This packages computes a Latent Dirichlet Allocation (LDA) model of single-cell RNA-seq data and builds a compact tree modelling the relationship between individual cells over time or space.

Chicago - A pipeline for analysing Capture Hi-C data.

chromPlot - Package designed to visualize genomic data along the chromosomes, where the vertical chromosomes are sorted by number, with sex chromosomes at the end.

CHRONOS - A package used for efficient unraveling of the inherent dynamic properties of pathways. MicroRNA-mediated subpathway topologies are extracted and evaluated by exploiting the temporal transition and the fold change activity of the linked genes/microRNAs.

CINdex - The CINdex package addresses important area of high-throughput genomic analysis. It allows the automated processing and analysis of the experimental DNA copy number data generated by Affymetrix SNP 6.0 arrays or similar high throughput technologies. It calculates the chromosome instability (CIN) index that allows to quantitatively characterize genome-wide DNA copy number alterations as a measure of chromosomal instability. This package calculates not only overall genomic instability, but also instability in terms of copy number gains and losses separately at the chromosome and cytoband level.

clustComp - clustComp is a package that implements several techniques for the comparison and visualisation of relationships between different clustering results, either flat versus flat or hierarchical versus flat. These relationships among clusters are displayed using a weighted bi-graph, in which the nodes represent the clusters and the edges connect pairs of nodes with non-empty intersection; the weight of each edge is the number of elements in that intersection and is displayed through the edge thickness. The best layout of the bi-graph is provided by the barycentre algorithm, which minimises the weighted number of crossings. In the case of comparing a hierarchical and a non-hierarchical clustering, the dendrogram is pruned at different heights, selected by exploring the tree by depth-first search, starting at the root. Branches are decided to be split according to the value of a scoring function, that can be based either on the aesthetics of the bi-graph or on the mutual information between the hierarchical and the flat clusterings. A mapping between groups of clusters from each side is constructed with a greedy algorithm, and can be additionally visualised.

ClusterSignificance - The ClusterSignificance package provides tools to assess if clusters have a separation different from random or permuted data. ClusterSignificance investigates clusters of two or more groups by first, projecting all points onto a one dimensional line. Cluster separations are then scored and the probability of the seen separation being due to chance is evaluated using a permutation method.

CONFESS - Single Cell Fluidigm Spot Detector.

consensusSeekeR - This package compares genomic positions and genomic ranges from multiple experiments to extract common regions. The size of the analyzed region is adjustable as well as the number of experiences in which a feature must be present in a potential region to tag this region as a consensus region.

contiBAIT - Using strand inheritance data from multiple single cells from the organism whose genome is to be assembled, contiBAIT can cluster unbridged contigs together into putative chromosomes, and order the contigs within those chromosomes.

CountClust - Fits grade of membership models (GoM, also known as admixture models) to cluster RNA-seq gene expression count data, identifies characteristic genes driving cluster memberships, and provides a visual summary of the cluster memberships.

CrispRVariants - CrispRVariants provides tools for analysing the results of a CRISPR-Cas9 mutagenesis sequencing experiment, or other sequencing experiments where variants within a given region are of interest. These tools allow users to localize variant allele combinations with respect to any genomic location (e.g. the Cas9 cut site), plot allele combinations and calculate mutation rates with flexible filtering of unrelated variants.

dada2 - The dada2 package provides "OTU picking" functionality, but instead of picking OTUs the DADA2 algorithm exactly infers samples sequences. The dada2 pipeline starts from demultiplexed fastq files, and outputs inferred sample sequences and associated abundances after removing substitution and chimeric errors. Taxonomic classification is also available via a native implementation of the RDP classifier method.

dcGSA - Distance-correlation based Gene Set Analysis for longitudinal gene expression profiles. In longitudinal studies, the gene expression profiles were collected at each visit from each subject and hence there are multiple measurements of the gene expression profiles for each subject. The dcGSA package could be used to assess the associations between gene sets and clinical outcomes of interest by fully taking advantage of the longitudinal nature of both the gene expression profiles and clinical outcomes.

debrowser - Bioinformatics platform containing interactive plots and tables for differential gene and region expression studies. Allows visualizing expression data much more deeply in an interactive and faster way. By changing the parameters, user can easily discover different parts of the data that like never have been done before. Manually creating and looking these plots takes time. With this system users can prepare plots without writing any code. Differential expression, PCA and clustering analysis are made on site and the results are shown in various plots such as scatter, bar, box, volcano, ma plots and Heatmaps.

DEFormats - Covert between different data formats used by differential gene expression analysis tools.

diffloop - A suite of tools for subsetting, visualizing, annotating, and statistically analyzing the results of one or more ChIA-PET experiments.

DNAshapeR - DNAhapeR is an R/BioConductor package for ultra-fast, high-throughput predictions of DNA shape features. The package allows to predict, visualize and encode DNA shape features for statistical learning.

doppelgangR - The main function is doppelgangR(), which takes as minimal input a list of ExpressionSet object, and searches all list pairs for duplicated samples.  The search is based on the genomic data (exprs(eset)), phenotype/clinical data (pData(eset)), and "smoking guns" - supposedly unique identifiers found in pData(eset).

DRIMSeq - The package provides two frameworks. One for the differential splicing analysis between different conditions and one for the sQTL analysis. Both are based on modeling the counts of genomic features (i.e., transcripts, exons or exonic bins) with Dirichlet-multinomial distribution. The package also makes available functions for visualization and exploration of the data and results.

EBSEA - Calculates differential expression of genes based on exon counts of genes obtained from RNA-seq sequencing data.

EGAD - The package implements a series of highly efficient tools to calculate functional properties of networks based on guilt by association methods.

EGSEA - This package implements the Ensemble of Gene Set Enrichment Analyses (EGSEA) method for gene set testing.

EmpiricalBrownsMethod - Combining P-values from multiple statistical tests is common in bioinformatics. However, this procedure is non-trivial for dependent P-values. This package implements an empirical adaptation of Brown’s Method (an extension of Fisher’s Method) for combining dependent P-values which is appropriate for highly correlated data sets found in high-throughput biological experiments.

epivizrData - Serve data from Bioconductor Objects through a WebSocket connection.

epivizrServer - This package provides objects to manage WebSocket connections to epiviz apps. Other epivizr package use this infrastructure.

epivizrStandalone - This package imports the epiviz visualization JavaScript app for genomic data interactive visualization. The 'epivizrServer' package is used to provide a web server running completely within R. This standalone version allows to browse arbitrary genomes through genome annotations provided by Bioconductor packages.

EWCE - Used to determine which cell types are enriched within gene lists. The package provides tools for testing enrichments within simple gene lists (such as human disease associated genes) and those resulting from differential expression studies. The package does not depend upon any particular Single Cell Transcriptome dataset and user defined datasets can be loaded in and used in the analyses.

ExpressionAtlas - This package is for searching for datasets in EMBL-EBI Expression Atlas, and downloading them into R for further analysis. Each Expression Atlas dataset is represented as a SimpleList object with one element per platform. Sequencing data is contained in a SummarizedExperiment object, while microarray data is contained in an ExpressionSet or MAList object.

FamAgg - Framework providing basic pedigree analysis and plotting utilities as well as a variety of methods to evaluate familial aggregation of traits in large pedigrees.

flowAI - The package is able to perform an automatic or interactive quality control on FCS data acquired using flow cytometry instruments. By evaluating three different properties: 1) flow rate, 2) signal acquisition, 3) dynamic range, the quality control enables the detection and removal of anomalies.

garfield - GARFIELD is a non-parametric functional enrichment analysis approach described in the paper GARFIELD: GWAS analysis of regulatory or functional information enrichment with LD correction. Briefly, it is a method that leverages GWAS findings with regulatory or functional annotations (primarily from ENCODE and Roadmap epigenomics data) to find features relevant to a phenotype of interest. It performs greedy pruning of GWAS SNPs (LD r2 > 0.1) and then annotates them based on functional information overlap. Next, it quantifies Fold Enrichment (FE) at various GWAS significance cutoffs and assesses them by permutation testing, while matching for minor allele frequency, distance to nearest transcription start site and number of LD proxies (r2 > 0.8).

genbankr - Reads Genbank files.

GenoGAM - This package allows statistical analysis of genome-wide data with smooth functions using generalized additive models based on the implementation from the R-package 'mgcv'. It provides methods for the statistical analysis of ChIP-Seq data including inference of protein occupancy, and pointwise and region-wise differential analysis. Estimation of dispersion and smoothing parameters is performed by cross-validation. Scaling of generalized additive model fitting to whole chromosomes is achieved by parallelization over overlapping genomic intervals.

genphen - Given a set of genetic polymorphisms in the form of single nucleotide poylmorphisms or single amino acid polymorphisms and a corresponding phenotype data, often we are interested to quantify their association such that we can identify the causal polymorphisms. Using statistical learning techniques such as random forests and support vector machines, this tool provides the means to estimate genotype-phenotype associations. It also provides visualization functions which enable the user to visually inspect the results of such genetic association study and conveniently select the genotypes which have the highest strenght ofassociation with the phenotype.

GenRank - Methods for ranking genes based on convergent evidence obtained from multiple independent evidence layers. This package adapts three methods that are popular for meta-analysis.

GenVisR - Produce highly customizable publication quality graphics for genomic data primarily at the cohort level.

ggcyto - With the dedicated fority method implemented for flowSet, ncdfFlowSet and GatingSet classes, both raw and gated flow cytometry data can be plotted directly with ggplot. ggcyto wrapper and some customed layers also make it easy to add gates and population statistics to the plot.

Glimma - This package generates interactive visualisations of RNA-sequencing data based on output from limma, edgeR or DESeq2. Interactions are built on top of popular static displays from the limma package, providing users with access to gene IDs and sample information. Plots are generated using d3.js and displayed in HTML pages.

globalSeq - The method may be conceptualised as a test of overall significance in regression analysis, where the response variable is overdispersed and the number of explanatory variables exceeds the sample size.

GMRP - Perform Mendelian randomization analysis of multiple SNPs to determine risk factors causing disease of study and to exclude confounding variabels and perform path analysis to construct path of risk factors to the disease.

GSALightning - GSALightning provides a fast implementation of permutation-based gene set analysis for two-sample problem. This package is particularly useful when testing simultaneously a large number of gene sets, or when a large number of permutations is necessary for more accurate p-values estimation.

Harman - Harman is a PCA and constrained optimisation based technique that maximises the removal of batch effects from datasets, with the constraint that the probability of overcorrection (i.e. removing genuine biological signal along with batch noise) is kept to a fraction which is set by the end-user.

HDF5Array - This package implements the HDF5Array class for convenient access and manipulation of HDF5 datasets. In order to reduce memory usage and optimize performance, operations on an HDF5Array object are either delayed or executed using a block processing mechanism. The delaying and block processing mechanisms are independent of the on-disk backend and implemented via the DelayedArray class. They even work on ordinary arrays where they can sometimes improve performance.

iCARE - An R package to compute Individualized Coherent Absolute Risk Estimators.

iCOBRA - This package provides functions for calculation and visualization of performance metrics for evaluation of ranking and binary classification (assignment) methods. It also contains a shiny application for interactive exploration of results.

IHW - Independent hypothesis weighting (IHW) is a multiple testing procedure that increases power compared to the method of Benjamini and Hochberg by assigning data-driven weights to each hypothesis. The input to IHW is a two-column table of p-values and covariates. The covariate can be any continuous-valued or categorical variable that is thought to be informative on the statistical properties of each hypothesis test, while it is independent of the p-value under the null hypothesis.

ImmuneSpaceR - Provides a convenient API for accessing data sets within ImmuneSpace (www.immunespace.org), the data repository and analysis platform of the Human Immunology Project Consortium (HIPC).

InteractionSet - Provides the GInteractions, InteractionSet and ContactMatrix objects and associated methods for storing and manipulating genomic interaction data from Hi-C and ChIA-PET experiments.

ISoLDE - This package provides ISoLDE a new method for identifying imprinted genes. This method is dedicated to data arising from RNA sequencing technologies. The ISoLDE package implements original statistical methodology described in the publication below.

isomiRs - Characterization of miRNAs and isomiRs, clustering and differential expression.

JunctionSeq - A Utility for Detection and Visualization of Differential Exon or Splice-Junction Usage in RNA-Seq data.

kimod - This package allows to work with mixed omics data (transcriptomics, proteomics, microarray-chips, rna-seq data), introducing the following improvements: distance options (for numeric and/or categorical variables) for each of the tables, bootstrap resampling techniques on the residuals matrices for all methods, that enable perform confidence ellipses for the projection of individuals, variables and biplot methodology to project variables (gene expression) on the compromise. Since the main purpose of the package is to use these techniques to omic data analysis, it includes an example data from four different microarray platforms (i.e.,Agilent, Affymetrix HGU 95, Affymetrix HGU 133 and Affymetrix HGU 133plus 2.0) on the NCI-60 cell lines.NCI60_4arrays is a list containing the NCI-60 microarray data with only few hundreds of genes randomly selected in each platform to keep the size of the package small. The data are the same that the package omicade4 used to implement the co-inertia analysis. The references in packages follow the style of the APA-6th norm.

Linnorm - Linnorm is an R package for the analysis of RNA-seq, scRNA-seq, ChIP-seq count data or any large scale count data. Its main function is to normalize and transform these datasets for parametric tests. Examples of parametric tests include using limma for differential expression analysis or differential peak detection, or calculating Pearson correlation coefficient for gene correlation study. Linnorm can work with raw count, CPM, RPKM, FPKM and TPM. Additionally, Linnorm provides the RnaXSim function for the simulation of RNA-seq raw counts for the evaluation of differential expression analysis methods. RnaXSim can simulate RNA-seq dataset in Gamma, Log Normal, Negative Binomial or Poisson distributions.

lpsymphony - This package was derived from Rsymphony_0.1-17 from CRAN. These packages provide an R interface to SYMPHONY, an open-source linear programming solver written in C++. The main difference between this package and Rsymphony is that it includes the solver source code (SYMPHONY version 5.6), while Rsymphony expects to find header and library files on the users' system. Thus the intention of lpsymphony is to provide an easy to install interface to SYMPHONY. For Windows, precompiled DLLs are included in this package.

LymphoSeq - This R package analyzes high-throughput sequencing of T and B cell receptor complementarity determining region 3 (CDR3) sequences generated by Adaptive Biotechnologies' ImmunoSEQ assay.  Its input comes from tab-separated value (.tsv) files exported from the ImmunoSEQ analyzer.

MBttest - MBttest method was developed from beta t-test method of Baggerly et al(2003). Compared to baySeq (Hard castle and Kelly 2010), DESeq (Anders and Huber 2010) and exact test (Robinson and Smyth 2007, 2008) and the GLM of McCarthy et al(2012), MBttest is of high work efficiency,that is, it has high power, high conservativeness of FDR estimation and high stability. MBttest is suit- able to transcriptomic data, tag data, SAGE data (count data) from small samples or a few replicate libraries. It can be used to identify genes, mRNA isoforms or tags differentially expressed between two conditions.

Mergeomics - The Mergeomics pipeline serves as a flexible framework for integrating multidimensional omics-disease associations, functional genomics, canonical pathways and gene-gene interaction networks to generate mechanistic hypotheses. It includes two main parts, 1) Marker set enrichment analysis (MSEA); 2) Weighted Key Driver Analysis (wKDA).

metaCCA - metaCCA performs multivariate analysis of a single or multiple GWAS based on univariate regression coefficients. It allows multivariate representation of both phenotype and genotype. metaCCA extends the statistical technique of canonical correlation analysis to the setting where original individual-level records are not available, and employs a covariance shrinkage algorithm to achieve robustness.

MethPed - Classification of pediatric tumors into biologically defined subtypes is challenging and multifaceted approaches are needed. For this aim, we developed a diagnostic classifier based on DNA methylation profiles. We offer MethPed as an easy-to-use toolbox that allows researchers and clinical diagnosticians to test single samples as well as large cohorts for subclass prediction of pediatric brain tumors.  The current version of MethPed can classify the following tumor diagnoses/subgroups: Diffuse Intrinsic Pontine Glioma (DIPG), Ependymoma, Embryonal tumors with multilayered rosettes (ETMR), Glioblastoma (GBM), Medulloblastoma (MB) - Group 3 (MB_Gr3), Group 4 (MB_Gr3), Group WNT (MB_WNT), Group SHH (MB_SHH) and Pilocytic Astrocytoma (PiloAstro).

miRNAmeConverter - Package containing an S4 class for translating mature miRNA names to different miRBase versions, checking names for validity and detecting miRBase version of a given set of names (data from http://www.mirbase.org/).

MMDiff2 - This package detects statistically significant differences between read enrichment profiles in different ChIP-Seq samples. To take advantage of shape differences it uses Kernel methods (Maximum Mean Discrepancy, MMD).

multiClust - Whole transcriptomic profiles are useful for studying the expression levels of thousands of genes across samples. Clustering algorithms are used to identify patterns in these profiles to determine clinically relevant subgroups. Feature selection is a critical integral part of the process. Currently, there are many feature selection and clustering methods to identify the relevant genes and perform clustering of samples. However, choosing the appropriate methods is difficult as recent work demonstrates that no method is the clear winner. Hence, we present an R-package called `multiClust` that allows researchers to experiment with the choice of combination of methods for gene selection and clustering with ease. In addition, using multiClust, we present the merit of gene selection and clustering methods in the context of clinical relevance of clustering, specifically clinical outcome. Our integrative R- package contains: 1. A function to read in gene expression data and format appropriately for analysis in R. 2. Four different ways to select the number of genes a. Fixed b. Percent c. Poly d. GMM 3. Four gene ranking options that order genes based on different statistical criteria a. CV_Rank b. CV_Guided c. SD_Rank d. Poly 4. Two ways to determine the cluster number a. Fixed b. Gap Statistic 5. Two clustering algorithms a. Hierarchical clustering b. K-means clustering 6. A function to calculate average gene expression in each sample cluster 7. A function to correlate sample clusters with clinical outcome Order of Function use: 1. input_file, a function to read-in the gene expression file and assign gene probe names as the rownames. 2. number_probes, a function to determine the number of probes to select for in the gene feature selection process. 3. probe_ranking, a function to select for gene probes using one of the available gene probe ranking options. 4. number_clusters, a function to determine the number of clusters to be used to cluster genes and samples. 5. cluster_analysis, a function to perform Kmeans or Hierarchical clustering analysis of the selected gene expression data. 6. avg_probe_exp, a function to produce a matrix containing the average expression of each gene probe within each sample cluster. 7. surv_analysis, a function to produce Kaplan-Meier Survival Plots of selected gene expression data.

MultiDataSet - Implementation of the BRGE's (Bioinformatic Research Group in Epidemiology from Center for Research in Environmental Epidemiology) MultiDataSet and MethylationSet. MultiDataSet is designed for integrating multi omics data sets and MethylationSet to contain normalized methylation data. These package contains base classes for MEAL and rexposome packages.

normalize450K - Precise measurements are important for epigenome-wide studies investigating DNA methylation in whole blood samples, where effect sizes are expected to be small in magnitude. The 450K platform is often affected by batch effects and proper preprocessing is recommended. This package provides functions to read and normalize 450K '.idat' files. The normalization corrects for dye bias and biases related to signal intensity and methylation of probes using local regression. No adjustment for probe type bias is performed to avoid the trade-off of precision for accuracy of beta-values.

nucleoSim - This package can generate a synthetic map with reads covering the nucleosome regions as well as a synthetic map with forward and reverse reads emulating next-generation sequencing. The user has choice between three different distributions for the read positioning: Normal, Student and Uniform.

odseq - Performs outlier detection of sequences in a multiple sequence alignment using bootstrap of predefined distance metrics. Outlier sequences can make downstream analyses unreliable or make the alignments less accurate while they are being constructed. This package implements the OD-seq algorithm proposed by Jehl et al (doi 10.1186/s12859-015-0702-1) for aligned sequences and a variant using string kernels for unaligned sequences.

OncoScore - OncoScore is a tool to measure the association of genes to cancer based on citation frequency in biomedical literature. The score is evaluated from PubMed literature by dynamically updatable web queries.

oppar - The R implementation of mCOPA package published by Wang et al. (2012). Oppar provides methods for Cancer Outlier profile Analysis. Although initially developed to detect outlier genes in cancer studies, methods presented in oppar can be used for outlier profile analysis in general. In addition, tools are provided for gene set enrichment and pathway analysis.

PanVizGenerator - PanViz is a JavaScript based visualisation tool for functionaly annotated pangenomes. PanVizGenerator is a companion for PanViz that facilitates the necessary data preprocessing step necessary to create a working PanViz visualization. The output is fully self-contained so the recipient of the visualization does not need R or PanVizGenerator installed.

pbcmc - The pbcmc package characterizes uncertainty assessment on gene expression classifiers, a. k. a. molecular signatures, based on a permutation test. In order to achieve this goal, synthetic simulated subjects are obtained by permutations of gene labels. Then, each synthetic subject is tested against the corresponding subtype classifier to build the null distribution. Thus, classification confidence measurement can be provided for each subject, to assist physician therapy choice. At present, it is only available for PAM50 implementation in genefu package but it can easily be extend to other molecular signatures.

pcaExplorer - This package provides functionality for interactive visualization of RNA-seq datasets based on Principal Components Analysis. The methods provided allow for quick information extraction and effective data exploration. A Shiny application encapsulates the whole analysis.

PCAN - Phenotypes comparison based on a pathway consensus approach. Assess the relationship between candidate genes and a set of phenotypes based on additional genes related to the candidate (e.g. Pathways or network neighbors).

pqsfinder - The main functionality of the this package is to detect DNA sequence patterns that are likely to fold into an intramolecular G-quadruplex (G4). Unlike many other approaches, this package is able to detect sequences responsible for G4s folded from imperfect G-runs containing bulges or mismatches and as such is more sensitive than competing algorithms.

profileScoreDist - Regularization and score distributions for position count matrices.

psygenet2r - Package to retrieve data from PsyGeNET database (www.psygenet.org) and to perform comorbidity studies with PsyGeNET's and user's data.

PureCN - This package estimates tumor purity, copy number, loss of heterozygosity (LOH), and status of short nucleotide variants (SNVs). PureCN is designed for hybrid capture next generation sequencing (NGS) data, integrates well with standard somatic variant detection pipelines, and has support for tumor samples without matching normal samples.

QuaternaryProd - QuaternaryProd is an R package that performs causal reasoning on biological networks, including publicly available networks such as String-db. QuaternaryProd is a free alternative to commercial products such as Quiagen and Inginuity pathway analysis. For a given a set of differentially expressed genes, QuaternaryProd computes the significance of upstream regulators in the network by performing causal reasoning using the Quaternary Dot Product Scoring Statistic (Quaternary Statistic), Ternary Dot product Scoring Statistic (Ternary Statistic) and Fisher's exact test. The Quaternary Statistic handles signed, unsigned and ambiguous edges in the network. Ambiguity arises when the direction of causality is unknown, or when the source node (e.g., a protein) has edges with conflicting signs for the same target gene. On the other hand, the Ternary Statistic provides causal reasoning using the signed and unambiguous edges only. The Vignette provides more details on the Quaternary Statistic and illustrates an example of how to perform causal reasoning using String-db.

QUBIC - The core function of this R package is to provide the implementation of the well-cited and well-reviewed QUBIC algorithm, aiming to deliver an effective and efficient biclustering capability. This package also includes the following related functions: (i) a qualitative representation of the input gene expression data, through a well-designed discretization way considering the underlying data property, which can be directly used in other biclustering programs; (ii) visualization of identified biclusters using heatmap in support of overall expression pattern analysis; (iii) bicluster-based co-expression network elucidation and visualization, where different correlation coefficient scores between a pair of genes are provided; and (iv) a generalize output format of biclusters and corresponding network can be freely downloaded so that a user can easily do following comprehensive functional enrichment analysis (e.g. DAVID) and advanced network visualization (e.g. Cytoscape).

R4RNA - A package for RNA basepair analysis, including the visualization of basepairs as arc diagrams for easy comparison and annotation of sequence and structure.  Arc diagrams can additionally be projected onto multiple sequence alignments to assess basepair conservation and covariation, with numerical methods for computing statistics for each.

recoup - recoup calculates and plots signal profiles created from short sequence reads derived from Next Generation Sequencing technologies. The profiles provided are either sumarized curve profiles or heatmap profiles. Currently, recoup supports genomic profile plots for reads derived from ChIP-Seq and RNA-Seq experiments. The package uses ggplot2 and ComplexHeatmap graphics facilities for curve and heatmap coverage profiles respectively.

RGraph2js - Generator of web pages which display interactive network/graph visualizations with D3js, jQuery and Raphael.

RImmPort - The RImmPort package simplifies access to ImmPort data for analysis in the R environment. It provides a standards-based interface to the ImmPort study data that is in a proprietary format.

ROTS - Calculates the Reproducibility-Optimized Test Statistic (ROTS) for differential testing in omics data.

SC3 - Interactive tool for clustering and analysis of single cell RNA-Seq data.

scater - A collection of tools for doing various analyses of single-cell RNA-seq gene expression data, with a focus on quality control.

scde - The scde package implements a set of statistical methods for analyzing single-cell RNA-seq data. scde fits individual error models for single-cell RNA-seq measurements. These models can then be used for assessment of differential expression between groups of cells, as well as other types of analysis. The scde package also contains the pagoda framework which applies pathway and gene set overdispersion analysis to identify and characterize putative cell subpopulations based on transcriptional signatures. The overall approach to the differential expression analysis is detailed in the following publication: "Bayesian approach to single-cell differential expression analysis" (Kharchenko PV, Silberstein L, Scadden DT, Nature Methods, doi: 10.1038/nmeth.2967). The overall approach to subpopulation identification and characterization is detailed in the following pre-print: "Characterizing transcriptional heterogeneity through pathway and gene set overdispersion analysis" (Fan J, Salathia N, Liu R, Kaeser G, Yung Y, Herman J, Kaper F, Fan JB, Zhang K, Chun J, and Kharchenko PV, Nature Methods, doi:10.1038/nmeth.3734).

scran - This package implements a variety of low-level analyses of single-cell RNA-seq data. Methods are provided for normalization of cell-specific biases, assignment of cell cycle phase, and detection of highly variable and significantly correlated genes.

sevenbridges - R client and utilities for Seven Bridges platform API, from cancer genomics cloud to other Seven Bridges supported platforms.

SMITE - This package builds on the Epimods framework which facilitates finding weighted subnetworks ("modules") on Illumina Infinium 27k arrays using the SpinGlass algorithm, as implemented in the iGraph package. We have created a class of gene centric annotations associated with p-values and effect sizes and scores from any researchers prior statistical results to find functional modules.

SpidermiR - The aims of SpidermiR are : i) facilitate the network open-access data retrieval from GeneMania data, ii) prepare the data using the appropriate gene nomenclature, iii) integration of miRNA data in a specific network, iv) provide different standard analyses and v) allow the user to visualize the results. In more detail, the package provides multiple methods for query, prepare and download network data (GeneMania), and the integration with  validated and predicted miRNA data (mirWalk, miR2Disease,miRTar, miRandola,Pharmaco-miR,DIANA, Miranda, PicTar and TargetScan) and the use of standard analysis (igraph) and visualization methods (networkD3).

splineTimeR - This package provides functions for differential gene expression analysis of gene expression time-course data. Natural cubic spline regression models are used. Identified genes may further be used for pathway enrichment analysis and/or the reconstruction of time dependent gene regulatory association networks.

sscu - The package can calculate the selection in codon usage in bacteria species. First and most important, the package can calculate the strength of selected codon usage bias (sscu) based on Paul Sharp's method. The method take into account of background mutation rate, and focus only on codons with universal translational advantages in all bacterial species. Thus the sscu index is comparable among different species. In addition, detainled optimal codons (selected codons) information can be calculated by optimal_codons function, so the users will have a more accurate selective scheme for each codons. Furthermore, we added one more function optimal_index in the package. The function has similar mathematical formula as s index, but focus on the estimates the amount of GC-ending optimal codon for the highly expressed genes in the four and six codon boxes. The function takes into account of background mutation rate, and it is comparable with the s index. However, since the set of GC-ending optimal codons are likely to be different among different species, the index can not be compared among different species.

SwathXtend - It contains utility functions for integrating spectral libraries for SWATH and statistical data analysis for SWATH generated data.

tofsims - This packages offers a pipeline for import, processing and analysis of ToF-SIMS 2D image data. Import of Iontof and Ulvac-Phi raw or preprocessed data is supported. For rawdata, mass calibration, peak picking and peak integration exist. General funcionality includes data binning, scaling, image subsetting and visualization. A range of multivariate tools common in the ToF-SIMS community are implemented (PCA, MCR, MAF, MNF). An interface to the bioconductor image processing package EBImage offers image segmentation functionality.

transcriptR - The differences in the RNA types being sequenced have an impact on the resulting sequencing profiles. mRNA-seq data is enriched with reads derived from exons, while GRO-, nucRNA- and chrRNA-seq demonstrate a substantial broader coverage of both exonic and intronic regions. The presence of intronic reads in GRO-seq type of data makes it possible to use it to computationally identify and quantify all de novo continuous regions of transcription distributed across the genome. This type of data, however, is more challenging to interpret and less common practice compared to mRNA-seq. One of the challenges for primary transcript detection concerns the simultaneous transcription of closely spaced genes, which needs to be properly divided into individually transcribed units. The R package transcriptR combines RNA-seq data with ChIP-seq data of histone modifications that mark active Transcription Start Sites (TSSs), such as, H3K4me3 or H3K9/14Ac to overcome this challenge. The advantage of this approach over the use of, for example, gene annotations is that this approach is data driven and therefore able to deal also with novel and case specific events. Furthermore, the integration of ChIP- and RNA-seq data allows the identification all known and novel active transcription start sites within a given sample.

tximport - Imports transcript-level abundance, estimated counts and transcript lengths, and summarizes into matrices for use with downstream gene-level analysis packages. Average transcript length, weighted by sample-specific transcript abundance estimates, is provided as a matrix which can be used as an offset for different expression of gene-level counts.

Uniquorn - This packages enables users to identify cancer cell lines. Cancer cell line misidentification and cross-contamination reprents a significant challenge for cancer researchers. The identification is vital and in the frame of this package based on the locations/ loci of somatic and germline mutations/ variations. The input format is vcf/ vcf.gz and the files have to contain a single cancer cell line sample (i.e. a single member/genotype/gt column in the vcf file). The implemented method is optimized for the Next-generation whole exome and whole genome DNA-sequencing technology.

NEWS from new and existing packages
===================================

Package maintainers can add NEWS files describing changes to their
packages since the last release. The following package NEWS is available:


XXX

Packages removed since the last release
=================================

No packages were removed in this release.
17 packages were marked as deprecated, to be removed in the next release.
