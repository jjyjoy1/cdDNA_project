Cell-free DNA (cfDNA) and cell-free RNA (cfRNA) detection are widely used in non-invasive diagnostics. cfDNA is primarily utilized in liquid biopsies for cancer detection, minimal residual disease (MRD) monitoring, prenatal testing, and transplant rejection assessment. cfRNA, being more dynamic, provides insights into gene expression changes and is used for disease biomarker discovery, immune response profiling, and monitoring treatment effects in real time. The combination of cfDNA and cfRNA enhances the sensitivity and specificity of detecting disease-related molecular signatures in blood and other biofluids.

In this project, I try to integrate cfDNA and cfRNA detection together. 

Assay Focus Comparison
Rare Variant Analysis:
Objective:
Detect extremely frequency sequence variants.
Advantages:
Ultra-sensitive detection enabled by high-depth sequencing and UMI-based error correction.
Useful for identifying low-level signals amidst high background noise or contamination.
Challenges:
Requires rigorous error-correction and consensus building because false positives from PCR or sequencing errors are a significant concern in cfDNA.
High sequencing depth is needed, which can be costly, and the approach may be more prone to false positives if not carefully controlled.

Reads Count Analysis:
Objective:
Quantify the gene expression read when low amount material available, low microbial DNA, such as rare pathogens or low abundance of different bacterial taxa based on 16S read counts, akin to standard microbiome profiling.
Advantages:
Suitable for profiling the gene expression, and overall microbial community or microbiome in blood.
Less demanding in terms of sequencing depth compared to rare variant detection.
For 16S target sequecing, often leverages established bioinformatics pipelines (like DADA2, QIIME, or mothur) to generate operational taxonomic units (OTUs) or amplicon sequence variants (ASVs).
Challenges:
With limited starting material, obtaining robust and representative counts can be difficult.
Environmental and reagent contamination is a common issue in low-biomass samples like blood, so distinguishing true signals from noise requires stringent controls and normalization.
