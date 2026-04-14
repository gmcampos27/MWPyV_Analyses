# MWPyV_Analyses
Pipeline designed for the genomic characterization of __MW Polyomavirus (MWPyV)__. This workflow automates from raw NGS data preprocessing to consensus sequence generation and quality metrics reporting.


Quality Control: Raw data processing using __fastp v0.23.4__, retaining bases with Q > 20 and discarding reads < 15 bp. 

Automated Reference Selection: Reads are queried against a custom database of 24 _Deltapolyomavirus decihominis_ (MWPyV) genomes using __VAPOR v1.0.2__ to identify the most biologically similar reference.

Read Mapping: Alignment with __BWA-MEM v0.7.17__ followed by processing with __SAMtools v1.19__.

Genome Polishing: Correction of gaps and small indels using __Pilon v1.24__ (parameters: depth $\ge 5\times$, base quality $\ge 20$, mapping quality $\ge 10$).

Consensus Generation: Final consensus sequences are called using __iVar v1.4.2__ (thresholds: $Q20$, $5\times$ depth). Regions below the depth threshold are masked with "N".

Phylogenetic Readiness: Only genomes with $>80\%$ coverage breadth are recommended for downstream analysis (as calculated in the final metrics report).

Beyond the assembly pipeline, this study utilized a phylogenetic workflow for both whole-genome and VP1 gene sequences:

Phylogenetic trees were reconstructed for three different datasets (Complete Genome, VP1 Gene, and Complete Genome alongside JCPyV and BKPyV). 

**Methodology:**
- **Alignment:** Sequences were aligned using `MAFFT v7.520`.
- **Tree Construction:** Maximum-Likelihood (ML) trees were inferred using `IQ-TREE2`.
- **Command:**
  ```bash
  iqtree2 -s input_alignment.fasta -fast -m MFP -alrt 1000 -nt AUTO

Visualization & Annotation: * Initial visualization via FigTree v1.4.4.

Advanced graphical editing and annotation using the __ggtree package__ (v3.16.3) in R v4.5.0.
