## TE-pipeline

TE-pipeline is an integrated and automated workflow for genome-wide transposable element (TE) annotation in plant genomes.  
The pipeline combines **EDTA**, **RepeatModeler**, **RepeatMasker**, **DeepTE**, and **TEsorter** to generate a high-confidence, well-classified TE annotation set, suitable for downstream genome annotation, methylation analysis, and comparative genomics.

This workflow is designed for **HPC environments** and supports large, complex genomes.

## Overview

The pipeline performs TE annotation in multiple complementary stages:

- Structure-based and homology-based TE annotation using **EDTA**
- De novo repeat discovery using **RepeatModeler**
- Genome-wide repeat masking using **RepeatMasker**
- Unknown TE re-annotation using **DeepTE** and **TEclassify**
- Protein-domain–based refinement using **TEsorter**
- Standardization of TE classification and final GFF generation

The goal is to maximize sensitivity for TE discovery while improving classification accuracy, especially for previously unclassified elements.

## Workflow

The TE-annotation-pipeline consists of four major stages:

### 1. Initial TE annotation
- EDTA (structure-based + curated library)
- RepeatModeler (de novo library construction)
- RepeatMasker (homology-based annotation)

### 2. Library integration and genome masking
- Merge EDTA and RepeatModeler libraries
- Generate soft-masked and hard-masked genomes
- Produce combined GFF and BED annotations

### 3. Re-annotation of unknown TEs
- Extract unknown TE sequences
- DeepTE classification
- TEclassify-based genome-wide statistics

### 4. Protein-domain validation and refinement
- TEsorter domain-based classification
- Resolve conflicts and standardize TE nomenclature

## Installation

### Required software

The pipeline relies on the following tools:

- EDTA
- RepeatModeler
- RepeatMasker
- DeepTE
- TEsorter
- BEDTools
- Perl (with standard bioinformatics modules)
- Python (≥ 3.7)

It is recommended to manage dependencies using **conda environments**.

### Conda environments
All required conda environments are provided in the `Installation/` directory as YAML files.

```bash
conda env create -f Installation/TEanno1.yaml
conda env create -f Installation/TEanno2.yaml
conda env create -f Installation/TEanno3.yaml
```
Activate environments as needed during execution.

```bash
conda activate TEanno1
conda activate TEanno2
conda activate TEanno3
```
- TEanno1 → EDTA + RepeatMasker environment
- TEanno2 → DeepTE environment
- TEanno3 → TEsorter environment

## Usage on HPC

The pipeline is designed for PBS-based HPC systems.

### Job submission header

```bash
#PBS -j oe
#PBS -q queue_name
#PBS -V
#PBS -l nodes=1:ppn=40
```

### Input parameters

Edit the main shell script before submission:

```bash
fa=/path/to/genome.fasta
genome_size=334768680
name=sample.hapA
rp_species=Carica_papaya
cpu=40
```

### Parameter description:

- fa → genome FASTA file
- genome_size → total genome size (bp), required for TEclassify
- name → sample prefix
- rp_species → RepeatMasker species parameter
- cpu → number of CPU threads

### Running the pipeline

Submit the job using:

```bash
qsub TE-pipeline.sh
```
The pipeline runs sequentially and automatically switches conda environments as required.

### Directory structure

Typical output structure:
```md
.
├── 00.EDTA/
│   ├── genome.fasta.mod.EDTA.raw/
│   └── EDTA annotation results
│
├── 01.rpanno/
│   ├── rmout_sample/
│   ├── rmout_denovo_sample/
│   └── rmout_combine_sample/
│       ├── sample.out
│       ├── sample.out.gff
│       ├── sample.sm.fa
│       ├── sample.rm.fa
│       └── Unknown_TE/
│           ├── DeepTE results
│           ├── TEclassify results
│           └── results.txt
│
├── TEsorter/
│   ├── *.dom.gff3
│   └── TEsorter.log
│
└── formatted.gff
```

## Output files

Key output files include:

### GFF files
- `formatted.gff` → final standardized TE annotation  
- `update.out.gff` → DeepTE-updated annotations  

### Masked genomes
- `*.sm.fa` → soft-masked genome  
- `*.rm.fa` → hard-masked genome  

### Statistics
- `results.txt` → TEclassify genome-wide TE composition  

### Intermediate files
- RepeatMasker output files (`*.out`, `*.gff`, `*.bed`)  
- Protein-domain–based annotations from **TEsorter**

## TE classification standard

The pipeline standardizes TE names into a unified hierarchical format, for example:

- `Class_I/LTR/Ty3_gypsy`
- `Class_I/LTR/Ty1_copia`
- `Class_I/LINE`
- `Class_I/SINE`
- `Class_II/subclass_1/TIR/hAT`
- `Class_II/subclass_1/TIR/MuDR/Mutator`
- `Class_II/subclass_2/Helitron`

Non-TE annotations (e.g. tRNA, rRNA, simple repeats) are filtered out from the final results.

## Notes

### RepeatMasker and RepeatModeler configuration

RepeatMasker and RepeatModeler are **not fully configured by conda alone**.  
After installation, users must manually configure these tools following their official documentation, including:

- Setting the correct paths to alignment engines (e.g. RMBlast)
- Configuring RepeatMasker with required dependencies
- Verifying RepeatModeler database and executable paths

Please refer to the official manuals for detailed setup instructions:
- RepeatMasker: https://github.com/Dfam-consortium/RepeatMasker
- RepeatModeler: https://github.com/Dfam-consortium/RepeatModeler

### Dfam database

RepeatMasker requires a properly configured **Dfam database**, which must be downloaded and indexed manually.

The pipeline was tested using **Dfam release 3.9**.  
Users should download the following files and configure them according to RepeatMasker guidelines:

- https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.0.h5.gz  
- https://www.dfam.org/releases/Dfam_3.9/families/FamDB/dfam39_full.5.h5.gz  

After downloading, make sure the Dfam database is correctly linked and recognized by RepeatMasker.

### MIPS PlantDB (curated repeat library)

This pipeline uses the **MIPS PlantDB repeat library** as a curated TE reference database for EDTA.

MIPSPlantsDB is a comprehensive plant genome resource for integrative and comparative genomics:

Spannagl M, Noubibou O, Haase D, et al.  
*MIPSPlantsDB – plant database resource for integrative and comparative plant genome research.*  
Nucleic Acids Research, 2007, 35(Database issue).  
PMID: 17202173

Users should obtain the MIPS PlantDB library independently and configure the path in the EDTA step:

```bash
--curatedlib /path/to/MIPS_PlantDB_library.fa
```

## Citation
If you use this pipeline in your research, please cite the following tools and their respective publications:

1. **EDTA (Extensive de novo TE Annotator)**  
   Ou S, Chen J, Jiang N, et al. EDTA: automated de novo annotation of transposable elements in eukaryotic genomes. *Genome Biology* (2019).  
   DOI: [10.1186/s13059-019-1921-0](https://doi.org/10.1186/s13059-019-1921-0)

2. **RepeatModeler2**  
   Flynn JA, Hubley R, Goubert C, et al. RepeatModeler2 for automated genomic discovery of transposable element families. *PNAS* (2020).  
   DOI: [10.1073/pnas.1921046117](https://doi.org/10.1073/pnas.1921046117)

3. **RepeatMasker**  
   Smit AFA, Hubley R, Green P. RepeatMasker Open-4.0. *2013-2015*.  
   Available at: [http://www.repeatmasker.org](http://www.repeatmasker.org)

4. **DeepTE**  
   Yan H, Wang J, Li Y, et al. DeepTE: a deep learning-based tool for transposable element annotation. *Briefings in Bioinformatics* (2021).  
   DOI: [10.1093/bib/bbab404](https://doi.org/10.1093/bib/bbab404)

5. **TEsorter**  
   Zhang H, Wang X, Wang J, et al. TEsorter: an efficient and flexible tool for classifying transposable elements based on protein domains. *Bioinformatics* (2021).  
   DOI: [10.1093/bioinformatics/btab708](https://doi.org/10.1093/bioinformatics/btab708)
