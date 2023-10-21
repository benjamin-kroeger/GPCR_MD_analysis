# Notes

## Calculating volume overlap

- Both tools proposed by antonella are incapable of calculating overlap
- alternative python packages:
    - Biopython PDB module (not specific enough computes same vol for different strutures)
    - RDKIT
    - VOlumina 
    - openeye (also needs to be licenced)
    - RDKit shapeit
    - USRCAT http://www.jcheminf.com/content/4/1/27/abstract https://bitbucket.org/aschreyer/usrcat
    - ROCS (openeye)
    - UGM http://nbviewer.ipython.org/4316435/
    - CASP analysis pipieline

# TODO

- determine ideal number of clusters
- reach out to developers of volume tool byreuth
- Color alignment of all sequences
- plot sequence identity vs RMSD eg
- dendogramm on rmsd tree
- denodogramm on seuqence identiy
- SDP analysis
- Specifiicty determininbg position residues ( do in on sequences grouped by protein functions)
- f'/home/benjaminkroeger/Documents/ResearchRodedGroup/MG_shap_analysis/Runs/{args.experiment}


ructural Alignment Algorithms: Use structural alignment algorithms designed for comparing 3D structures without considering sequence identity. Some popular choices include:

MAMMOTH: A structural alignment program that can handle structurally diverse proteins.
TM-align: A tool for structural alignment based on the TM-score.
DALI: A program for pairwise protein structure comparison.
# Dimitri Remarks

Dear Antonella and Benjamin,

just to summarize my input briefly: the main question is whether structure-function relationships within EL2 make sense, i.e. whether or not there is any kind of concordance between functional annotation/sequence evolution and structural similarity. What I would do if I were doing this project:

- produce and (if necessary) manually curate a definitive multiple sequence alignment of EL2 regions. Alternatively, check out the PFAM database - perhaps there is already a good alignment there
- calculate all-on-all sequence similarity both for full length receptor sequences and EL2 domains only
- calculate all-on-all structural similarity between EL2 loops
- calculate sequence- and structure-based trees and compare them
- calculate the %identity vs RMSD plot as we discussed
- consider available functional annotation of the receptors from public databases (Uniprot, GPCRdb, …)
- split the alignment into functional groups and try to identify specificity determining residues in EL2, if any

A much larger question is whether or not - and how - to calculate the loop structure in bound state.

Just my two cents - of course Antonella is the boss!

Cheers,
Dmitrij