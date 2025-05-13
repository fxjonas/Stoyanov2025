This Repository contains the scripts and data to reproduce the analysis and figures presented in Stoyanov 2025.
borisAna.m: Code to create all the analysis
borisColl-1: Parsed data from AID study and ChEC-seq profiles of the Barkai lab
borisColl-2: Parsed data from Deleteome and IDEA study and ChEC-seq profiles of the Mahendrawada et al. (2025)
-- both these data sets are parsed so that the columns correpsond to a specifc TF (can be seen as Column GTitle) and each row (6701) corresponds to a specifc gene/promoter. THe order of promoters can be found in GP.gene_infoR64.nameNew

AFoutDatabig: Alphafold analysis for TF-Co-factor interactions
AFoutDatasmall: Alphafold analyssi for TF-CoFActor interactions (with less cofcactors)
allChECdata: more ChECseq data in the same format (6701 rows)
group_impGitHub.mat: information on the genes according to the general order (6701 rows)
TFstats.mat: summary of the analysis for each TFs (it is created durign borsAna.m but cna be loaded directly for further analysis)
deAnalysis: Collected meta data on TF's from other sources
TFdesc: TF descriptions from SGD and correpsonding classificaiton into Activtors (1) and repressors (-1)
AF-Cyc8.pdb : Aphafold model for CYC8
bestPBMonScer: Motif occurences of all TF motifs in the yeast genome
