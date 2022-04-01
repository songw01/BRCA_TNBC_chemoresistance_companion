### BRCA TNBC chemoresistance companion codes ###

*Correspondence: Won-Min Song (won-min.song@mssm.edu or wonmin1984@gmail.com)*

### 
This Github page houses the companion R codes to reproduce results and figures for the study. The input data can be downloaded at https://www.synapse.org/#!Synapse:syn28711096, or by using the following R codes: 

library(synapser) 
library(synapserutils) 
 
synLogin('synapse_username', 'password') 
files <- synapserutils::syncFromSynapse('syn28707254')

# Code descriptions
Upon the download, the data must be unzipped under folder named *"Data".*

Then, the R codes in '/scripts' folder, should be executed in the following order: 

trajectory_inference_slingshot.R - runs trajectory analysis single-cell cluster marker calculations. Reproduces Figure 2. 

bulk_GSVA_analysis.R - runs GSVA analysis with single-cell cluster markers, and identify correlations to the survival. Reproduces Figure 3. 

generate_consensus_signature.R - analyzes multiple TNBC chemo-resistance gene signatures to derive consensus signature. Reproduces Figure 4. 

combo_signature_analysis.R - analyzes RNA-seq of drug treated MDA-MB-231 cells to reproduce Figure 5. 



