# Introduction to implementation

This file serves to briefly explain which concepts can be called up 
and how. It also explains what the output files stand for.

## Pipeline

### Format of File containing reads

The program input is modelled according to the file 
5V_allPairs.txt, 5VirusMixReference.fasta and 5V_ref_seq.fasta. 
These file names were hard coded in the implementation.  

### Correction of non-variant call positions

In a first step, the reads were corrected using Shannon 
Entropy. To achieve this, the main function of the class 
entropyErrorCorrection can be executed. The reads are 
then saved in a file called Documents\\EntropyFile.txt. 
The positions that are relevant respectively have not been corrected 
are saved in the file Documents\\filteredPositions.txt. 
Documents\\checkCorrection.npy lists which reads were changed 
during the correction. These must not be used in the 
scoring system.

Important methods:

entropyCorrection: applies the correction


### Correction using clustering approach

To correct those positions which were not considered before, 
the clustering concept was used. It requires the EntropyFile.txt from before. 
To execute it, the main function of the clusteringErrorCorrection class can be 
used. Using threads, it tries to parallelize the correction. For each window 
created, a file is saved in the Data folder. A file lists for a read to which 
haplotype it has been assigned the most. After all windows have been created, 
the reads are corrected using the majority rule. They are then saved in the file 
Documents\\ClusteringFile.txt. In our version, the correction was done 3 times: 
with a window of lengths 500, 600 and 700. After the correction with the 500-windows, 
the ClusteringFile.txt was manually renamed to EntropyFile.txt. The correction was 
then repeated on the new EntropyFile with windows of length 600. The same was done 
for the 700-windows. To change the window length, the variable windowLength must be 
changed in the pipe() method. 

Important methods:

pipe: executes all relevant steps for a window with a certain length

sampleHaplotype: samples the nucleotides of a haplotype with its read set

assignToClass: selects a class for a read.

updateGamma: updates gamma

updateTheta: updates theta

### Generation of the global haplotypes
After the reads have been corrected and the final ClusteringFile has been created, 
the global haplotypes are constructed. To do this, the main function of the class 
haplotypeReconstruction_DBG can be executed. It generates and uses several important 
data structures. The reads from the ClusteringFile are stored in the clusteringTable. 
Since new non-variant call positions may be corrected during error correction, the 
relevant positions must be recalculated. They are stored in variantPositions and  
the file Documents\\variantPositions.txt. The reads are then transformed into the compact 
form and stored in the file Documents\\PairedEndReadList.txt. The split reads are saved 
in Documents\\ReadList.txt. In addition, we also calculate the compact form of the raw 
reads. They are listed in Documents\\filteredReadTable. They will be needed for the 
scoring system. In the file Documents\\tracker.txt it is noted for each paired-end read 
which are its predecessor raw reads. This is relevant because they will later be important 
for the scoring system.

After these data structures have been generated, the algorithm for generating 
the global haplotypes is executed. When generating a read graph, the positions where branches 
exist are noted in parallel. All these positions are stored in the file Documents\\checkPositions.txt. 
They may be relevant for the consistency check and the scoring system.  For each subregion, 
a file is created in the Haplotype folder. It contains all its local haplotypes. After the global 
haplotypes have been calculated, they are stored in the file Haplotype\\GlobalHaplotypes.txt. 

Important methods:

pipe: executes the entire program.

getLocalHaplotypes: controls the computation of the read graph in a given region. Here, k from the 
k-mer structure can be defined.

removeInconsistentHaplotypes: removes possible local/global haplotypes by checking the flow of 
extending reads. It is also used in the consistency check.

deBruijnAlgorithm: computes for a given region its read graph.

### Optimization

To find the optimal set of global haplotypes, the globalHaplotypeReconstruction class can 
be used. The files mentioned above must exist. To run the algorithm, the main 
function can be run. It is important to specify the border and thus the interval of the global 
region in the compact form in the method pipe(). The pipe calculates the different h#. A file is 
created for each h# in the Frequencies folder. Then the optimal result is calculated and saved 
in Solutions\\Solution.txt. In mode 1 (which can be defined in pipe) the consistency check method is 
run first to reduce the number of global haplotypes. The class also includes the implementation 
of the scoring system. To do this you need the solution file.

Important methods:

pipe: executes the pipeline

calculateBestFit: computes the best h#

consistencyCheck: removes false haplotypes before optimization

defineScore: scoring system

### Helper classes
StageWiseRegression: computes the regression

LagrangianRegression: computes the regression

Experiment: computes the variants for the experiments

HelperClasses: contains datastructures relevant for the main implementations

Statistic: contains tests to examine the error correction
