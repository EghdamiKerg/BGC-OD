# BGC-OD
BGC-OD: Block-Boosted GC-aware Overlap Detector

Input:

FASTQ format for input reads

Output:

Overlap candidates in CAN format
Pairwise alignments in M4 format (in future)

Options:

-j [task] (0 or 1): Select between detecting overlapping candidates only (0) or outputting pairwise alignments in M4 format (1). Pairwise alignment will be completed in future. 
-d [input] : Input file.
-w [working folder] : Directory for temporary files.
-t [# of threads] : Number of CPU threads.
-o [output] : Output file name.
-n [# of candidates] : Number of candidates to consider for gapped extension.
-g [0/1] : Whether to output gapped extension start points.
-x [0/1] : Sequencing platform (0 for PacBio, 1 for Nanopore).

Note: The -n parameter is crucial for optimizing FODI's performance. It's recommended to adjust this parameter based on genome size and read length.
