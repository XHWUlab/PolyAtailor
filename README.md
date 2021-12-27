# PolyAtailor
Multifunctional R package for poly(A) tail analysis.
## About
Polyatailor starts with the original sequencing data, first pre-processed the original data of different sequencing technology, and converted into FASTQ format files. Then enter the TAIL_SCAN process, perform tail extraction, tail filtering, and tail classification to get TAIL_SCAN TAILS. If the user provides a reference genome, PolyAtailor will continue to use TAIL_MAP to fix the preliminary results of Tail_scan to get TAIL_MAP TAILS (Figure below). After obtaining accurate tail data, the user can perform the recognition, annotation of the PA bit, and the visualization analysis of the Poly (A) tail base. In addition, users can also analyze significant analysis of PLs of different conditions via Polyatailor.   
<img src="./overview.png" width="100%" />
