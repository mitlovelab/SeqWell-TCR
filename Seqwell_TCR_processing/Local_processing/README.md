# CDR3 processing of TCR sequencing results
Here contains the scripts necessary to analyze the raw TCR sequencing results. 

## Before starting
The goal is to produce a summary file of the sequencing results that contain the CDR3 call
associated with each cellular barcode and UMI, as well as quality of the call.

To speed up processing, it is highly recommended that a list of cellular barcodes be generated
from the single-cell transcriptomic data. Providing this list of barcodes will limit the analysis
to just the cells with those barcodes.