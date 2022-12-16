# The Protein Analysis pipeline using Nextflow
Pipeline created by [Bio2Byte](https://bio2byte.be) research group intended to teach how to use Nextflow.

![pipeline](https://github.com/agdiaz/sequencesAnalysis/blob/main/pipeline.png)

## Usage

### Local environment

`nextflow run agdiaz/sequencesAnalysis -r main --targetSequences /path/to/file.fasta -profile standard,withsingularity`

### VUB-HPC (Hydra) context

`nextflow run agdiaz/sequencesAnalysis -r main --targetSequences /path/to/file.fasta -profile hydra,withsingularity`
