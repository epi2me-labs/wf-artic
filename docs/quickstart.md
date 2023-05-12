## Quickstart

The workflow uses [Nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such Nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-artic --help
```

to see the options for the workflow.

**Workflow outputs**

The primary outputs of the workflow include:

* a [FASTA](https://en.wikipedia.org/wiki/FASTA) file containing the consensus sequence for all samples,
* a [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file sample all samples,
* an HTML report document detailing QC metrics and the primary findings of the workflow.
