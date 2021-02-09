# Artic SARS-CoV-2 Workflow

This repository contains a Nextflow workflow and associated Docker
container build for performing the Artic SARS-CoV-2 workflow on 
multiplexed MinION, GridION, and PromethION runs.

The workflow also supports using conda environments as an alternative
software isolation method to Docker.

## Quickstart

### Building the container

> This step is not necessary if you intend to run the workflow using
> conda environments.

The Docker container image can be built with the following command:

```bash
CONTAINER_TAG=epi2melabs/artic
docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=epi2melabs/base-workflow-image:latest \
    .
```

The `BASEIMAGE` argument here can be changed to use an alternative image.

### Running the workflow

The source-code repository contains two datasets to test the workflow:
a) a source set of files containing reads with a mixture of barcodes and,
b) a pre-demultiplexed set of .fastq files (derived from the source set).

**Running the workflow with Docker containers**

To run the workflow using Docker containers supply the `-profile standard`
argument to `nextflow run`:

```
OUTPUT=workflow-test
nextflow run workflow.nf \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --fastq test_data/sars-samples-mixed/  \
    --out_dir ${OUTPUT}
```

The output of the pipeline will be found in `./workflow-test` for the above
example. This directory contains the nextflow working directories alongside
the two primary outputs of the pipeline.

**Using conda environments**

To run the workflow backed by conda environments, simply provide the
`-profile conda` argument to `nextflow run`.

```
OUTPUT=workflow-test
nextflow run workflow.nf \
    -w ${OUTPUT}/workspace \
    -profile conda \
    --fastq test_data/sars-samples-mixed/  \
    --out_dir ${OUTPUT}
```

This will create a conda environment with all required software within the
workspace directory. When running multiple analyses on distinct datasets
it may not be desirable to have Nextflow create a conda environment for each
analysis. To avoid the situation editing the file `nextflow.config` will
be necessary. Search for the term `cacheDir` and set this to a directory
where you wish the conda environment to be placed.
