# ARTIC SARS-CoV-2 Workflow

This repository contains a Nextflow workflow and associated Docker
container build for running the ARTIC SARS-CoV-2 workflow on 
multiplexed MinION, GridION, and PromethION runs.

The workflow also supports using conda environments as an alternative
software isolation method to Docker.

## Quickstart

`nextflow run epi2me-labs/wf-artic --help` 

SHOULD WE INCLUDE AN EXAMPLE WITH SOME TOY DATA HERE?

`nextflow run epi2me-labs/wf-artic ...`

## Supported installations and GridION devices

Installation of the software on a GridION can be performed using the command

`sudo apt install ont-nextflow`

This will install a java runtime, Nextflow and docker. If *docker* has not already been
configured the command below can be used to provide user access to the *docker*
services. Please logout of your computer after this command has been typed.

`sudo usermod -aG docker $USER`


## Installation and dependencies

> Nextflow requires a java runtime (JRE)

`sudo apt install default-jre`

> Nextflow may be downloaded from https://www.nextflow.io or through conda

`curl -s https://get.nextflow.io | bash`

THIS PLACES THE NEXTFLOW BINARY IN CWD ... SHOULD WE SAY ANY MORE HERE?

> Docker is recommended 

```
sudo apt install docker.io
sudo usermod -aG docker $USER
```

## Running the workflow

The `wf-artic` workflow can be controlled by the following parameters. The `fastq` parameter
is the most important parameters because it is required to identify the location of the
sequence files to be analysed. It is also important to note the `scheme_version` parameter
that should be changed from the default **V3** value if e.g. the LoCost, ECO or midnight 
variants of the ARTIC protocol are being used.

Parameters:
- `fastq` specifies a *directory* path to FASTQ files (required)
- `samples` locates a CSV file with columns named `barcode` and `sample_name`
   (or simply a sample name for non-multiplexed data) - this is used to replace
   barcode identifiers with other sample identifiers
- `out_dir` the path for output (default: output)
- `medaka_model` the medaka model name (default: r941_min_high_g360) to use during
   the consensus sequence polishing.
- `min_len` Minimum read length (default: set by scheme)
- `max_len` Maximum read length (default: set by scheme)
- `scheme_version` Primer scheme ([V1, V2, V3, V1200]

## Running the workflow with Conda

I HAVEN'T TESTED WORKFLOW WITH CONDA - PLEASE CAN WE HAVE A COUPLE OF STATEMENTS HERE TO
REASSURE?


## Updating the workflow

`nextflow pull epi2me-labs/wf-artic`

## Configuration and tuning

The default settings for the workflow are described in the configuration file `nextflow.config`.
This file defines an *executor* that can use a maximum of four CPU cores and eight gigabytes of
RAM. If the workflow is being run on a device other than a GridION, the available memory and
number of CPUs may be adjusted to the available number of CPU cores.

## Building the docker container from source

```
CONTAINER_TAG=ontresearch/wf-artic

git clone https://github.com/epi2me-labs/wf-artic
cd wf-artic

docker build \
    -t ${CONTAINER_TAG} -f Dockerfile \
    --build-arg BASEIMAGE=ontresearch/base-workflow-image:v0.1.0 \
    .
```

