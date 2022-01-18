#!/usr/bin/env bash

# Get nextclade data from their github repository

data_dir_location=$1

cd ${data_dir_location}

curl https://codeload.github.com/nextstrain/nextclade_data/tar.gz/master | tar -xz --strip=2 nextclade_data-master/data/datasets/sars-cov-2
