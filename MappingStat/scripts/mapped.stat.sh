#!/usr/bin/env bash

# @Author: jsgounot
# @Date:   2020-12-09 11:18:39
# @Last Modified by:   jsgounot
# @Last Modified time: 2025-03-11 15:06:07

mapped=$(samtools view -c -F 4 "${snakemake_input[0]}")
unmapped=$(samtools view -c -f 4 "${snakemake_input[0]}")
sumvalue=$(($mapped+$unmapped))
prc=$(bc <<< "scale=4;$mapped/$sumvalue")
echo "mapped:$mapped" > "${snakemake_output[0]}"
echo "unmapped:$unmapped" >> "${snakemake_output[0]}"
echo "sum:$sumvalue" >> "${snakemake_output[0]}"
echo "prc:$prc" >> "${snakemake_output[0]}"
