#!/bin/bash

awk 'FNR > 1 || NR == 1' /scratch/vitor/str/sharetable_{1..1000}.tsv > sharetable_all.tsv

rm /scratch/vitor/str/sharetable_{1..1000}.tsv
