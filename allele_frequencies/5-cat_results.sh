#!/bin/bash

awk 'FNR > 1 || NR == 1' /scratch/vitor/str/matchtbl_{1..100}.tsv > matchtbl_all.tsv

#rm /scratch/vitor/str/sharetbl_{1..100}.tsv
