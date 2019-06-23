#!/bin/bash

awk 'FNR > 1 || NR == 1' sharetable_{1..100}.tsv > sharetable_all.tsv

rm sharetable_{1..100}.tsv

awk '$3 > 0.9' sharetable_all.tsv > dupstable.tsv
