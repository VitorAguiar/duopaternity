#!/bin/bash

awk 'FNR > 1 || NR == 1' sharetable_{1..100}.tsv > sharetable_all.tsv

rm sharetable_{1..100}.tsv
