The Duo paternity project
================

### Percentage of exclusion per locus

| marker   |   perc|
|:---------|------:|
| D12S391  |  26.46|
| D2S1338  |  25.86|
| D1S1656  |  25.57|
| PentaE   |  25.48|
| D18S51   |  24.74|
| FGA      |  24.35|
| D21S11   |  22.78|
| PentaD   |  22.61|
| D19S433  |  22.16|
| D8S1179  |  20.62|
| vWA      |  20.16|
| D7S820   |  19.83|
| D13S317  |  19.32|
| D16S539  |  19.31|
| TH01     |  19.21|
| D2S441   |  18.92|
| D10S1248 |  18.89|
| D3S1358  |  18.36|
| D22S1045 |  16.88|
| CSF1PO   |  16.45|
| D5S818   |  16.42|
| TPOX     |  14.93|

### Trios sets

We rearranged the trios in different sets of loci, giving the following sample sizes:

| marker\_set |   total|  exclusion|
|:------------|-------:|----------:|
| codis       |  38,320|     13,105|
| identifiler |  39,211|     13,352|
| pp16        |  78,459|     24,797|

After reducing the original number of loci (up to 22) to 18 or 15, the paternity status change in a few cases. Below we see the number of cases which changed from the status in column `before` to the status in column `after`.

| before    | after        | marker\_set |    n|
|:----------|:-------------|:------------|----:|
| exclusion | inclusion    | codis       |    0|
| exclusion | inclusion    | identifiler |    1|
| exclusion | inclusion    | pp16        |    0|
| exclusion | inconclusive | codis       |    2|
| exclusion | inconclusive | identifiler |    5|
| exclusion | inconclusive | pp16        |    8|
| inclusion | inconclusive | codis       |   12|
| inclusion | inconclusive | identifiler |  326|
| inclusion | inconclusive | pp16        |  344|

Definitions:

-   Exclusion: mismatch between AF and child at ≥3 loci.

-   Inclusion: Less than 3 mismatches and PI ≥ 10,000.

-   Inconclusive: Less than 3 mismatches and PI &lt; 10,000.

### False inclusion in Duos

Then, we sought to evaluate the rate of false inclusions in Duo cases. We took the Trios rearranged in the Codis, Identifiler and PP16 marker sets, we removed the mother, and we recalculated to probability of paternity.

Given the mutation models below, we observe the following false inclusions:

#### LR = 0.001 in case of mismatch between AF and child

| marker\_set |  case\_no| trio        |  exclusions|        cpi|
|:------------|---------:|:------------|-----------:|----------:|
| codis       |    250992| M1\_F1\_SP1 |           1|  144542.54|
| identifiler |    253626| M1\_F1\_SP1 |           0|   38484.42|
| pp16        |    117770| M1\_F1\_SP1 |           1|   10892.31|
| pp16        |    155180| M1\_F1\_SP1 |           1|   10661.89|

#### Step-wise mutation model, mutation rates from STRbase, mutate range = 0.1, point mutation rate = 0.00001

| marker\_set |  case\_no| trio        |       cpi|
|:------------|---------:|:------------|---------:|
| identifiler |    253626| M1\_F1\_SP1 |  38113.29|
| pp16        |    117770| M1\_F1\_SP1 |  13829.57|
