The Duo paternity project
================

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

We compared the observed rate of false inclusion in the real data to that observed in 1 million simulated Duo cases.

Given the mutation models below, we observe the following false inclusions:

| marker\_set | model   |       obs|  simulation|
|:------------|:--------|---------:|-----------:|
| codis       | r001    |  7.63e-05|     1.0e-06|
| identifiler | r001    |  7.49e-05|     1.1e-05|
| pp16        | r001    |  8.07e-05|     9.0e-06|
| identifiler | STRbase |  7.49e-05|     1.3e-05|
| pp16        | STRbase |  4.03e-05|     9.0e-06|

Mutation models:

-   r001: LR = 0.001 is assigned in case of mismatch between AF and child

-   STRbase: Step-wise mutation model, mutation rates from STRbase, mutate range = 0.1, point mutation rate = 0.00001
