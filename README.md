The Duo paternity project
================

From 128,140 Trios, we calculated the percentage of exclusion per locus:

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

Then, we rearranged the trios in different sets of loci, giving the following sample sizes:

| marker\_set |  total|  exclusion|
|:------------|------:|----------:|
| codis       |  38320|      13105|
| identifiler |  39211|      13352|
| pp16        |  78459|      24797|

After reducing the original number of loci (up to 22) to 18 or 15, we observe that the paternity status would change in a few cases. Below we show the number of cases which change from the status in column `before` to the status in column `after`.

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

Then, we sought to evaluate the rate of false inclusions in Duo cases. We took the Trios rearranged in the Codis, Identifiler and PP16 marker sets, we removed the mother, and recalculated to probability of paternity.
