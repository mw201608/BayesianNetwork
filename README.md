## Bayesian network analysis pipeline

This repository provides a simple example pipeline to perform Bayesian probabilistic network inference using the [RIMBANet](https://labs.icahn.mssm.edu/zhulab/?s=rimbanet) software.

Usage:

1. Download folders `bin`, `script`, and `Perls`. Place them in a directory say `$HOME/user/RIMBANet/`. Add eXecute permission to files in the two folders. Lastly, either add `$HOME/user/RIMBANet/bin` to your system `PATH` variable or move `bin/testBN` to a place included in `PATH`.

2. Prepare discretized data matrix (discretized continuous values into odrinal `0`, `1`, and `2`) and save in a tab-delimited text file named `data.discretized.txt` in your working directory, say `/path/to/working/directory`. In this file, each row represents a variable, and each column represents a sample except the first columns which consists of variable names. However, there should be no header row.

3. Call RIMBANet to construct networks with different random starting structures. A wraper script `caller.bsh` is provided to create all required files from the data and then submit parallel jobs on a LSF environment. To change any of the default behabivors, please edit `caller.bsh`, `BNJob.bsh`, and `BNJobbase.bsh`. To add additional prior information, please refer to Perl scripts in folder `Perls`.

```
bash $HOME/user/RIMBANet/script/caller.bsh /path/to/working/directory
```

4. When LSF jobs are completed successfully, create a consensus network structure (default final output `result.links3.links.txt`) consists of de-looped links present in, by default, at least 30% of the individual structure.
```
bash $HOME/user/RIMBANet/script/createStructure.bsh bn.param.txt
```
