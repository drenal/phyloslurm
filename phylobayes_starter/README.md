# Phylobayes(-MPI) starter and governor script

The scripts in this folder let you start 4 types of analysis:

- phylobayes single thread as single slurm job
- phylobayes single thread as slurm job array
- phylobayes-mpi as single slurm job
- phylobayes-mpi as slurm job array

## 1 Prerequisities

- [Phylobayes](http://megasun.bch.umontreal.ca/People/lartillot/www/download.html)
- [Phylobayes MPI](https://github.com/bayesiancook/pbmpi)
- [Slurm](https://slurm.schedmd.com/) (will be optional in the future)
- [Python 3 and above](https://www.python.org/) 

## 2 Workflow

First unzip / untar your data in a folder and save the path to the directory:
```shell
export PHYLODATA="/absolute/path/to/data/directory"
```

Note the file extension of your sequence files. The default is `phylip`, save this in an environmental variable as well (without the preceding point):
```shell
export PHYLOEXTENSION="phylip"
```

Next, decide what you want to run: 

- one single slurm job
- job array or
- multiple job arrays.

### 2.1 Running one single slurm job

To run one single slurm job, use the `start_phylobayes.sbatch.sh`.

First open it in a text editor and change the `#SBATCH` settings according to your needs and your cluster configuration.

Also change the bottom part of the file, to adjust the waiting time between checking for convergence, number of chains, etc.

If everything is set, start the slurm job (we assume that there is a `sequence1.phylip` file in the previously uncompressed directory):
```shell
sbatch start_phylobayes.sbatch.sh $PHYLODATA/sequence1.phylip
```
This will start the `phylobayes_starter.py` Python script which will start Phylobayes chains and check them periodically.

### 2.2 Running one slurm job array

To run one a slurm job array, use the `start_phylobayes.sbatch_array.sh`.

First open it in a text editor and change the `#SBATCH` settings according to your needs and your cluster configuration.

Also change the bottom part of the file, to adjust the waiting time between checking for convergence, number of chains, etc.

If everything is set, start the slurm job:
```shell
sbatch start_phylobayes.sbatch_array.sh $PHYLODATA $PHYLOEXTENSION
```

This will start the `phylobayes_starter.py` Python script which will start Phylobayes chains and check them periodically.


### 2.3 Running multiple slurm job arrays

To run multiple slurm job arrays, use the `multiple_array_starter.bash.sh` (which uses the `start_phylobayes.sbatch_array.sh` in the background).

First open the `start_phylobayes.sbatch_array.sh` in a text editor and change the `#SBATCH` settings according to your needs and your cluster configuration.

Also change the bottom part of the file, to adjust the waiting time between checking for convergence, number of chains, etc.

You don't have to change the `multiple_array_starter.bash.sh`.

If everything is set, start the slurm job:
```shell
nohup bash ./multiple_array_starter.bash.sh $PHYLODATA $PHYLOEXTENSION &
```

This will loop through the directories found under `$PHYLODATA` and call `start_phylobayes.sbatch_array.sh` on them one by one.



