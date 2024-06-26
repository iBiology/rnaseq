# rnaseq

A pipeline for processing RNA-Seq data to identify gene differential expression.

# Example
```shell
source /work/08944/fuzzy/ls6/software/virtual-screening/venv/bin/activate

vs \
  --sdf /work/08944/fuzzy/share/chemical_library/docking/MolPort.ligand.core.sdf.gz \
  --pdb /work/08944/fuzzy/share/receptor/dc1_f4293/dc1_f4293.pdb \
  --center 5.860 4.130 6.280 \
  --outdir /work/08944/fuzzy/share/DC1 \
  --residue 476 \
  --time 50 \
  --project CHE23039
```

Successfully run the above commands will submit a job named `vs` to job queue for 
performing virtual screening. Output files will save to `/work/08944/fuzzy/share/DC1`. 
To fine tune the screening, see the detailed usage and options below.

## Usage
```shell
source /work/08944/fuzzy/ls6/software/virtual-screening/venv/bin/activate

$ vs -h
usage: vs [-h] --sdf SDF --pdb PDB [--center CENTER CENTER CENTER] [--filters FILTERS] [--consensus CONSENSUS] 
          [--outdir OUTDIR] [--flexible FLEXIBLE] [--size SIZE SIZE SIZE] [--top TOP] [--residue [RESIDUE ...]] 
          [--clusters CLUSTERS] [--poses_per_cluster POSES_PER_CLUSTER] [--time TIME] [--separate] [--nodes NODES] 
          [--name NAME] [--project PROJECT] [--email EMAIL] [--hold] [--quiet] [--verbose] [--debug]

A pipeline for perform virtual screening in an easy and smart way

options:
  -h, --help            show this help message and exit
  --sdf SDF             (Path, required) path to an SDF file containing prepared ligands
  --pdb PDB             (Path, required) path to a PDB file containing prepared receptor
  --center CENTER CENTER CENTER
                        (tuple[float, float, float], default=None) the X, Y, and Z coordinates of the center
  --filters FILTERS     (Path, default=None) path to a JSON file containing filters
  --consensus CONSENSUS
                        (float, default=20)
  --outdir OUTDIR       (Path, default=None) path to a directory for saving outputs
  --flexible FLEXIBLE   (str, default=) path to a flexible PDB file
  --size SIZE SIZE SIZE
                        (tuple[int, int, int], default=(15, 15, 15)) the size in the X, Y, and Z dimension (Angstroms)
  --top TOP             (float, default=20) fraction or number of top poses need to retain
  --residue [RESIDUE ...]
                        (list[int], default=None) list of integers representing the residue index
  --clusters CLUSTERS   (int, default=100) number of clusters to clustering
  --poses_per_cluster POSES_PER_CLUSTER
                        (int, default=5) maximum number of poses per cluster
  --time TIME           (float, default=50) time for molecule dynamics simulation in nanosecond
  --nodes NODES         (int, default=8) number of nodes need to request
  --name NAME           (str, default=vs) name of the virtual screening job
  --project PROJECT     (str, default=MCB23087) nmme of project you would like to be charged
  --email EMAIL         (str, default=) your email address for sending emails
  --separate            (bool, default=False) whether to separate the docking and molecule dynamics jobs in to submissions
  --hold                (bool, default=False) flag will only generate submit script but hold for submitting
  --quiet               (bool, default=False) whether to enable quiet mode to suppress debug messages
  --verbose             (bool, default=False) whether to enable verbose mode to print debug messages
  --debug               (bool, default=False) whether to enable debug mode

```
