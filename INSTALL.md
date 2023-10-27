# GAinSAW Installation and Setup

## Requirements

  * GAinSAW should run on any platform that supports [Singularity](https://apptainer.org/) (Linux, MacOS), although we have only tested on Linux.

  * GAinSAW can execute on a single processor machine, but realistically you would want to have 10-20 cores available.

  * GAinSAW requires a considerable amount of memory (>32 GB) and free disk space (~50 GB for large genomes).


## Installation as a singularity container

Assuming _git_ and  _singularity_ are installed on your system, you can get the
GAinSAW code from GitHub and the container from our
[Singularity Hub](http://BrendelGroup.org/SingularityHub) and run _xgainsaw_ as follows:

```bash
git clone https://github.com/vpbrendel/GAinSAW
cd GAinSAW
wget https://BrendelGroup.org/SingularityHub/GAinSAW.sif
alias rws="singularity exec -e -B${PWD}  ${PWD}/GAinSAW.sif"
rws xgainsaw -h
```

For a gentle introduction to singularity, see our group
[handbook article](https://github.com/BrendelGroup/bghandbook/blob/master/doc/06.2-Howto-Singularity-run.md).
Alternatively, install auxillary programs and the _gainsaw_ python package following the instructions below. 

## System-wide installation of auxiliary programs

**GAinSAW** use via the singularity container is highly recommended, with no known
drawbacks.
However, if desired, you can of course install all the required software and
packages individually on your computer system.
The Singularity [definition file](./GAinSAW.def) in this repository should serve as
a guide to perform such an installation and is replicated in the script [0Record-system](./0Record-system)
(requires _sudo_ privileges).


## Setting up a conda environment for _gainsaw_

A conda environment for _gainsaw_ can be set up with the script [0Record-setup](./0Record-setup)
in the _GAinSAW_ directory.
Then, please activate the _gainsaw_ environment and prepare data for testing as follows:
```bash
conda activate gainsaw


# 1. Get some tools we need (please read tools/0README for full instructions!):
#
cd tools
bash 0README
cd ..
export PATH=$PATH:${PWD}/tools

# 2. Download data we need for the usage examples:
#
#
xgainsaw download -i data_sources.txt -D downloads
sleep 30

# 3. Use xgainsaw to populate the data directory:
#
xgainsaw prepare -q mm10 \
	-f downloads/mm10/mm10.fa.gz \
	-z downloads/mm10/mm10.chrom.sizes \
	-A downloads/mm10/chromAlias.txt.gz \
	-C downloads/mm10 \
	-d data -c work.conf
xgainsaw prepare -q mm39 \
	-f downloads/mm39/mm39.fa.gz \
	-z downloads/mm39/mm39.chrom.sizes \
	-a downloads/mm39/GCF_000001635.27_GRCm39_genomic.gff.gz \
	-A downloads/mm39/chromAlias.txt.gz \
	-C downloads/mm39 \
	-d data -c work.conf
xgainsaw prepare -q rn7 \
	-f downloads/rn7/rn7.fa.gz \
	-z downloads/rn7/rn7.chrom.sizes \
	-a downloads/rn7/GCF_015227675.2_mRatBN7.2_genomic.gff.gz \
	-A downloads/rn7/rn7.chromAlias.txt \
	-C downloads/rn7 \
	-d data -c work.conf

# 4. Pull the chr5 and chr7 data from mm39ToRn7.over.chain.gz for use in the
#    test examples:
#
chainFilter -t=chr5,chr7 data/liftovers/mm39ToRn7.over.chain.gz > data/liftovers/mm39chr5and7ToRn7.over.chain

# 5. Make the tools accessible. You may want to include the following into your
#    ~/.bashrc script (or equivalent):
#
echo -e "\n\nPlease execute the following in your working shell to make the downloaded"
echo -e "tools accessible. You may want to add this command to your ~/.bashrc"
echo -e "(or equivalent) file for future use:\n"
echo -e 'export PATH=$PATH:${PWD}/tools'
echo -e "\n\n"
```


## Finally

proceed to the [HOWTO](./HOWTO.md) document for usage examples.
