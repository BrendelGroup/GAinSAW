# GAinSAW - *G*enome *A*lignment cha*in* file *S*election and *A*nnotation *W*orkbench

Continuing improvements in DNA sequencing technologies and whole genome assembly strategies are enabling large-scale comparative genomics studies that probe how regions in one ("query") genome map to other ("target") genomes.
Because of the large size of the genomes (compared to sequences in traditional multiple sequence alignments) and potentially complex genome rearrangements, including translocations, inversions, and duplications, such studies pose difficult computational problems.
[Progressive Cactus](https://www.nature.com/articles/s41586-020-2871-y) provides a state-of-the-art solution.

Results of pairwise genome alignments are typically represented in the [UCSC Chain Format](https://genome.ucsc.edu/goldenPath/help/chain.html).
The [UCSC LiftOver Tool](https://genome.ucsc.edu/cgi-bin/hgLiftOver)
can be used to map query genome coordinates to target genome coordinates, as long as the query coordinates fall into an alignment block in the input chain file.
The figure below depicts points that can be lifted (long arrows) and points that fall into query-unique segments (short arrows).

![PoinSetConservation](./PointSetConservation.png)

While successful as a tool to discover genomic synteny, chained alignments are [not recommended for SNP liftover](https://genome.ucsc.edu/FAQ/FAQreleases.html#snpConversion).
Necessarily, chained alignments represent best-guess global alignments.
Determination of orthologous points will only be confident if the local alignments around the points are unambiguous.

**GAinSAW** was developed as a tool to liftover sets of points from a query genome to a target genome and then filter the points with respect to annotation and alignment confidence criteria.
The code conforms to our [RAMOSE](https://brendelgroup.github.io/)
philosophy: it generates __reproducible__, __accurate__, and __meaningful__
results; it is __open__ (source) and designed to be __scalable__ and
__easy__ to use.


## Quick Start

The simplest way to get going is to use the GAinSAW
[Singularity](https://apptainer.org/) container available from our
[Singularity Hub](http://BrendelGroup.org/SingularityHub/) site; e.g.:


```bash
cd
git clone https://github.com/vpbrendel/GAinSAW
cd GAinSAW
wget https://BrendelGroup.org/SingularityHub/GAinSAW.sif
alias rws="singularity exec -e -B ~/GAinSAW  ~/GAinSAW/GAinSAW.sif"
rws xgainsaw -h
```

In the above example, you clone this repository into your Linux home directory,
go into the thus created GAinSAW directory, download the GAinSAW Singularity
container, define the bash alias _rws_ ("run with singularity"), and check that
everything works by showing help information for _xgainsaw_.

Of course this assumes that you have [Apptainer/Singularity](https://apptainer.org/) installed on your system.
Check whether there is a package built for your system.
Otherwise, follow the instructions to [install Singularity from source code](https://apptainer.org/user-docs/master/quick_start.html#quick-installation-steps).


## Realistic Start

Please find detailed installation instructions and options in the
[INSTALL](./INSTALL.md) document.
Once all preparatory steps are taken care of, see the [HOWTO](./HOWTO.md)
document for usage examples.


## Reference

__Wenshu Chen and Volker P. Brendel (2023)__
  **GAinSAW**: **G**enome **A**lignment cha**in** file **S**election and **A**nnotation **W**orkbench;
  _in preparation_


## Contact

Please direct all comments and suggestions to [Wenshu Chen](<mailto:wenschen@indiana.edu>) or [Volker Brendel](<mailto:vbrendel@indiana.edu>) at [Indiana University](http://brendelgroup.org/).
