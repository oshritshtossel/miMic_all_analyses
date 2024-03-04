# miMic - overall analyses code 

This code is attached to the paper "miMic - a novel multi-layer statistical test for microbiome disease".
**MiMic (Mann-Whitney iMage Microbiome)**, a novel approach for microbiome differential abundance analysis, tackles the key challenges of such statistical tests: 
a large number of tests, sparsity, varying abundance scales, and taxonomic relationships. Mimic first converts microbial counts to cladogram of means. 
It then applies an apriori tests on the upper levels of the cladogram to detect overall relationships. Finally, it performs a Mann Whitney test on paths 
that are consistently significant along the tree or on the leaves. 
MiMic has much higher true to false positives ratios than existing tests, as measured by a new real to shuffle positive score.

## How to apply miMic?

miMic's code is available at a separate [GitHub](https://github.com/oshritshtossel/miMic), a [PyPi](https://pypi.org/project/mimic-da/), and a [website](https://micros.math.biu.ac.il/Home).

## Comparison to other state-of-the-art tools

During our work, we compared our results to other state-of-the-art tools in the field. All the models' codes are available in separate Google colabs.

### ANCOM and DeSeq2

The implementation of the model (in python) over the datasets tested in the manuscript is available in in the following [notebook](https://colab.research.google.com/drive/1K4lS3mtCCeGYz-XpwcEdNxPX6OCDWF0b?usp=sharing).

### ANCOM-BC2, ALDEx2, structSSI, and ada-ANCOM
The implementation of the model (in R) over the datasets tested in the manuscript is available in in the following [notebook](https://colab.research.google.com/drive/1KRH4eFUW59KuXxmxJOAp_uKgAQguftFG?usp=sharing).

### LINDA
The implementation of the model (in R) over the datasets tested in the manuscript is available in in the following [notebook](https://colab.research.google.com/drive/1nzAEIZ27FGpRwylBeJG00Ci_K7VULPAr?usp=sharing).

### LEfSe
LEfSe was applied via its [Galaxy module](https://huttenhower.sph.harvard.edu/lefse/).

## Simulations

To evaluate miMic, 2 types of simulations are available:
- General hierarcial simulations (see <code style="background-color: lightgrey;">Simulations/hierarchical_simulations.py</code>).
- Microbiome-oriented simulations (see <code style="background-color: lightgrey;">Simulations/microbiome_oriented_simulations.py</code>).
