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

## Statistical analysis, comparisons and visualizations

The analyses are presented according to the figures in the manuscript.

### Fig_1 - Challenges in Difference Analysis of Microbiome
1. **"fig1_pie_chart.py"**-
   - Pie diagram illustrating the usage of differential analysis methods within the field during 2023. The colors represent different methods: LefSe (yellow), DeSeq2 (orange), ANCOM (sea-green), ANCOM-BC2 (green), ALDEx2 (grey), and LINDA (light blue),**(A)**.
     
2. **"fig1_rp_vs_sp.py"**-
   - Scatter plot comparing the performance (SP vs. RP) of popular methods across more than 20 different microbial datasets. Each shape represents a distinct cohort with colors matching those indicated in the pie chart **A**. Dark colors indicate methods with FDR corrections, while light colors represent methods without any corrections. The dashed grey line represents the y = x line, where the SP rate is similar to the RP rate.Pink and purple colors in the upper left corner represent miMic and miMic-relative, respectively. In many cases (excluding miMic), the RP and SP rates are similar **(B)**.
   - Note that the  <code style="background-color: lightgrey;">csv</code> for creating this plot is  <code style="background-color: lightgrey;">all_corrected_fp_tp.csv</code> in the <code style="background-color: lightgrey;">Fig_1/Results</code> folder.

3. **"fig1_sis_correlations.py"**-
   - Inner sisters'-labels SCCs over the different datasets. The stars represent the significance of the correlations, such that *-p<0.05, **-p<0.01, ***-p<0.001 **(C)**.
   - Note that the  <code style="background-color: lightgrey;">csv</code> for creating this plot is <code style="background-color: lightgrey;">all_corrs.csv</code> in the <code style="background-color: lightgrey;">Fig_1/Results</code> folder.

### Fig_2 - Validation of miMics' assumptions on analytical models and simulations.
1. **"fig2_mother_daughter_distributions.py"**
   - Daughter's distribution (pink) vs. mother's distribution (dark pink) in the regime of (0,0,0). The mother's distribution is noticeably narrower, with approximately half that of the daughter's distribution **(A)**.

2. **"fig2_leaf_confidence_vs_mimic_confidence.py"**
   - Comparison of the leaf confidence with miMic's confidence based on analytical integral calculations over 3 different regimes: regime (0,0,0) **(B)**, regime (0,0,μ) **(C)**, regime (0,μ,α*μ) **(D)**. The lines represent estimated slopes, denoted as S. In D, different colors represent varying levels of connection between the sisters, controlled by α values (0.25, 0.5, 1, 2, 4).

    
