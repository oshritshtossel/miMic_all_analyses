[![DOI](https://zenodo.org/badge/766421687.svg)](https://zenodo.org/doi/10.5281/zenodo.10952811)

# miMic - overall analyses code 

This code is attached to the paper "mi-Mic: a novel multi-layer statistical test for microbiota-disease associations".
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
1. **"fig2_mother_daughter_distributions.py"**-
   - Daughter's distribution (pink) vs. mother's distribution (dark pink) in the regime of (0,0,0). The mother's distribution is noticeably narrower, with approximately half that of the daughter's distribution **(A)**.

2. **"fig2_leaf_confidence_vs_mimic_confidence.py"**-
   - Comparison of the leaf confidence with miMic's confidence based on analytical integral calculations over 3 different regimes: regime (0,0,0) **(B)**, regime (0,0,μ) **(C)**, regime (0,μ,α*μ) **(D)**. The lines represent estimated slopes, denoted as S. In D, different colors represent varying levels of connection between the sisters, controlled by α values (0.25, 0.5, 1, 2, 4).
   - Note that for the creation of:
      -  **B:** <code style="background-color: lightgrey;">STATE = '000'</code>.
      -  **C:** <code style="background-color: lightgrey;">STATE = '00mu'</code>.
      -  **D:** <code style="background-color: lightgrey;">STATE = '0mumu'</code>.

   - Note that the data for this simulation can be created by <code style="background-color: lightgrey;">Simulations/hierarchical_simulations.py</code>.

        
3. **"fig2_hist_sisters_distribution.py"**-
   -  Histogram illustrating the distribution of inner sisters-label SCCs across different cohorts. The black line represents the zero line, and the dashed pink line represents the average of the distribution, indicating a right-skewed distribution, with most sisters showing a consistent positive correlation with the label **(E)**.

4. **"fig2_fp_tp_bar_plots_hierarchial_simulations.py"**-
   - A comparison between the number of samples and the number of FPs for miMic and Mann-Whitney leaf simulations based on the regime of (0, 0, 0) **(F)** and the regime of (0, 0, μ), where μ was set to 1. The lightest pink color represents the Leaf-C model, the middle pink represents the Leaf model, and the darkest pink represents the miMic model.
   - **H.** Comparison between the number of samples and the number of TPs for miMic and Leaf Mann-Whitney simulations based on the regime of (0, 0, μ). The color coding is similar to the color coding of **F-G**.
   - Note that the data for this simulation can be created by <code style="background-color: lightgrey;">Simulations/hierarchical_simulations.py</code>.


5. **"fig2_surface_pos_neg_hierarchial_simulations.py"**-
   - **I - J.** Comparison between the FPs **(I)** and TPs **(J)** of the miMic model and the Leaf model in the simulation of the regime (0, 0, μ) over different numbers of samples and different values of μ. The Leaf model shows a higher number of FPs compared to miMic, while miMic's TPs are similar to Leaf's TPs.

   - **K - L.** Comparison between the FPs **(K)** and TPs **(L)** of the miMic model and the Leaf model in the simulation of the regime (0, μ, α*μ) over different numbers of samples and different values of α. The Leaf model exhibits a higher number of FPs compared to miMic, while miMic's TPs are similar to the Leaf's TPs.
     
   - Note that for the creation of:
      - **I:** <code style="background-color: lightgrey;">TASK = 'A'</code>
      - **J:** <code style="background-color: lightgrey;">TASK = 'B'</code>
      - **K:** <code style="background-color: lightgrey;">TASK = 'C'</code>
      - **L:** <code style="background-color: lightgrey;">TASK = 'D'</code>

   - Note that the data for this simulation can be created by <code style="background-color: lightgrey;">Simulations/hierarchical_simulations.py</code>.

6. **"fig2_microbiome_oriented_f1_scores.py"**-
   - The F1 scores of different differential abundance (DA) methods are depicted across three distinct setups of microbiome-oriented simulations **(M)**.
   - Note that the data for this simulation can be created by <code style="background-color: lightgrey;">Simulations/microbiome_oriented_simulations.py</code>.

### Fig_3 - Validation of miMic vs. SOTA models on real-world datasets.
1. **"fig3_plotly_RSP_comparison.py"**-
   - Comparative analysis of different Differential Analysis (DA) methods as a function of RSP(β) over 16S cohorts **(B)** and Whole Genome Sequencing (WGS) cohorts **(C)**. Each color represents a specific model: orange for DeSeq2 (light without FDR correction and dark with FDR correction, denoted as DeSeq2-C), yellow for LefSe, green for ANCOM (light without FDR correction and dark with FDR correction, referred to as ANCOM-C), blue for LINDA (light without FDR correction and dark with FDR correction, denoted as LINDA-C), brown for ada-ANCOM, and pink for miMic (pink for log SUB-PCA MIPMLP preprocessing and purple for relative mean MIPMLP preprocessing). Each line illustrates the average RSP(β) score across all cohorts (12 16S in **B** and 8 WGS in **C**). The light shadows surrounding each line represent standard errors calculated over 10 simulations of the shuffled models across all cohorts (12 16S in **B** and 8 WGS in **C**).

2. **"fig3_plotly_RSP_taxonomy.py"**-
   - Comparison of different starting taxonomy levels of the miMic test and their corresponding RSP(β) scores over 16S cohorts **(D)** and WGS cohorts **(E)**. Each taxonomy level is indicated by a different line style (1 for kingdom, 2 for phylum, 3 for class, 4 for order, 5 for family, 6 for genus, and 7 for species). Typically, the best RSP scores are achieved in the first two taxonomy levels. The inner bar plot presents the number of cohorts in which miMic is deemed significant when commencing the Mann-Whitney test at each taxonomy level.

3. **"fig3_scatter_sister_scc_rsp.py"**-
   - Scatter plot of the sister's-labels SCC vs. the RSP(1) score. A significant positive correlation of 0.588 is observed between the SCCs of sister labels and the model's performance **(F)**.

### Fig_4 - Consistency and robustness analysis of miMic
1. **"fig4_plot_consistency_within_models.py"**-
   - Within-study differential abundance consistency analysis across multiple tools. The percentage of total significant features is plotted against the number of tools that identified the feature as significant. Results are shown for the miMic model (pink) and the average of all state-of-the-art (SOTA) models (white). Refer to Supp. Mat. Fig. S3 for a detailed analysis of all 13 tools. The total number of significant features identified by each tool is provided in the legend. miMic demonstrates slightly higher consistency compared to the average of all SOTA models **(A)**.

2. **"fig4_plotexpected_vs_observed_consistency_within_datasets.py"**-
   - The percentage of significant species is plotted against the number of studies where each species was identified as significant, conducted on five inflammatory bowel disease (IBD) cohorts. Results for miMic are depicted in pink **(B)**, while those for the average state-of-the-art (SOTA) model are shown in white **(C)**. The expected results are presented in black (see Methods). Additionally, a parallel analysis on shuffled labels is provided for the ANCOM-BC2 model (green) within **(C)**. The models' performance exceeds that of the expected random model. However, certain tools, such as ANCOM-BC2, exhibit artificially consistent results, as indicated in the inner plot of **(C)**. For a comprehensive analysis of all 13 tools, refer to Supplementary Material Fig. S4.

   - Note that in order to find the intersection, one should use <code style="background-color: lightgrey;">Additional evaluations tools/consistency_within_datasets.py</code>.
    
   - Note that in order to compute the expected simulations, one should use <code style="background-color: lightgrey;">Additional evaluations tools/expected_simulations.py</code>.

3. **"fig4_generic_features_robustness.py"**-
   - Sensitivity robustness assessment. The heatmap illustrates Spearman correlation coefficients (SCCs) between each generic dataset characteristic and the percentage of significant taxa identified by each tool per dataset. Positive correlations are depicted in red, while negative correlations are shown in blue. Stars indicate a significant correlation (p-value < 0.05). miMic demonstrates robustness across all tested generic features in 16S datasets. For parallel analyses conducted on 16S and whole-genome sequencing (WGS) cohorts, detailing the percentage of significant taxa identified by each tool per dataset and RSP score, refer to Supplementary Material Fig. S5.

### Fig_5 - Differential abundance analysis results are visualized on a cladogram for the IBD cohort.
   - This plot is created by the [miMic PyPi](https://pypi.org/project/mimic-da/).
     
   - Each color represents the sign of the Mann-Whitney score (blue for positive scores, red for negative scores, and grey for non-significant taxa). The node size corresponds to -log10(p-value) from the Mann-Whitney test in miMic. The node shape represents its origin of significance: spheres were identified by both miMic and the Mann-Whiteny test on leaves, circles were identified by miMic only, and squares were identified by only the Mann-Whitney test. The colors represent the taxonomic family of each node.

### Fig_6 - miMic's plots - example on IBD cohort.
   - These plots is created by the [miMic PyPi](https://pypi.org/project/mimic-da/).

   - **A.** Bar plot illustrating the taxonomy levels in the miMic test vs. the number of significant findings in a real run (RP) shown in blue, and in a shuffled run (SP) shown in red. The highest bar plot represents the actual RP vs SP of the selected taxonomy level of miMic combined with the leaves test as explained in the Methods. Taxonomy levels used for the a priori nested ANOVA test are shaded in grey. The number of RPs significantly exceeds the number of SPs.

   - **B.** Interaction between significant taxa found in miMic. Each taxon is colored according to its significant family color, similar to Fig. 5 above. Each node shape represents the taxon's order. An edge is drawn between two nodes if their Spearman correlation coefficient (SCC) is above 0.3 (user-adjustable) and its p-value < 0.05 (user-adjustable). The width of the edge corresponds to its SCC. A blue edge represents a positive relation, while a red edge represents a negative one.

   - **C.** Analysis of significant positive and negative relations within taxonomic families. The y-axis displays significant families in the cohort (defined by a family that has at least 1 significant descendant), while the x-axis shows the count of positive relations within a family in blue or the count of negative relations within a family in red. Each family is colored according to its color in the interaction network in **(B)** and the cladogram of correlations in Fig. 5 above.


## Cite us

Shtossel, Oshrit, and Yoram Louzoun. "miMic-a novel multi-layer statistical test for microbiome disease." (2023).

## Contact us

[Oshrit Shtossel](oshritvig@gmail.com)





    
