[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)


# Efficient Nested Simulation Experiment Design via the Likelihood Ratio Method

This archive is distributed in association with the [INFORMS Journal on Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [Efficient Nested Simulation Experiment Design via the Likelihood Ratio Method](XXX Some DOI link) by Ben Mingbin Feng and Eunhye Song.

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2022.0392

https://doi.org/10.1287/ijoc.2022.0392.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{feng2024nestedsimLR,
  author =        {Feng, Ben Mingbin and Song, Eunhye},
  publisher =     {INFORMS Journal on Computing},
  title =         {Efficient Nested Simulation Experiment Design via the Likelihood Ratio Method},
  year =          {2024},
  doi =           {10.1287/ijoc.2022.0392.cd},
  url =           {Available for download at https://github.com/mbfeng/2022.0392/},
}
```

## Description
The goal of this software is to demonstrate the effectiveness of the optimal likelihood-ratio nested simulation design, as proposed in the paper, comparing to other nested simulation designs in various applications.

## Building
R is required to run source codes included. Please see requirements.txt and install the required R libraries. In May 2024, the authors used the following version of R:

R version 4.3.3 (2024-02-29 ucrt) -- "Angel Food Cake"

The library versions are the ones used by the authors in May, 2024.

For the single-asset enterprise risk management (ERM) example in Section 7.1 of the paper, run the following scripts:
```
src\Section_7_1_SingleAssetERM\ERM_Single_Asset_1.R # produces Figure 1 & Figure 2 in Section 7.1
src\Section_7_1_SingleAssetERM\ERM_Single_Asset_2.R # produces Table 1 in Section 7.1
```

For the multi-asset enterprise risk management (ERM) example in Section 7.2 of the paper, the user should take three steps:
1. Set the directory in line 20 of src\Section_7_2_MultiAssetERM\MultiAssetERM_main.R
2. Run the following shell script
```
src\Section_7_2_MultiAssetERM\runfile.sh
```
3. Run the script src\Section_7_2_MultiAssetERM\MultiAssetERM_postprocessing.R in the directory set in step 1, which produces Table 2 and Table 3 in Section 7.2 of the paper.

For the newsvendor example in Section 7.3 of the paper, run the following scripts:
```
src\Section_7_3_Newsvendor\Newsvendor_NestedSim.R # produces Table 4 in Section 7.3
src\Section_7_3_Newsvendor\Newsvendor_VarRatio.R # produces Figure 3(a) in Section 7.3
src\Section_7_3_Newsvendor\Newsvendor_SimBudgetGrowth.R # produces Figure 3(b) in Section 7.3
```

## Results
Sample results are in the results sub-directory. For instance,
```
results\Section_7_1_SingleAssetERM\OuterDist.pdf # Figure 1(a) in Section 7.1
results\Section_7_1_SingleAssetERM\ConditionalMean.pdf # Figure 1(b) in Section 7.1
results\Section_7_1_SingleAssetERM\ERM_All_RedGreen.pdf # Figure 2(a) in Section 7.1
results\Section_7_1_SingleAssetERM\ERM_Zoomed_RedGreen.pdf # Figure 2(b) in Section 7.1
results\Section_7_3_Newsvendor\VarRatio.pdf # Figure 3(b) in Section 7.3
results\Section_7_3_Newsvendor\SimBudgetGrowth.pdf # Figure 3(b) in Section 7.3
```

## Replicating
Please follow the instructions in the "Build" section to replicate results in the paper.
