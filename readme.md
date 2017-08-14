
# Single arm two-stage studies: Improved designs for molecularly targeted agents

Contained within this repository are programs for running single arm two-stage trials with two subgroups defined by a biomarker with an assumed ordering into marker-positive and marker negative patients. Files in the Programs folder is simulation and utilities code to run these simulations. Files at the top level were used to create the tables and figures in the paper.

## Dependencies

Some of the simulation rely on the package trialSim (v1.0) which is available on github only. This package is a modular system for simulating trials. All other dependencies are available on CRAN. Further instructions can be found in the file: 0_instructions_and_runPaper.R


## Acknowledgements

We would like to thank Yong Zang and Ying Yuan who kindly provided code to generate designs for their sequential enrichment design.

## Author

Programmed by: Peter Dutton (Centre for Statistic in Medicine, University of Oxford)

## References

Zang Y, Yuan Y. Optimal sequential enrichment designs for phase II clinical trials. Statistics in Medicine. 2017;36(1):54-66.

Pusztai L, Anderson K, Hess KR. Pharmacogenomic Predictor Discovery in Phase II Clinical Trials for Breast Cancer. American Association for Cancer Research. 2007;13(20):6080-6086.

Parashar D, Bowden J, Starr C, Wernisch L, Mander A. An optimal stratified Simon two-stage design. Pharm Stat. 2016;15(4):333-340.

Jones CL, Holmgren E. An adaptive Simon Two-Stage Design for Phase 2 studies of targeted therapies. Contemporary clinical trials. 2007;28(5):654-661.

Whitehead J, Valdés-Márquez E, Johnson P, Graham G. Bayesian sample size for exploratory clinical trials incorporating historical data. Statistics in Medicine. 2008;27(13):2307-2327.