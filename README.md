# ebci_matlab

Matlab code for computing robust empirical Bayes confidence intervals (EBCIs)

**Reference:**
Armstrong, Timothy B., Michal Kolesár, and Mikkel Plagborg-Møller (2022), "[Robust Empirical Bayes Confidence Intervals](https://doi.org/10.3982/ECTA18597)", *Econometrica* 90(6), 2567-2602 (also available as [arXiv:2004.03448](https://arxiv.org/abs/2004.03448))

See also the [ebci](https://github.com/kolesarm/ebci) R package by Michal Kolesár (who deserves all intellectual credit for the code structure) and the [ebciStata](https://github.com/kolesarm/ebciStata) Stata package

Tested in Matlab R2021a on Windows 10 PC, 64-bit

## Contents

**[ebci](ebci):** EBCI routines
- [cva.m](ebci/cva.m): robust EBCI critical value
- [ebci.m](ebci/ebci.m): main function for computing parametric or robust EBCIs
- [parametric_ebci_maxnoncov.m](ebci/parametric_ebci_maxnoncov.m): worst-case non-coverage probability of parametric EBCI
- Additional auxiliary functions

**[vignette](ebci):** example code with empirical illustration
- [cz.csv](vignette/cz.csv): Chetty & Hendren (2018) data set (obtained from https://opportunityinsights.org/data/)
- [run_ebci.m](vignette/run_ebci.m): annotated script with examples

## Requirements

Matlab's Optimization Toolbox must be installed.

## Acknowledgements

This software package is based upon work supported by the National Science Foundation under grant numbers SES-2049765 (Armstrong), SES-22049356 (Kolesár), and SES-1851665 (Plagborg-Møller), and by work supported by the Alfred P. Sloan Research Fellowship (Kolesár).
