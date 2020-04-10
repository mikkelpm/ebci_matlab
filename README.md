# ebci_matlab

Matlab code for computing robust empirical Bayes confidence intervals (EBCIs)

**Reference:**
Armstrong, Timothy B., Michal Kolesár, and Mikkel Plagborg-Møller (2020)
"Robust Empirical Bayes Confidence Intervals"
https://scholar.princeton.edu/mikkelpm/ebci

See also R package "ebci": https://github.com/kolesarm/ebci

This version: 2020-04-09
Tested in Matlab R2019b on Windows 10 PC, 64-bit

## Content

**[ebci](ebci):** EBCI routines
- [cva.m](ebci/cva.m): robust EBCI critical value
- [ebci.m](ebci/ebci.m): main function for computing parametric or robust EBCIs
- [parametric_ebci_maxnoncov.m](ebci/parametric_ebci_maxnoncov.m): worst-case non-coverage probability of parametric EBCI
- Additional auxiliary functions

**[vignette](ebci):** example code with empirical illustration
- [cz.csv](vignette/cz.csv): Chetty & Hendren (2018) data set (obtained from https://opportunityinsights.org/data/)
- [run_ebci.m](vignette/run_ebci.m): annotated script with examples
