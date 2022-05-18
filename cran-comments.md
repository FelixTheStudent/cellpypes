## Resubmission
This is a resubmission (of the previous resubmission). In this version I fixed:

  * DESCRIPTION file:
    - removed "in R" from package title
    - use single quotes in description (e.g. 'cellpypes')
    - added references describing the methods in this package
      (doi and https)
  * made messages suppressible with 'if(verbose) cat(...)'.
    Exception: function 'pype_code_template' still uses cat because its sole
    purpose is to print a code template to the console during interactive use.
  * no longer set a seed to a specific number within function 'find_knn'.


## Previous Resubmission (0.1.2)
These are the things I fixed due to Uwe Ligges' feedback:

* removed the VignetteBuilder field as there are no vignettes in this package
* updated links to CRAN packages, now using the canonical URL
* updated link to Bioconductor package, now starting with 'https'


## Test environments
* local R installation (ubuntu), R 4.1.3
* win-builder (devel)
* win-builder (release)
* win-builder (oldrelease)
* macOS 10.13.6 High Sierra, R-release, CRAN's setup (rhub::check())


## R CMD check results

0 errors | 0 warnings | 1 notes

   
* Note1: New submission.


## CRAN incoming feasibility

* Note: Possibly misspelled words in DESCRIPTION:
    cellpypes  (that's the name of this package, no spelling error)
    pseudobulks (that's jargon/terminology, no spelling error)
