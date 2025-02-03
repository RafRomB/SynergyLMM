## Resubmission

This is a resubmission. In this version I have:

* Converted the DESCRIPTION title to title case.

* Reduced the title length to less than 65 characters.

* Added more details about the package functionality and implemented methods in the Description text of the DESCRIPTION file.

* Added references describing the methods in the package in the Description field of the DESCRIPTION file.

* Removed \dontrun and unwrapped the corresponding examples.

* Substituted print()/writeLines() with message()/warning(), or used an additional 'verbose' argument in the functions to allow users to suppress printing in the following functions: R/APrioriPwr.R; R/CookDistance.R; R/lmmModel.R; R/lmmSynergy.R; R/logLikSubjectDisplacements.R; R/ObsvsPred.R; R/PwrSampleSize.R; R/PwrTime.R; R/ranefDiagnostics.R; R/residDiagnostics.R


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
