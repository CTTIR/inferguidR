# Guidance provider integration

The numerical and canonical JSON implementations have one owner in infercalcR. Guidance retains planning, its historical document readers, and local structural validation. `R/suite-providers.R` contains explicit calls to named public calc providers. There is no runtime source loading or namespace discovery.

Apply `NAMESPACE.proposed` to the coordinating guid package (or remove the registrations listed in `NAMESPACE.remove-S3`, append `NAMESPACE.imports`, retain the existing exports). The calc namespace owns S3 methods for transferred classes. No DESCRIPTION change is needed beyond the existing infercalcR import.

`removed-symbols.csv` records definitions removed from guidance. `retained-structural-symbols.csv` records schema constructors, condition helpers and argument validation retained here. Guidance's `.ig_to_plain` still dispatches its own power, advisor, envelope and report classes before calling calc's canonical encoder/hash provider; this is required for historical JSON compatibility.

`calc-tests/testthat/` is the proposed test transfer to calc. Copy those files into calc/tests/testthat, retain calc's base R test scripts, and add an ordinary testthat runner guarded by availability (testthat/withr remain Suggests). It contains the original independent references and numerical assertions with package prefixes changed to the implementation owner. `test-moves.csv` records transferred files/blocks. Two private metric edge snippets from `test-classif-review.R` are combined in `test-kernels-classif-review.R`; all surrounding public guidance integration checks remain here. The calculator and confidence-boundary reference fixtures move byte-for-byte with their numerical tests; no retained guidance test uses them.

Mixed tests use public providers for control-record equivalence and canonical scalar string encoding, preserving the same exact-value assertions. No provider API was added solely for tests.

The migration scripts are development-time records tied to the specified baseline, not package runtime components. Execution evidence and results are in the coordinating P07 handoff and ignored run directory.
