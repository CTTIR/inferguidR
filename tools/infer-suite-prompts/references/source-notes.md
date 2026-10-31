# Source baseline and primary references

These observations informed the stack on 31 October 2026. The execution agent must inspect the actual checkout and recheck changing external interfaces. Paths and findings are not an instruction to overwrite newer work.

## inferguidR baseline

Repository: https://github.com/CTTIR/inferguidR

Observed commit: `e91938f77eaa991b69989b121bfd32c124c95db0`, dated 29 October 2026. DESCRIPTION declares development version 0.0.0.9000 and R >= 4.1.0 in this snapshot.

| Area | Source | Consequence |
|---|---|---|
| Power method definitions | [core-method.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/core-method.R), [methods-t.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/methods-t.R) | methods-* includes forward power/design machinery; ownership must be decided per function |
| Planning solvers | [num-modes.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/num-modes.R) | Keep planning/inverse search in guid |
| Namespace providers | [core-registry.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/core-registry.R), [core-api.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/core-api.R) | Extraction needs explicit public interfaces |
| Advisor scripts | [advisor-code.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/advisor-code.R) | Teaching scripts are not a uniform production API |
| Method distinctions | [analysis-variants.json](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/inst/advisor/analysis-variants.json) | Preserve observed/planning variant identity, including Welch versus pooled t |
| Advisor planning bridge | [advisor-power.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/advisor-power.R) | Existing public bridge is specifically scoped to Boschloo variants |
| Boschloo test | [boschloo-test.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/boschloo-test.R) | Numerical p bounds and unresolved decisions are distinct from sampling intervals |
| Boschloo planning | [boschloo-method.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/boschloo-method.R), [boschloo-design.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/boschloo-design.R) | Common rejection semantics and discrete allocation search must survive |
| Result envelopes | [core-envelope.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/core-envelope.R) | Observed/power/classification kinds and hashes already exist |
| Portable reporting | [report-document.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/report-document.R), [report-export.R](https://github.com/CTTIR/inferguidR/blob/e91938f77eaa991b69989b121bfd32c124c95db0/R/report-export.R) | Preserve old documents and exports; observed-test path is Boschloo-specific in this snapshot |

## R package and verification guidance

- [Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html); [official manual mirror](https://stat.ethz.ch/R-manual/R-devel/doc/manual/R-exts.html). Consult the manual matching the actual runtime for package tests, dependencies, build/check, S3 registration and portable resources.
- [stopifnot](https://stat.ethz.ch/R-manual/R-devel/library/base/html/stopifnot.html).
- [all.equal](https://stat.ethz.ch/R-manual/R-devel/library/base/html/all.equal.html).
- [tools condition assertions](https://stat.ethz.ch/R-manual/R-devel/library/tools/html/assertCondition.html).

New tests can use standard R assertions; contributed testing frameworks are not a prerequisite for ordinary R package tests. Verify installed tarballs and actual exit/check results.

## Reporting and plots

- [DescrTab2 CRAN status](https://stat.ethz.ch/CRAN/web/packages/DescrTab2/index.html): retrieved status says archived 18 October 2025. Recheck before implementation.
- [DescrTab2 descr interface](https://imbi-heidelberg.github.io/DescrTab2/reference/descr.html) and [source](https://github.com/imbi-heidelberg/DescrTab2/blob/master/R/descr.R): descriptive and test functionality coexist. Printing controls are not calculation controls.
- [ggplot2 internals](https://ggplot2-book.org/internals.html): scales transform/censor layer data before statistical processing. Precomputed-result annotation separates presentation from inference.
- [ggpubr stat_compare_means](https://rpkgs.datanovia.com/ggpubr/reference/stat_compare_means.html): current comparisons path documents independent unadjusted pairwise tests.
- [ggpubr stat_pvalue_manual](https://rpkgs.datanovia.com/ggpubr/reference/stat_pvalue_manual.html): a possible optional existing renderer for precomputed results.
- [gtsummary precomputed summary interface](https://www.danieldsjoberg.com/gtsummary/reference/tbl_ard_summary.html) and [custom statistic hook](https://www.danieldsjoberg.com/gtsummary/reference/add_stat.html): optional integration opportunities, not mandatory core dependencies.

## GitHub operations

- [gh repo create](https://cli.github.com/manual/gh_repo_create).
- [GITHUB_TOKEN scope](https://docs.github.com/en/actions/concepts/security/github_token).
- [Workflow authentication](https://docs.github.com/en/actions/tutorials/authenticate-with-github_token).

Verify exact owner, visibility and permissions. Private cross-repository access needs an appropriate authorized credential or artifact route. Creating workflow configuration is not evidence that the workflow executed.

