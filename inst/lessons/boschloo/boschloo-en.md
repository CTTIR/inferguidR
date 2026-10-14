# Two independent binomial groups: Boschloo's test

This original lesson uses synthetic counts and the accompanying `boschloo.R` script. It teaches the explicitly named **two-sided minlike** variant. It does not establish full release certification of the package.

## Question and design

Two independent groups have fixed sizes. Each participant contributes one binary event/non-event observation. Is the event probability different between groups? Write the event definition, group order, two-sided alternative and significance level before inspecting the counts. Paired observations, repeated events within participants and survey weights need another analysis.

Here group 1 has 7 events among 8 participants; group 2 has 1 among 8. The observed risk difference, group 1 minus group 2, is 0.75. This is an estimate, not a future planning assumption.

## Why this test, and which version?

Boschloo orders possible tables using a Fisher conditional p-value, then maximises the null probability of equally or more extreme tables over the common event probability. This uses the independent-binomial sampling design without conditioning the analysis on the observed total number of events. See the [original paper](https://doi.org/10.1111/j.1467-9574.1970.tb00104.x).

The two-sided **minlike** and **central** variants differ. Minlike uses the two-sided Fisher ordering; central combines the separately maximised one-sided tests. This lesson fixes minlike throughout analysis and planning. For a genuinely predeclared directional question, `greater` means p1 > p2 and `less` means p1 < p2. Do not switch direction after seeing the result. The [exact2x2 documentation](https://search.r-project.org/CRAN/refmans/exact2x2/html/boschloo.html) also distinguishes the two-sided conventions.

## Predict before running

Before executing the script, predict whether the unconditional p-value and the conditional Fisher p-value must be equal. Then predict what happens to future power when p1 changes from 0.8 to 0.7 while p2 remains 0.2. Record the prediction before viewing the output.

Run in R with inferguidR installed:

```r
lesson <- new.env()
sys.source(inferguidR::learning_lesson("boschloo", "en")$script, envir = lesson)
answer <- lesson$run_boschloo_lesson(tempfile("boschloo-example-"))
answer
```

The script records package/R versions and session information. It writes only to the requested output directory. It creates English/German JSON and TXT reports and PDF reports when Cairo is available; it does not install a renderer. Planning and sensitivity tables and a separate planning plot are also written.

## Read the results correctly

The executed example gives Boschloo p approximately 0.004180908 and Fisher p approximately 0.01010101. Both tests reject at 0.05 here; this is not a claim that two-sided conventions are interchangeable. Compare the sampling assumptions and declared variant, not merely the smaller number.

`p_upper` is an upper bound from numerical optimisation. It is **not** a sampling confidence limit for the p-value, risk difference or event probability. The script separately gives two unadjusted marginal 95% Clopper–Pearson intervals, approximately [0.47349, 0.99684] and [0.00316, 0.52651]. These are neither a simultaneous interval pair nor a confidence interval for the difference. Their overlap is not a test of equal probabilities. This example does not compute a sampling interval for the risk difference.

## Plan a future study using the same variant

The future assumptions are p1=0.8, p2=0.2, alpha=0.05 and power=0.8, with equal group sizes. They are hypothetical design inputs, deliberately distinct from the observed proportions. The bounded search considers total N at most 40 and checks three subsequent designs for local stability.

The first qualifying design is n1=n2=10, total N=20, with power approximately 0.8053902. The reported stability window also passes. This is a minimum along the stated allocation sequence and search domain; the local stability window is not a proof for every larger design. Exact-test power can have discrete jumps. The scatterplot therefore shows evaluated designs, not an interpolated sample-size guarantee.

At the same 10+10 design, power is approximately 0.621632, 0.805390 and 0.936660 for p1=0.7, 0.8 and 0.9 respectively, with p2 fixed at 0.2. These are design scenarios, not observed-data confidence intervals or retrospective evidence for the hypothesis.

## Related teaching material and explicit gaps

The [CTTIR Fisher tutorial](https://cttir.github.io/tutorials/tutorials/inference/fisher-exact-test.html) is **related background only**, not a Boschloo lesson. Its displayed numerical output does not match its shown table (known issue SP10); rerun its code rather than quoting its numbers. [Pinned source fallback](https://github.com/CTTIR/tutorials/blob/7d8d01c4edabd294f88568c4e370ab83dc09d313/tutorials/inference/fisher-exact-test.qmd). No CTTIR Boschloo lesson was found in the audited source corpus. This original local lesson fills that specific gap; other teaching topics remain separate work.

A report should name the sampling design, event orientation, variant, observed effect, p-value, numerical status and missing interval explicitly. Share only reviewed output: a successful teaching run is not full scientific certification.
