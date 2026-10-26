# Advisor review, 26 October 2026

These tables make the scope of the roadmap's advisor review inspectable.
They describe the source hashes in `source-identity.json`; they are not a
scientific approval, a study-specific recommendation or a release certificate.

- `variant-semantics.csv`: all 69 variants, including 49 primary routes,
  19 override-only variants and one rejected-only variant. The review covers
  estimands, assumptions, inference, missingness, executable-template hashes,
  source relationships and unresolved judgments. All 68 executable templates
  were read. Corrections include the Boschloo interval promise, an implicit
  signed-rank null, paired-data missingness wording, regression ratio overflow
  and the mixed-data missingness assumption.
- `bibliography.csv`: all 30 reference identities reconciled with publisher,
  DOI or official records, including the 26 whose retained historical labels
  say they had not been re-verified online. Editions, original dates, preprints
  and book reviews were distinguished. Bibliographic existence does not show
  that every rule follows from the work; full-content review is explicitly
  not established for the whole bibliography.
- `learning-gaps.csv`: all 19 variants without a registered CTTIR lesson.
  Eighteen have scoped external method, formula or API reference candidates;
  the variance route is intentionally rejected-only. Candidates have not been
  registered as exact lessons. Their read scope and caveats remain visible.

The full [Delacre correction](https://doi.org/10.5334/irsp.661) was compared with
the current Welch/Student rule fields. Those fields do not repeat the corrected
claims. Its simulation and generalization limitations still matter; a moderate
sample size does not guarantee robustness under every unequal-size/skewed
population combination.

Two candidate-source cautions illustrate why filling every lesson slot is not
a qualification criterion. The [DescTools Stuart-Maxwell documentation](https://andrisignorell.github.io/DescTools/reference/StuartMaxwellTest.html)
uses symmetry language in its Details, while this test's target is marginal
homogeneity. The [NIST Levene page](https://www.itl.nist.gov/div898/handbook/eda/section3/eda35a.htm)
promotes variance pretesting, which is not the advisor's routing policy.
Neither source should be adopted without these qualifications.

Historical source labels and source snapshots remain unchanged. Availability,
bibliographic identity, content review, numerical validation and independent
scientific acceptance are separate evidence levels.
