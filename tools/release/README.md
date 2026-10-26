# Installed release preflight

Run against a freshly installed source archive, with its library prepended to
`R_LIBS` and recorded in the required `INFERGUIDR_EXPECTED_LIBRARY` variable:

```sh
Rscript --vanilla tools/release/check-package.R /path/to/new/evidence
```

The preflight checks exported help aliases, namespace bindings, the installed
licence, startup state preservation, core calculations and portable JSON/text
reports. It retains a JSON receipt and session information even on a failed
check. Its output directory must be new or empty.

For the optional-dependency absence profile, install only inferguidR, digest
and jsonlite in an isolated library, hide site and user libraries, and append
`--minimal`. Recommended packages shipped with R can remain available. The
profile requires Shiny, ggplot2, jsonvalidate, pROC, Exact, exact2x2 and WebPower
to be absent; it verifies typed refusals for unavailable display functionality
and successful core/report operations. A full profile and a minimal profile
must use the same source archive.

These checks do not certify numerical correctness, all filesystem/network
behavior, graphical interfaces, hosted platforms or readiness for submission.
Record external network tracing or offline-browser evidence separately.

Audit the same source archive without extraction or execution:

```sh
python3 tools/release/audit-source.py inferguidR_VERSION.tar.gz /path/to/new/archive-evidence
```

The inventory records every regular file's SHA-256 and size. The audit rejects
unsafe paths, links/special members, duplicate names, retained private or
development directories, machine-specific references, and unexpected binaries
in this pure-R package. It also checks the package identity and licence files.
This is distribution hygiene; a clean result does not replace a review of
content ownership and third-party licences.

The optional online resource audit performs HTTP reads only when explicitly
invoked; it is never run on package load or during ordinary package checks:

```sh
python3 tools/release/audit-resources.py inst/metadata/learning-resources.json /path/to/new/resource-evidence
```

It records HTTP status, redirects, page hashes/titles and exact fragment
presence using two concurrent requests. Unavailable pages remain failures.
Successful retrieval does not verify the page's scientific claims or the
catalog's semantic mapping to a method.
