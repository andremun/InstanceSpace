# Security Policy

This is research software, not a security-critical service, but a
documented reporting policy costs little and is standard practice --
particularly given this toolkit's shared aspiration (with the Python port,
`pyInstanceSpace`) to eventually power the MATILDA web backend.

## Supported versions

Only the latest released version is supported. Fixes land on `master` and
are included in the next release; older releases are not patched
separately.

## Reporting a vulnerability

Please **do not** open a public GitHub issue for a suspected security
vulnerability (e.g. unsafe deserialization of a loaded `.mat`/model file,
path handling in `buildIS`/`exploreIS`, or a supply-chain concern with a
bundled binary such as the LIBSVM MEX files -- see #29).

Instead, report it privately through MATILDA's
[Queries and Feedback](http://matilda.unimelb.edu.au/matilda/contact-us)
page, or via GitHub's
[private vulnerability reporting](https://github.com/andremun/InstanceSpace/security/advisories/new)
if enabled for this repository.

We aim to acknowledge reports within a reasonable time and will credit
the reporter (unless anonymity is requested) once a fix ships. There is
no bug bounty.

## Scope

This toolkit runs locally, on data you provide, inside your own MATLAB
session -- there is no hosted service or network-facing component in this
repository. The main areas worth a security-minded look are:

- Loading untrusted `.mat` model files (`InstanceSpace.load`,
  `ISAmigrateModel`) -- MATLAB's `.mat` format can embed executable
  content, so only load models from sources you trust.
- The bundled LIBSVM MEX binaries, used only for migrating pre-v1.7
  legacy models (see #29 for their provenance status).
