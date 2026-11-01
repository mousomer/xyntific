# Governance-only harness adoption

This adoption starts from published Xyntific source `ef34759504741ce8ee3d45b38c84dde480f27267`. `README.md` remains the legacy product description. The immutable workspace policy-kit distribution is 1.4.0 and its runtime is 0.7.0; their release identity and content lock are separate.

The declared command is `make check-governance`. It validates the governance interpreter, manifest/authority references, immutable runtime bytes, compatibility scripts and bounded source sanitation. `scripts/verify.sh` is an alias for this governance-only command. A successful result proves this harness boundary only. It does not certify the C algorithms, CLI compatibility, scientific claims, buildability, model/device behavior, or product acceptance.

The published tree has no `test` or `test-sanitize` Make targets and no automated product test contract. Its existing default `all` target references missing historical `xynRawVid` source. The original Makefile bytes and default/build semantics are preserved; the harness appends an explicit separate target and does not repair, replace, suppress or claim success for `make all`.

No unpublished research/product commits, corrected core implementation, experiment results, new science authority, or private data are included. Product repairs and any new test/acceptance contract require a separate approved scope. Existing inherited safety and environment constraints continue to apply.

Select an existing approved governance interpreter through `GOVERNANCE_PYTHON` or an ignored local `.governance/workspace.local.json`; no environment creation, package installation or network operation is part of verification. Machine choices and collected evidence remain private. The preserved 1.1.3 compatibility scripts/version record are legacy interfaces, not the active runtime identity.


## Continuous verification

`.github/workflows/governance.yml` runs the same `make check-governance` contract
for pull requests and pushes to the main branches or the approved adoption
branch. The hosted Linux job selects Python 3.14.8 explicitly and has a ten-minute
time limit. Checkout and Python setup actions use immutable commit pins; checkout
credentials are not persisted and the workflow token has only `contents: read`.

This lane checks published source and governance metadata only. It uses the
standard-library governance tools profile, so it does not install Python packages
or build the legacy C product. It reads no external corpus or machine overlay,
uses no repository secrets, uploads no artifacts, and performs no data collection
or model download. Existing private verification evidence is not a CI input.
Hosted runner setup is CI infrastructure; local verification continues to use an
already-approved interpreter without creating or changing an environment.

A green CI job proves the declared governance boundary only. The missing legacy
`xynRawVid` dependency and absent product acceptance contract remain open work.
