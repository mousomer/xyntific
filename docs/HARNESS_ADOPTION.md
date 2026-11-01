# Governance-only harness adoption

This adoption starts from published Xyntific source `ef34759504741ce8ee3d45b38c84dde480f27267`. `README.md` remains the legacy product description. The immutable workspace policy-kit distribution is 1.4.0 and its runtime is 0.7.0; their release identity and content lock are separate.

The declared command is `make check-governance`. It validates the governance interpreter, manifest/authority references, immutable runtime bytes, compatibility scripts and bounded source sanitation. `scripts/verify.sh` is an alias for this governance-only command. A successful result proves this harness boundary only. It does not certify the C algorithms, CLI compatibility, scientific claims, buildability, model/device behavior, or product acceptance.

The published tree has no `test` or `test-sanitize` Make targets and no automated product test contract. Its existing default `all` target references missing historical `xynRawVid` source. The original Makefile bytes and default/build semantics are preserved; the harness appends an explicit separate target and does not repair, replace, suppress or claim success for `make all`.

No unpublished research/product commits, corrected core implementation, experiment results, new science authority, or private data are included. Product repairs and any new test/acceptance contract require a separate approved scope. Existing inherited safety and environment constraints continue to apply.

Select an existing approved governance interpreter through `GOVERNANCE_PYTHON` or an ignored local `.governance/workspace.local.json`; no environment creation, package installation or network operation is part of verification. Machine choices and collected evidence remain private. The preserved 1.1.3 compatibility scripts/version record are legacy interfaces, not the active runtime identity.
