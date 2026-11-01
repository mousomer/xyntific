---
name: governance-context
description: Resolve a governed repository's task context and required evidence before implementation or review. Use for manifest-pack context, route selection, provenance, or verification decisions; follow the consumer's existing dispatcher.
---

Use the consumer's adopted pack to obtain bounded context and a defensible
verification decision. The pack resolves facts; the declared human authorities
determine semantics. This skill guides the workflow and does not enforce gates.

## Resolve context

1. Read the applicable repository dispatcher and inherited instructions. Locate
   its workspace membership, project manifest and locked pack from declarations;
   do not assume a directory layout or a language toolchain. A location without
   a local launcher resolves its declared central fleet record through inherited
   dispatch. Respect read-only/parent ownership and unresolved adoption gaps;
   central preservation does not certify installation or product acceptance.
2. Run the declared governance launcher, normally `./gov check`. Use
   `./gov resolve` and `./gov explain <key>` for facts and their provenance.
   Context validation checks declared adapter body/target/release relationships;
   global discovery must resolve to the active immutable release.
   `./gov explain --task "<description>"` probes representative scenarios.
   Scenario matching has declared overlap rules; a fallback is not proof that
   the task belongs only to its default route. Select every route the actual
   work requires and report ambiguities rather than resolving them by file order.
3. Choose `LOCAL_FIX`, `FEATURE` or `STRUCTURAL_CHANGE` for exploration breadth.
   Escalate discovery when ownership is unclear or crossed, implementations
   disagree, generated authority changes, or packaging/release behaviour changes.
   Execution mode does not remove evidence or authority-update duties.
4. Consume stable resolved facts with provenance. Reopen canonical sources for
   omitted semantics, ambiguity/conflict, boundary or escalation decisions,
   provenance audits, and edits to that authority. Read only the owning design,
   adjacent contracts, implementation and direct tests needed by the task.

## Choose evidence

Compose evidence from all selected routes, the task, actual diff, affected
layers, claims, contracts and risk. A repository change needs nonempty evidence;
record why any declared default is inapplicable. Cross-layer work needs provider,
consumer and integration evidence. Distinguish behavioural, deterministic,
compatibility, named-platform, packaging and human visual claims.

Use the consumer's declared interpreter and runtime profile. Run `./gov doctor`
when execution identity is uncertain; governance-runtime success does not certify
an application toolchain. Run focused checks before the required full gate and
serialize full gates. Report exact commands against the checked tree.

When logging observations matter, read `gov telemetry status --json` for the
last observed append/readback readiness and its time/store identity. An explicit
`gov telemetry selftest --json` is a synthetic private write probe, never product
evidence. Disabled logging creates no observations. Preserve deferred telemetry
gates and do not infer present writability or complete agent activity from status.

A partial mechanical check establishes only its documented boundary. Preserve
remaining human review and platform/device acceptance. A build, metric, lock,
resolved route or skill invocation cannot supply evidence outside that boundary.

## Finish the task

Compare the intended change with its owning design and state acceptance criteria.
Update authorities, documentation and open work when their semantics or validity
change, or when explicit task instructions require it. Preserve unrelated work,
stronger applicable rules and declared external-data constraints.

Report routes, affected layers, evidence, omissions, exact results and remaining
acceptance. Distinguish implementation, automated verification, inspection,
human acceptance and publication. Do not claim absent observations or failed
checks as evidence. Retain explicit unknowns when reconstructing past work.
