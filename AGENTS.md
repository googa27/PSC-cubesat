# AGENTS.md — PSC-cubesat

Purpose: Project #24 legacy preservation. This repository is classified as `Hardware / Docs Legacy` with `legacy` profile and Advisory enforcement. Do not make unsupported maturity, security, maintenance, or production-readiness claims.

Canonical docs:
- `README.md` root preservation notice
- `docs/ARCHITECTURE.yaml` machine-readable source of truth
- `docs/ARCHITECTURE.md` rationale and maintained-library/revival notes

Provenance and attribution:
- Origin owner: `googa27`; issue: https://github.com/googa27/PSC-cubesat/issues/1
- Configured canonical origin: https://github.com/googa27/PSC-cubesat.git
- Separate `upstream` remote: none configured locally
- Preserve history, existing public names, authorship, copyright notices, and file contents. Do not delete, rewrite, or hide inherited material in this preservation change.

Safety boundaries:
- License/provenance: No root LICENSE detected; PDFs/DOCX/datasheets may be third-party controlled documents and need rights/provenance review before reuse.
- Data posture: Telemetry/document sources only; schema/provenance, hardware-revision, unit, and operational-safety plans are required before revival.
- Private-data rule: Do not add mission-private telemetry, ground-station credentials, keys, frequencies subject to restriction, unpublished ICDs, or proprietary hardware files.
- Security/hardware warning: Hardware/spacecraft documents and scripts are advisory only; do not use for flight, operations, command generation, or safety-critical work without formal engineering review.

Exact commands:
- Setup: no supported automated setup is declared; treating runtime setup as a revival gate is required.
- Tests: no inherited runtime test suite is claimed; run the architecture checker only.
- Lint/format: no lint/format command is declared.
- Architecture: `python scripts/check_portfolio_architecture.py`

Implementation rules for future work:
- Research upstream/current maintained libraries, standards, datasets, licenses, and security posture before changing runtime code.
- Prefer maintained libraries; custom code must be limited to domain semantics, adapters, composition, or genuinely missing algorithms with oracle/reference tests.
- Avoid invasive refactors of inherited code. Record exact no-growth exceptions and compatibility risks before structural changes.
- Do not introduce generated caches, secrets, private identifiers, restricted data, or fabricated outputs.
- Keep AI-facing contracts deterministic and local. Add Hermes skills for recurring workflows only; plugin/MCP needs stable public contracts, measured multi-client need, least privilege, and separate verification.
- Human/notebook interface: No notebook/package API mandate until revival; preserve hardware semantics and reproducible scripts.
- Core posture: No finance core coupling.

Revival gates:
- Inventory every PDF/DOCX/datasheet/script for rights, source, document revision, and hardware applicability.
- Classify telemetry/command data and exclude mission-private or controlled operational material from public fixtures.
- Add public test vectors with units, frames, endian/encoding, and authoritative expected outputs.
- Obtain formal engineering/security review before using any artifact for operations, command generation, or flight/safety-critical decisions.
- Pin compiler/interpreter dependencies and add deterministic build/check commands before claiming maintained status.

Definition of done for preservation edits:
- README, AGENTS, `docs/ARCHITECTURE.yaml`, `docs/ARCHITECTURE.md`, and tests agree.
- `python scripts/check_portfolio_architecture.py` passes.
- Only advisory governance files are changed unless a separate reviewed revival task authorizes runtime edits.

### AI-assisted change controls

- Treat agent output as untrusted until a human reviews it and executable repository gates verify it. The human author remains accountable.
- Keep agent changes small, single-purpose, and completely reviewable. Generated tests are not a sufficient sole oracle for generated implementation.
- New dependencies require human approval plus package-existence, maintenance, API, license, vulnerability, and typosquat checks; lock reproducibly.
- Security-sensitive code (authentication, cryptography, parsers, serialization, SQL, filesystem, subprocess, network, permissions, or private data) requires dedicated human review.
- Use least privilege: workspace-scoped writes, network/secret access only when approved, no autonomous merge/deploy, and exact command/result provenance.
- Measure AI impact with lead time, review time, CI failures, reverts, escaped defects, and churn; do not infer productivity from self-report.

### Semantic source-tree hierarchy

- Do **not** balance source folders like AVL/B-trees. Package boundaries follow information hiding, cohesion, coupling, public contracts, ownership, and change patterns; naturally heavy-tailed sizes are expected.
- Empty marker packages and speculative folder scaffolds are forbidden unless an exact, dated structural-role exception exists. Keep future plans in architecture/roadmap documents.
- `__init__.py` is a compatibility/public facade only: imports, re-exports, `__all__`, metadata, and bounded lazy hooks. Domain classes and business functions belong in cohesive modules.
- Severe branch concentration is a review trigger, not a command to redistribute files. Fix it only when dependency, churn, ownership, or comprehension evidence shows a bad boundary.

- AI/hierarchy policy: `python3 scripts/check_ai_hierarchy_policy.py`
