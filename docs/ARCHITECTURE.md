# Architecture — PSC-cubesat

## Project #24 preservation profile

Source of truth: `docs/ARCHITECTURE.yaml`. Tracking issue: https://github.com/googa27/PSC-cubesat/issues/1. Profile: `legacy`; enforcement: Advisory.

This repository is preserved as `Hardware / Docs Legacy`. The governance files are intentionally advisory and additive: they document provenance, risks, and revival gates without refactoring inherited code or claiming active maintenance.

## Archival, supersession, and provenance

- Archival notice: historical/reference preservation only; not production-ready, maintained, secure, or operationally validated.
- Supersession notice: prefer the configured canonical origin repository or maintained libraries for new work.
- Canonical/provenance: historical hardware/docs repository; configured canonical origin is https://github.com/googa27/PSC-cubesat.git.
- Upstream remote: No separate `upstream` remote is configured locally.
- License/provenance warning: No root LICENSE detected; PDFs/DOCX/datasheets may be third-party controlled documents and need rights/provenance review before reuse.
- Security/private-data warning: Hardware/spacecraft documents and scripts are advisory only; do not use for flight, operations, command generation, or safety-critical work without formal engineering review. Do not add mission-private telemetry, ground-station credentials, keys, frequencies subject to restriction, unpublished ICDs, or proprietary hardware files.

## Research-backed defaults

| Decision | Evidence | Repository application |
|---|---|---|
| Agent context | Hermes context files; AGENTS.md convention | Root `AGENTS.md`; progressive detail in this architecture document. |
| AI tool escalation | MCP tools specification | Stable local contracts first; no repo-specific plugin/MCP during preservation. |
| Python source layout | PyPA src-layout guidance | No forced migration for legacy/fork/hardware preservation. |
| Test layout | pytest good practices | Unit/integration/e2e/architecture directories exist; empty suites declare activation triggers. |
| Module budget | Pylint too-many-lines rationale plus AI review locality | 500-line default is a no-growth ratchet where runtime source roots are activated. |
| Evolution | Evolutionary architecture | Revival requires executable fitness functions and explicit exceptions. |
| Data layers | Medallion architecture | Applied only if revived with real data; current posture is advisory. |
| Python protocols | Python data model; NumPy dispatch | Dunders are not decoration; API/protocol redesign waits for revival. |

## Maintained-library decision table

| Capability | Selected route | Alternatives | Boundary / custom-code rule |
|---|---|---|---|
| telemetry/document parsing | python-docx, pdfplumber/pypdf, construct/bitstruct after license/security review | Ad hoc binary/text parsing for new formats | Any parser needs schema, endian/unit, provenance, and test-vector evidence. |
| orbital/attitude math | NumPy/SciPy; Orekit/poliastro only after suitability review | Unvalidated operational algorithms | Research/demo only unless validated against mission-approved references. |
| C++ build/test | CMake/CTest or existing make after revival | Untracked manual compile commands | No operational binary claims without compiler/platform and test-vector evidence. |
| architecture bootstrap | Python standard-library json over JSON-subset YAML | Hand-written YAML parser | Dependency-free advisory gate only. |

## Data, security, and privacy posture

Telemetry/document sources only; schema/provenance, hardware-revision, unit, and operational-safety plans are required before revival.

Do not add mission-private telemetry, ground-station credentials, keys, frequencies subject to restriction, unpublished ICDs, or proprietary hardware files.

Hardware/spacecraft documents and scripts are advisory only; do not use for flight, operations, command generation, or safety-critical work without formal engineering review.

## AI and human interface

- AI interface: Minimal AGENTS documenting hardware/document safety and runnable scripts; no MCP/plugin.
- Human/notebook interface: No notebook/package API mandate until revival; preserve hardware semantics and reproducible scripts.
- Core posture: No finance core coupling.

## Revival gates

- Inventory every PDF/DOCX/datasheet/script for rights, source, document revision, and hardware applicability.
- Classify telemetry/command data and exclude mission-private or controlled operational material from public fixtures.
- Add public test vectors with units, frames, endian/encoding, and authoritative expected outputs.
- Obtain formal engineering/security review before using any artifact for operations, command generation, or flight/safety-critical decisions.
- Pin compiler/interpreter dependencies and add deterministic build/check commands before claiming maintained status.

## Research anchors

- https://hermes-agent.nousresearch.com/docs/user-guide/features/context-files
- https://agents.md/
- https://modelcontextprotocol.io/specification/2025-06-18/server/tools
- https://packaging.python.org/en/latest/discussions/src-layout-vs-flat-layout/
- https://docs.pytest.org/en/stable/explanation/goodpractices.html
- https://docs.python.org/3/reference/datamodel.html
- https://numpy.org/doc/stable/user/basics.dispatch.html
- https://evolutionaryarchitecture.com/precis.html
- https://learn.microsoft.com/en-us/azure/databricks/lakehouse/medallion
