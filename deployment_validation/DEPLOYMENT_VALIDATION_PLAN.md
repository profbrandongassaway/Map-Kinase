# Map-Kinase OpenShift Deployment Validation Plan

Purpose: validate Map-Kinase security controls in BYU Chemistry's OpenShift deployment in response to OIT's Final Security Assessment.

Application entry point: `app.py`

Default security settings to confirm:

- Per-file upload limit: `10 MB` (`MAPKINASE_MAX_UPLOAD_SIZE_MB`)
- Per-run total upload limit: `30 MB` (`MAPKINASE_MAX_TOTAL_UPLOAD_SIZE_MB`)
- Minimum accepted-run spacing: `5 seconds` (`MAPKINASE_MIN_SECONDS_BETWEEN_RUNS`)
- Rolling run-attempt limit: `10 per minute` (`MAPKINASE_MAX_RUNS_PER_MINUTE`)
- Session workspace root: `MAPKINASE_SESSION_TMP_ROOT`, or default Python temp directory plus `mapkinase_sessions`
- Session workspace TTL fallback: `24 hours` (`MAPKINASE_SESSION_TMP_TTL_HOURS`)

## Pre-Validation Setup

1. Record deployed URL:
   - `https://DOMAIN_HERE`

2. Confirm deployment environment:
   - Python version:
   - Image/tag:
   - Git commit:
   - OpenShift project/namespace:
   - Pod name:
   - Validation date:

3. Generate safe synthetic upload files:

   ```powershell
   python deployment_validation/generate_synthetic_uploads.py
   ```

4. Run deployed direct-URL probes:

   ```powershell
   python deployment_validation/check_direct_url_protections.py https://DOMAIN_HERE
   ```

5. Use application logs during testing. In OpenShift, typical commands are:

   ```bash
   oc get pods
   oc logs <pod-name> --since=30m
   oc exec <pod-name> -- printenv | grep -E 'MAPKINASE|M5_'
   oc exec <pod-name> -- python - <<'PY'
   import tempfile
   from pathlib import Path
   print(Path(tempfile.gettempdir()).resolve())
   print((Path(tempfile.gettempdir()) / "mapkinase_sessions").resolve())
   PY
   ```

## 1. Upload Validation

Confirm in the deployed app:

| Test | File | Upload Area | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- | --- | --- |
| Valid CSV accepted | `valid_protein.csv` | Protein upload | Accepted, preview/validation succeeds |  |  |
| Valid TSV accepted | `valid_ptm.tsv` | PTM upload | Accepted, preview/validation succeeds |  |  |
| Valid tab-delimited TXT accepted | `valid_metabolite.txt` | Metabolite upload | Accepted, preview/validation succeeds |  |  |
| JSON rejected for analysis upload | `analysis_payload.json` | Protein/PTM/metabolite upload | Rejected with safe message |  |  |
| JSON accepted only for custom pathway import | `custom_pathway_valid.json` | Custom pathway import | Accepted by custom import gate only |  |  |
| Executable extension rejected | `payload.exe` | Any analysis upload | Rejected with safe message |  |  |
| Script extension rejected | `payload.py`, `payload.js`, `payload.sh`, `payload.bat` | Any analysis upload | Rejected with safe message |  |  |
| HTML rejected | `payload.html` | Any analysis upload | Rejected with safe message |  |  |
| Archive/office rejected | `payload.zip`, `payload.xlsx` | Any analysis upload | Rejected with safe message |  |  |
| Double extension rejected | `data.csv.exe` | Any analysis upload | Rejected with safe message |  |  |
| Binary renamed to CSV rejected | `binary_renamed.csv` | Any analysis upload | Rejected with safe message |  |  |
| Malformed text rejected | `malformed_non_tabular.txt` | Analysis upload | Rejected with safe message |  |  |
| Path traversal filename rejected | Local smoke test or controlled multipart tool | Upload filename `../escape.csv` | Rejected with safe message |  |  |

Safe rejection means the page must not show Python tracebacks, server filesystem paths, environment variables, uploaded file contents, or raw payloads.

Log evidence to collect:

- `mapkinase.upload_validation` rejection log lines.
- Confirm logs include session/role/filename/size/reason only.
- Confirm logs do not include uploaded file contents.

## 2. Upload Size Limits

Confirm in deployment:

| Test | File(s) | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- | --- |
| Under per-file limit accepted | `under_10mb_valid_protein.csv` | Accepted if valid |  |  |
| Over per-file limit rejected before parsing | `over_10mb_valid_shape.csv` | Rejected with message referencing 10 MB limit |  |  |
| Combined uploads over per-run total rejected | Four `combined_under_10mb_*.csv` files, or deployment-equivalent multi-upload set | Rejected with message referencing 30 MB per run |  |  |
| Rejection messages safe | Oversized files | No traceback/server path/raw data |  |  |
| Logs safe | Oversized files | Rejection logged without file contents |  |  |

Evidence to collect:

- Browser screenshot of safe rejection message.
- OpenShift log line showing size rejection.
- Confirm no preview/session artifact was generated for rejected oversized/throttled attempts.

## 3. Rate Limiting and Repeated Run Throttling

Confirm in deployment:

| Test | Steps | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- | --- |
| Rapid repeated Run clicks | Upload valid input, click Run twice within 5 seconds | Second run is blocked |  |  |
| Minimum accepted-run spacing | Wait 5 seconds, click Run again | Run is allowed |  |  |
| Rolling minute attempt limit | Attempt more than 10 runs within one minute in one browser session | Later attempts are blocked |  |  |
| Throttled attempt has no expensive output | Inspect logs/session workspace before and after throttled attempt | No new parsing/rendering/generated output |  |  |
| User message safe | Trigger throttle | Message is understandable and contains no internals |  |  |

Log evidence:

- `mapkinase.run_guard` entries with session, reason, total size, and attempt count.
- No raw uploaded data in logs.

## 4. Temporary File Cleanup

Confirm in the deployed OpenShift container:

| Test | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- |
| Accepted upload temp cleanup | Framework temp upload file removed after validation/parsing |  |  |
| Rejected upload temp cleanup | Rejected temp file removed |  |  |
| Exception cleanup | Validation/parsing exception still removes temp file |  |  |
| Session end cleanup | Session workspace removed when browser session ends |  |  |
| TTL fallback cleanup | Workspaces older than 24 hours are removed on app startup/import cleanup path |  |  |
| Logs safe | Cleanup logs do not contain uploaded contents |  |  |

Directories to inspect:

- Shiny/framework upload temp location, if visible in container temp directory.
- `MAPKINASE_SESSION_TMP_ROOT`, if set.
- Default session workspace root: Python `tempfile.gettempdir()/mapkinase_sessions`.
- App-local paths that should not receive user/session output:
  - `MapKinase_WebApp/JSONfiles/`
  - `MapKinase_WebApp/cache/`
  - `output/`

Suggested commands:

```bash
oc exec <pod-name> -- printenv MAPKINASE_SESSION_TMP_ROOT
oc exec <pod-name> -- python - <<'PY'
import tempfile
from pathlib import Path
root = Path(tempfile.gettempdir()) / "mapkinase_sessions"
print(root)
for path in sorted(root.glob("*"))[:25]:
    print(path, path.stat().st_mtime)
PY
```

## 5. Session Isolation

Manual test:

1. Open Map-Kinase in Browser A.
2. Open Map-Kinase in Browser B using a separate browser profile or incognito session.
3. Upload different valid datasets in each session.
4. Run analysis in Browser A.
5. Run a different analysis in Browser B.
6. Confirm Browser A results do not change when Browser B runs.
7. Confirm Browser B cannot see Browser A generated outputs.
8. Inspect logs/container session root and confirm distinct session workspace names.
9. Confirm shared global paths are not modified by user runs:
   - `MapKinase_WebApp/JSONfiles/latest_preview.json` should not exist in the deployment image or be rewritten.
   - `MapKinase_WebApp/cache/global_protein_catalog.json` should not exist in the deployment image or be rewritten.

Expected implementation behavior:

- Preview JSON is session-scoped.
- CST overlay state is session-scoped or in-memory.
- Custom pathway export payloads are generated through controlled download handlers.
- User-derived global protein catalogs use in-memory/session flow, not shared global cache files.

## 6. Direct URL Protections

Run:

```powershell
python deployment_validation/check_direct_url_protections.py https://DOMAIN_HERE
```

Manual URL list:

| URL | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- |
| `https://DOMAIN_HERE/JSONfiles/latest_preview.json` | 404/403/not routed |  |  |
| `https://DOMAIN_HERE/cache/global_protein_catalog.json` | 404/403/not routed |  |  |
| `https://DOMAIN_HERE/output/testing_file_001/` | 404/403/not routed |  |  |
| `https://DOMAIN_HERE/MapKinase_WebApp/JSONfiles/latest_preview.json` | 404/403/not routed |  |  |
| `https://DOMAIN_HERE/MapKinase_WebApp/cache/global_protein_catalog.json` | 404/403/not routed |  |  |
| `https://DOMAIN_HERE/stored_pathways/` | 404/403/not routed unless intentionally exposed reference-only assets |  |  |
| `https://DOMAIN_HERE/.env` | 404/403/not routed |  |  |
| `https://DOMAIN_HERE/requirements.txt` | 404/403/not routed |  |  |
| `https://DOMAIN_HERE/DEPENDENCY_SECURITY_REVIEW.md` | 404/403/not routed |  |  |

Downloads to validate through the UI:

- Download handlers serve only intended JSON/sample exports.
- No UI or URL accepts arbitrary server file paths.

## 7. Debug/Admin/Developer Surface

Confirm in deployment:

| Test | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- |
| Debug UI controls | Not visible in production |  |  |
| App debug mode | Off by default |  |  |
| Terminal file logging | Off unless explicitly configured |  |  |
| Debug file output | Off unless explicitly configured |  |  |
| Persistent CST save | Off unless explicitly configured |  |  |
| User-facing errors | No tracebacks, env vars, server paths, or secret values |  |  |
| Debug/admin/test routes | Not exposed |  |  |

Environment variables to inspect:

- `MAPKINASE_ENV`
- `MAPKINASE_PRODUCTION`
- `MAPKINASE_ENABLE_DEBUG_UI`
- `MAPKINASE_ENABLE_DEBUG_EXPORT`
- `MAPKINASE_ENABLE_DEBUG_FILE_OUTPUT`
- `MAPKINASE_ENABLE_CST_SIDECAR_SAVE`
- `M5_TERMINAL_LOG`
- `M5_TERMINAL_LOG_FILE`

## 8. OpenShift Admin Checklist for CSR/OIT

To be completed by OpenShift/CSR/OIT administrators:

| Control | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- |
| CPU limits | Configured on deployment/pod |  |  |
| Memory limits | Configured on deployment/pod |  |  |
| Least-privilege container | Runs non-root or with restricted SCC/profile |  |  |
| Service account | Minimal permissions, no broad project/cluster privileges |  |  |
| Filesystem permissions | Write access limited to temp/session directories as needed |  |  |
| Secrets | No unnecessary secrets mounted |  |  |
| Data access | No access to unrelated BYU systems/data |  |  |
| Network egress | Limited or explicitly documented |  |  |
| Logs/monitoring | OpenShift logs, monitoring, and alerting enabled per CSR policy |  |  |
| Runtime monitoring | Falcon or equivalent coverage clarified |  |  |
| Image provenance | Image tag/digest recorded |  |  |
| Rollback plan | Previous known-good image/deployment available |  |  |

Suggested commands:

```bash
oc get deploy <deployment-name> -o yaml
oc get pod <pod-name> -o yaml
oc describe pod <pod-name>
oc adm policy who-can get pods -n <namespace>
oc get rolebindings -n <namespace>
oc get networkpolicy -n <namespace>
```

## 9. External Service Integration Checks

Document each external integration:

| Service | Purpose | Runtime or build/cache time | User dataset transmitted? | Data sent | PASS/FAIL |
| --- | --- | --- | --- | --- | --- |
| KEGG | Pathway layout/image/reference retrieval | Runtime and/or cache warmup depending workflow | No | Public pathway/species identifiers |  |
| WikiPathways | Pathway listing/layout/reference retrieval | Runtime and/or cache warmup depending workflow | No | Public pathway/species identifiers |  |
| CST/reference assets | Reference pathway metadata/assets where configured | Runtime and/or prebuilt cache depending workflow | No | Public pathway/reference identifiers |  |
| UniProt/reference lookup helpers | Reference annotation lookup where used | Runtime or maintenance scripts depending workflow | No user-uploaded dataset contents | Public IDs only |  |

Validation:

- Review logs during a user analysis and confirm no uploaded datasets are sent to external services.
- Confirm external requests use public pathway/reference identifiers only.
- If egress is unrestricted, document why or restrict per CSR/OIT policy.

## 10. Logging Checks

Confirm logs include:

- Upload validation failures.
- Per-file and total-size failures.
- Throttled run attempts.
- Session workspace creation and cleanup.
- Application errors, without sensitive payloads.

Confirm logs do not include:

- Uploaded file contents.
- Full raw datasets.
- Raw custom pathway JSON payloads.
- Sensitive environment variables or secrets.
- Server filesystem paths in user-facing messages.
- Internal temp paths shown to users.

OpenShift/CSR should confirm:

- Log retention and rotation policy.
- Alerting/monitoring integration.
- Access controls for logs that may include filenames or session identifiers.

## Manual Deployed Website Test Steps

1. Generate test files with `python deployment_validation/generate_synthetic_uploads.py`.
2. In the deployed site, upload `valid_protein.csv`, `valid_ptm.tsv`, and optionally `valid_metabolite.txt`; confirm valid files are accepted.
3. Upload each blocked-extension file to analysis uploads; confirm rejection and safe messages.
4. Upload `analysis_payload.json` to analysis uploads; confirm rejection.
5. Upload `custom_pathway_valid.json` only through the custom pathway import flow; confirm JSON is accepted there and rejected elsewhere.
6. Upload `binary_renamed.csv` and `malformed_non_tabular.txt`; confirm rejection.
7. Upload `over_10mb_valid_shape.csv`; confirm rejection before parsing.
8. Use valid inputs and click Run twice within 5 seconds; confirm throttling.
9. Attempt more than 10 Run clicks within one minute in a single session; confirm throttling.
10. Open two separate browser/incognito sessions, run different inputs, and confirm results do not overwrite or leak.
11. Run `python deployment_validation/check_direct_url_protections.py https://DOMAIN_HERE`.
12. Review OpenShift logs for validation, size-limit, throttling, and cleanup events without uploaded contents.

## OIT-Facing Results Template

Validation date:

Deployment URL:

Image/tag/digest:

Git commit:

Validator:

| Category | Result | Evidence | Notes/Exceptions |
| --- | --- | --- | --- |
| Upload validation | PASS/FAIL |  |  |
| Upload size limits | PASS/FAIL |  |  |
| Rate limiting/throttling | PASS/FAIL |  |  |
| Temporary file cleanup | PASS/FAIL |  |  |
| Session isolation | PASS/FAIL |  |  |
| Direct URL protections | PASS/FAIL |  |  |
| Debug/admin/developer surface | PASS/FAIL |  |  |
| OpenShift resource/isolation controls | PASS/FAIL |  |  |
| External service integration | PASS/FAIL |  |  |
| Logging safety | PASS/FAIL |  |  |

Summary statement:

```text
Map-Kinase deployment validation was performed on DATE against DEPLOYMENT_URL running image/tag IMAGE_TAG and commit COMMIT. Upload validation, file size limits, rate limiting, session-scoped temporary storage, cleanup behavior, direct URL protections, debug-surface gating, external service handling, and logging behavior were reviewed. Results: PASS/FAIL. Exceptions and remediation actions are listed below.
```

Exceptions/remediation:

| Finding | Risk | Compensating Controls | Owner | Target Date |
| --- | --- | --- | --- | --- |
|  |  |  |  |  |

## Current Known Gaps

No application code gap was identified from the existing local smoke tests and source review for the requested controls. Deployment validation still requires OpenShift-specific confirmation for:

- Actual route/static-file exposure behavior.
- Runtime environment variables.
- Container temp/session directory cleanup in the deployed pod.
- CPU/memory/security-context/service-account controls.
- Egress policy and runtime monitoring coverage.

The browser-only manual test cannot reliably submit a path traversal filename because browsers normally strip path components from file uploads. Keep the existing local smoke test evidence for path traversal filename rejection, or have OIT use a controlled multipart testing tool if they know the Shiny upload endpoint for the deployed session.
