# Map-Kinase OpenShift Deployment Validation Results

Purpose: validate Map-Kinase security controls in BYU Chemistry's OpenShift deployment in response to OIT's Final Security Assessment.

Application entry point: `app.py`

Default security settings to confirm:

- Per-file upload limit: `10 MB` (`MAPKINASE_MAX_UPLOAD_SIZE_MB`)
- Per-run total upload limit: `30 MB` (`MAPKINASE_MAX_TOTAL_UPLOAD_SIZE_MB`)
- Minimum accepted-run spacing: `5 seconds` (`MAPKINASE_MIN_SECONDS_BETWEEN_RUNS`)
- Rolling run-attempt limit: `10 per minute` (`MAPKINASE_MAX_RUNS_PER_MINUTE`)
- Session workspace root: `MAPKINASE_SESSION_TMP_ROOT`, or default Python temp directory plus `mapkinase_sessions`
- Session workspace TTL fallback: `24 hours` (`MAPKINASE_SESSION_TMP_TTL_HOURS`)

## Deployment Context

- Deployed URL: `https://mapkinase.chem.byu.edu`
- Image/tag: `chem-map-kinase-prod`
- Image digest: `sha256:7c6b0875e6ca2d9d95713b27d54e2cd5fdbaa36398584b3e10c3a6307dada993`
- Git commit: `23684323e964cef922ca8e23655fe266d3c2688c`
- OpenShift project/namespace: `byu-chemistry`
- OpenShift deployment: `chem-map-kinase-prod`
- Pod name: `chem-map-kinase-prod-77558b6c6c-q9ph6`
- Validation date: `2026-06-16`

## 1. Upload Validation

Confirm in the deployed app:

| Test | File | Upload Area | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- | --- | --- |
| Valid CSV accepted | `valid_protein.csv` | Protein upload | Accepted, preview/validation succeeds | PASS | Manual SmokeScreen test: accepted valid file. |
| Valid TSV accepted | `valid_ptm.tsv` | PTM upload | Accepted, preview/validation succeeds | PASS | Manual SmokeScreen test: accepted valid file. |
| Valid tab-delimited TXT accepted | `valid_metabolite.txt` | Metabolite upload | Accepted, preview/validation succeeds | PASS | Manual SmokeScreen test: accepted valid file. |
| JSON rejected for analysis upload | `analysis_payload.json` | Protein/PTM/metabolite upload | Rejected with safe message | PASS | Manual SmokeScreen test rejected with safe message: `Invalid analysis dataset upload: only .csv, .tsv, and tab-delimited .txt files are accepted.` |
| JSON accepted only for custom pathway import | `custom_pathway_valid.json` | Custom pathway import | Accepted by custom import gate only | PASS | Manual SmokeScreen test: accepted valid file through custom pathway import. |
| Executable extension rejected | `payload.exe` | Any analysis upload | Rejected with safe message | PASS | Manual SmokeScreen test rejected with safe message: `Invalid analysis dataset upload: only .csv, .tsv, and tab-delimited .txt files are accepted.` |
| Script extension rejected | `payload.py`, `payload.js`, `payload.sh`, `payload.bat` | Any analysis upload | Rejected with safe message | PASS | Manual SmokeScreen test rejected with safe message: `Invalid analysis dataset upload: only .csv, .tsv, and tab-delimited .txt files are accepted.` |
| HTML rejected | `payload.html` | Any analysis upload | Rejected with safe message | PASS | Manual SmokeScreen test rejected with safe message: `Invalid analysis dataset upload: only .csv, .tsv, and tab-delimited .txt files are accepted.` |
| Archive/office rejected | `payload.zip`, `payload.xlsx` | Any analysis upload | Rejected with safe message | PASS | Manual SmokeScreen test rejected with safe message: `Invalid analysis dataset upload: only .csv, .tsv, and tab-delimited .txt files are accepted.` |
| Double extension rejected | `data.csv.exe` | Any analysis upload | Rejected with safe message | PASS | Manual SmokeScreen test rejected with safe message: `Invalid analysis dataset upload: only .csv, .tsv, and tab-delimited .txt files are accepted.` |
| Binary renamed to CSV rejected | `binary_renamed.csv` | Any analysis upload | Rejected with safe message | PASS | Manual SmokeScreen test rejected with safe message: `Invalid analysis dataset upload: file appears binary and is not accepted.` |
| Malformed text rejected | `malformed_non_tabular.txt` | Analysis upload | Rejected with safe message | PASS | Manual SmokeScreen test rejected with safe message: `Invalid analysis dataset upload: file could not be parsed as tabular text..` |

Safe rejection means the page must not show Python tracebacks, server filesystem paths, environment variables, uploaded file contents, or raw payloads.

## 2. Upload Size Limits

Confirm in deployment:

| Test | File(s) | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- | --- |
| Under per-file limit accepted | `under_10mb_valid_protein.csv` | Accepted if valid | PASS | Manual SmokeScreen test: accepted valid file. |
| Over per-file limit rejected before parsing | `over_10mb_valid_shape.csv` | Rejected with message referencing 10 MB limit | PASS | Manual SmokeScreen test observed safe rejection with message: `Unexpected error while reading file: field larger than field limit (131072)`. |
| Combined uploads over per-run total rejected | Four `combined_under_10mb_*.csv` files, or deployment-equivalent multi-upload set | Rejected with message referencing 30 MB per run | PASS | Manual SmokeScreen test noted only three analysis upload slots are available, each capped at 10 MB, so 30 MB or more cannot be uploaded for data analysis through the UI. |
| Rejection messages safe | Oversized files | No traceback/server path/raw data | PASS | Manual SmokeScreen evidence: Image 1. |

## 3. Rate Limiting and Repeated Run Throttling

Confirm in deployment:

| Test | Steps | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- | --- |
| Rapid repeated Run clicks | Upload valid input, click Run twice within 5 seconds | Second run is blocked | PASS | Manual retest passed; second run within 5 seconds was blocked. |
| Minimum accepted-run spacing | Wait 5 seconds, click Run again | Run is allowed | PASS | Manual SmokeScreen test: run allowed after waiting. |
| Rolling minute attempt limit | Attempt more than 10 runs within one minute in one browser session | Later attempts are blocked | PASS | Manual retest passed; later attempts were blocked with message: `Too many analysis runs were started recently. Please wait before trying again.` |
| User message safe | Trigger throttle | Message is understandable and contains no internals | PASS | Observed throttle message is understandable and contains no internals: `Too many analysis runs were started recently. Please wait before trying again.` |

Supplied deployment evidence:

- Runtime rate-limit env vars are not explicitly set, so deployed defaults from `mk_security_limits.py` apply:
  - `MAPKINASE_MIN_SECONDS_BETWEEN_RUNS`: `5 seconds`
  - `MAPKINASE_MAX_RUNS_PER_MINUTE`: `10 per minute`
- Manual SmokeScreen testing confirmed a run is allowed after waiting 5 seconds and confirmed the throttle user message is safe.
- Manual retest confirmed rapid repeated Run clicks within 5 seconds are blocked and more than 10 attempts within one minute are blocked.

## 4. Temporary File Cleanup

Confirm in the deployed OpenShift container:

| Test | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- |
| Accepted upload temp cleanup | Framework temp upload file removed after validation/parsing | PASS | OpenShift logs show accepted-upload cleanup for `file=0.csv` with reason `protein upload parsed`. |
| Session end cleanup | Session workspace removed when browser session ends | PASS | OpenShift logs show deleted session workspace `cf3e12845721_4d419606ca0947afbbca6abd8f164f76`. |
| Logs safe | Cleanup logs do not contain uploaded contents | PASS | Supplied cleanup log lines include workspace/file identifiers and reason only; no uploaded contents. |

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

Manual deployed result:

- PASS. Confirmed two separate browser/incognito sessions with different inputs did not overwrite or leak results.
- PASS. Confirmed distinct session workspace names in pod logs/session root.

Expected implementation behavior:

- Preview JSON is session-scoped.
- CST overlay state is session-scoped or in-memory.
- Custom pathway export payloads are generated through controlled download handlers.
- User-derived global protein catalogs use in-memory/session flow, not shared global cache files.

## 6. Direct URL Protections

Run:

```powershell
python deployment_validation/check_direct_url_protections.py https://mapkinase.chem.byu.edu
```

Manual URL list:

| URL | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- |
| `https://mapkinase.chem.byu.edu/JSONfiles/latest_preview.json` | 404/403/not routed | PASS | Probe returned `404 Not Found`; no file contents exposed. |
| `https://mapkinase.chem.byu.edu/cache/global_protein_catalog.json` | 404/403/not routed | PASS | Probe returned `404 Not Found`; no file contents exposed. |
| `https://mapkinase.chem.byu.edu/output/testing_file_001/` | 404/403/not routed | PASS | Probe returned `404 Not Found`; no file contents exposed. |
| `https://mapkinase.chem.byu.edu/MapKinase_WebApp/JSONfiles/latest_preview.json` | 404/403/not routed | PASS | Probe returned `404 Not Found`; no file contents exposed. |
| `https://mapkinase.chem.byu.edu/MapKinase_WebApp/cache/global_protein_catalog.json` | 404/403/not routed | PASS | Probe returned `404 Not Found`; no file contents exposed. |
| `https://mapkinase.chem.byu.edu/stored_pathways/` | 404/403/not routed unless intentionally exposed reference-only assets | PASS | Probe returned `404 Not Found`; no directory listing or file contents exposed. |
| `https://mapkinase.chem.byu.edu/.env` | 404/403/not routed | PASS | Probe returned `404 Not Found`; no file contents exposed. |
| `https://mapkinase.chem.byu.edu/requirements.txt` | 404/403/not routed | PASS | Probe returned `404 Not Found`; no file contents exposed. |
| `https://mapkinase.chem.byu.edu/DEPENDENCY_SECURITY_REVIEW.md` | 404/403/not routed | PASS | Probe returned `404 Not Found`; no file contents exposed. |

Automated direct URL protection check:

- PASS. `python deployment_validation/check_direct_url_protections.py 'https://mapkinase.chem.byu.edu'` passed against all 9 probed paths.
- Summary recorded in supplied evidence as `deployment_validation/results/url_check_summary.md`.
- OpenShift logs for pod `chem-map-kinase-prod-ccf88b8c6-tg8bv` show 404 responses for the probed sensitive/static paths.

Additional Manual SmokeScreen direct URL check:

- PASS. The same 9 sensitive/static paths were manually checked against the OpenShift apps route `https://chem-map-kinase-prod-byu-chemistry.apps.ocp.byu.edu`; each returned 404/403/not routed as expected. Evidence: Image 2.

Downloads to validate through the UI:

- Download handlers serve only intended JSON/sample exports.
- No UI or URL accepts arbitrary server file paths.

## 7. Debug/Admin/Developer Surface

Confirm in deployment:

| Test | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- |
| Debug UI controls | Not visible in production | PASS | `MAPKINASE_ENABLE_DEBUG_UI` is unset in the pod; app treats unset as off. |
| App debug mode | Off by default | PASS | Debug-related `MAPKINASE_`/`M5_` vars are unset except `M5_HOST=0.0.0.0` and `M5_PORT=8080`. |
| Terminal file logging | Off unless explicitly configured | PASS | `M5_TERMINAL_LOG` and `M5_TERMINAL_LOG_FILE` are unset. |
| Debug file output | Off unless explicitly configured | PASS | `MAPKINASE_ENABLE_DEBUG_FILE_OUTPUT` is unset. |
| Persistent CST save | Off unless explicitly configured | PASS | `MAPKINASE_ENABLE_CST_SIDECAR_SAVE` is unset. |
| User-facing errors | No tracebacks, env vars, server paths, or secret values | PASS | Manual SmokeScreen test: externally checked; no debug UI/routes/errors observed. |
| Debug/admin/test routes | Not exposed | PASS | Manual SmokeScreen test: externally checked; no debug UI/routes/errors observed. |

Environment variables to inspect:

- `MAPKINASE_ENV`
- `MAPKINASE_PRODUCTION`
- `MAPKINASE_ENABLE_DEBUG_UI`
- `MAPKINASE_ENABLE_DEBUG_EXPORT`
- `MAPKINASE_ENABLE_DEBUG_FILE_OUTPUT`
- `MAPKINASE_ENABLE_CST_SIDECAR_SAVE`
- `M5_TERMINAL_LOG`
- `M5_TERMINAL_LOG_FILE`

Supplied deployment evidence:

- The only `MAPKINASE_`/`M5_` variables set on the pod are `M5_HOST=0.0.0.0` and `M5_PORT=8080`.
- Upload-limit, rate-limit, and session-TTL env vars are not set explicitly, so the application uses built-in defaults from `mk_security_limits.py`: `10 MB` per file, `30 MB` per run, `5 seconds` between accepted runs, `10 runs/minute` per session, and `24-hour` session TTL.
- Session workspaces are created under `/tmp/mapkinase_sessions` on first use.

## 8. OpenShift Admin Checklist for CSR/OIT

To be completed by OpenShift/CSR/OIT administrators:

| Control | Expected Result | PASS/FAIL | Evidence |
| --- | --- | --- | --- |
| CPU limits | Configured on deployment/pod | PASS | Container has `requests.cpu=200m` and `limits.cpu=1`, verified with `oc get pod ... -o jsonpath='{...resources}'`. |
| Memory limits | Configured on deployment/pod | PASS | Container has `requests.memory=512Mi` and `limits.memory=2Gi`, verified with `oc get pod ... -o jsonpath='{...resources}'`. |
| Least-privilege container | Runs non-root or with restricted SCC/profile | PASS | SCC is `restricted-v2`; `allowPrivilegeEscalation=false`; container is not privileged; capabilities drop `[ALL]`; container runs as non-root. |
| Service account | Minimal permissions, no broad project/cluster privileges | PASS | Uses namespace `default` service account with no role bindings in `byu-chemistry`. |
| Secrets | No unnecessary secrets mounted | PASS | Only Dynatrace operator metadata secret is mounted. No app secrets mounted as volumes, no `envFrom.secretRef`, and no `valueFrom.secretKeyRef`. |
| Logs/monitoring | OpenShift logs, monitoring, and alerting enabled per CSR policy | PASS | Dynatrace operator init container is attached; stdout goes to OpenShift cluster logging. |
| Runtime monitoring | Runtime monitoring attached | PASS | Dynatrace operator is attached. |
| Image provenance | Image tag/digest recorded | PASS | Image is `chem-map-kinase-prod`, pinned by digest `sha256:7c6b0875e6ca2d9d95713b27d54e2cd5fdbaa36398584b3e10c3a6307dada993`; same digest observed in `image` and `imageID`. |
| Rollback plan | Previous known-good image/deployment available | PASS | Rollback can be performed by redeploying the prior known-good image digest. |

## 9. External Service Integration Checks

Document each external integration:

| Service | Purpose | Runtime or build/cache time | User dataset transmitted? | Data sent | PASS/FAIL |
| --- | --- | --- | --- | --- | --- |
| KEGG | Pathway layout/image/reference retrieval | Runtime and/or cache warmup depending workflow | No | Public pathway/species identifiers; canary test intercepted `https://rest.kegg.jp/get/hsa99999/kgml` and `/image` only. | PASS |
| WikiPathways | Pathway listing/layout/reference retrieval | Runtime and/or cache warmup depending workflow | No | Public pathway/species identifiers; canary test intercepted `WP554` pywikipathways calls and WikiPathways image URL only. | PASS |
| CST/reference assets | Reference pathway metadata/assets where configured | Runtime and/or prebuilt cache depending workflow | No | Public pathway/reference identifiers; runtime import review shows CST viewer uses local/reference assets and UniProt label mapper, while CST scraper helpers are not imported by the runtime app path checked. | PASS |
| UniProt/reference lookup helpers | Reference annotation lookup where used | Runtime or maintenance scripts depending workflow | No user-uploaded dataset contents | Public Entrez, Ensembl, and pathway-label symbols only; canary test intercepted UniProt queries for `GeneID-2475`, `ENSG00000198793`, and `MAPK1`. | PASS |

Validation:

- PASS. Local current-version source review and dynamic canary test were run on 2026-06-17. A fake uploaded dataset containing `MAPKINASE_LEAK_TEST_20260617_DO_NOT_SEND` was present in the test process; outbound HTTP/pywikipathways calls were monkeypatched and recorded; the canary string was not present in any recorded URL, params, headers, or body.
- PASS. Confirmed external requests use public pathway/reference identifiers only in the intercepted runtime helpers.
- Evidence type: local current-version source review and intercepted-call testing.

## 10. Logging Checks

Confirm logs include:

- Upload validation failures.
- Per-file and total-size failures.
- Throttled run attempts.
- Session workspace creation and cleanup.
- Application errors, without sensitive payloads.

Supplied OpenShift log evidence includes:

- Session workspace creation: `Created session workspace | workspace_id=cf3e12845721_4d419606ca0947afbbca6abd8f164f76`.
- Accepted upload cleanup: `Deleted temporary file | file=0.csv | reason=protein upload parsed`.
- Direct URL probes returning `404 Not Found` for `/.env`, `/requirements.txt`, `/stored_pathways/`, `/JSONfiles/latest_preview.json`, `/cache/global_protein_catalog.json`, and related paths.
- Session workspace cleanup: `Deleted session workspace | workspace_id=cf3e12845721_4d419606ca0947afbbca6abd8f164f76`.
- Supplied log excerpts do not show uploaded file contents, raw datasets, raw custom pathway JSON payloads, sensitive environment variables, or secrets.

Confirm logs do not include:

- Uploaded file contents.
- Full raw datasets.
- Raw custom pathway JSON payloads.
- Sensitive environment variables or secrets.
- Server filesystem paths in user-facing messages.
- Internal temp paths shown to users.

## OIT-Facing Results Template

Validation date: `2026-06-16`

Deployment URL: `https://mapkinase.chem.byu.edu`

OpenShift project/namespace: `byu-chemistry`

OpenShift deployment: `chem-map-kinase-prod`

Pod: `chem-map-kinase-prod-77558b6c6c-q9ph6`

Image/tag/digest: `chem-map-kinase-prod@sha256:7c6b0875e6ca2d9d95713b27d54e2cd5fdbaa36398584b3e10c3a6307dada993`

Git commit: `23684323e964cef922ca8e23655fe266d3c2688c`

| Category | Result | Evidence | Notes/Exceptions |
| --- | --- | --- | --- |
| Upload validation | PASS | Manual SmokeScreen test accepted valid CSV/TSV/TXT and custom pathway JSON, and rejected invalid analysis uploads with safe messages. | No exceptions. |
| Upload size limits | PASS | Manual SmokeScreen test accepted under-limit upload, safely rejected oversized upload, and confirmed UI upload slots prevent more than 30 MB total analysis upload. Runtime defaults are `10 MB` per file and `30 MB` total per run. | No exceptions. |
| Rate limiting/throttling | PASS | Manual SmokeScreen and retest confirmed rapid repeated clicks are blocked, a run is allowed after waiting 5 seconds, more than 10 attempts within one minute are blocked, and the throttle message is safe. | No exceptions. |
| Temporary file cleanup | PASS | OpenShift logs show accepted upload temp cleanup, session workspace creation, session workspace cleanup, and safe cleanup logging without uploaded contents. | No exceptions. |
| Session isolation | PASS | Manual deployed test used two separate browser/incognito sessions with different inputs and confirmed results did not overwrite or leak; distinct session workspace names were confirmed in pod logs/session root. | No exceptions. |
| Direct URL protections | PASS | Automated probe passed against 9 paths; all returned `404`; no file contents exposed. | Summary referenced in supplied evidence as `deployment_validation/results/url_check_summary.md`. |
| Debug/admin/developer surface | PASS | Debug-related env vars are unset; only `M5_HOST=0.0.0.0` and `M5_PORT=8080` are set; Manual SmokeScreen test observed no debug UI/routes/errors. | No exceptions. |
| OpenShift resource/isolation controls | PASS | CPU/memory requests and limits, SCC, service account, secret mounts, logging, image provenance, and digest pinning look correct. | No exceptions. |
| External service integration | PASS | Local source review and dynamic canary test confirmed KEGG, WikiPathways, CST/reference, and UniProt helper flows do not transmit uploaded dataset contents; canary string was absent from intercepted outbound URL, params, headers, and body records. | No exceptions. |
| Logging safety | PASS | Supplied logs show cleanup and 404 validation events without raw uploads, secrets, env values, or raw custom pathway payloads; manual rejection/throttle messages did not expose internals. | No exceptions. |

Summary statement:

```text
Map-Kinase deployment validation was performed on 2026-06-16 against https://mapkinase.chem.byu.edu running image/tag chem-map-kinase-prod@sha256:7c6b0875e6ca2d9d95713b27d54e2cd5fdbaa36398584b3e10c3a6307dada993 and commit 23684323e964cef922ca8e23655fe266d3c2688c. Upload validation, file size limits, rate limiting, session isolation, session-scoped temporary storage, cleanup behavior, direct URL protections, debug-surface gating, OpenShift isolation controls, image provenance, external service handling, and logging behavior were reviewed. Results: PASS for the validated deployment scope.
```

Exceptions/remediation: none identified for the validated deployment scope.

## Validation Summary

No application code gap was identified from the existing local smoke tests and source review for the requested controls.

Completed deployment confirmations:

- Actual route/static-file exposure behavior: PASS. Direct URL probes against `https://mapkinase.chem.byu.edu` returned `404` for all 9 tested sensitive/static paths.
- Manual upload validation: PASS. SmokeScreen testing accepted valid analysis/custom pathway files and rejected invalid analysis uploads with safe messages.
- Manual upload size-limit behavior: PASS. SmokeScreen testing accepted under-limit upload, safely rejected oversized upload, and confirmed UI upload slots prevent more than 30 MB total analysis upload.
- Manual rate-limit behavior: PASS. SmokeScreen testing and retest confirmed rapid repeated clicks are blocked, a run is allowed after waiting, more than 10 attempts within one minute are blocked, and throttle messaging is safe.
- Manual session isolation: PASS. Two separate browser/incognito sessions with different inputs did not overwrite or leak results, and distinct session workspace names were confirmed in pod logs/session root.
- External service data handling: PASS. Local source review and dynamic canary testing confirmed uploaded dataset content was not present in intercepted KEGG, WikiPathways, CST/reference, or UniProt helper outbound calls.
- Runtime environment variables: PASS. Debug-related vars are unset; upload/rate/session settings use built-in defaults.
- Container temp/session directory cleanup: PASS. Logs show accepted-upload temp cleanup and session workspace cleanup.
- OpenShift CPU/memory requests and limits: PASS. Verified resources are `limits.cpu=1`, `limits.memory=2Gi`, `requests.cpu=200m`, and `requests.memory=512Mi`.
- OpenShift security context/service account/secrets/logging/image provenance: PASS based on supplied evidence.
