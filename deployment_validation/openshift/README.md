# Map-Kinase OpenShift Validation Toolkit

This directory contains scripts that turn the OpenShift-specific portions of
[`DEPLOYMENT_VALIDATION_PLAN.md`](../DEPLOYMENT_VALIDATION_PLAN.md) into a
roughly 10-minute, copy-paste runbook. The scripts only collect and summarize
information from the cluster — no application code is touched, no secret values
are printed, and no uploaded file contents are read.

## What gets produced

All outputs land under `deployment_validation/results/` (gitignored):

| File | Source script | What it covers in the plan |
| --- | --- | --- |
| `openshift_evidence/summary.txt` | `collect_openshift_evidence.sh` | Pre-Validation Setup block (URL, image, digest, commit, namespace, pod) |
| `openshift_evidence/deployment.yaml` | `collect_openshift_evidence.sh` | Section 8 raw evidence |
| `openshift_evidence/pod.yaml` | `collect_openshift_evidence.sh` | Section 8 raw evidence |
| `openshift_evidence/pod_describe.txt` | `collect_openshift_evidence.sh` | Section 8 raw evidence |
| `openshift_evidence/env_mapkinase.txt` | `collect_openshift_evidence.sh` | Section 2/3/7 runtime env |
| `openshift_evidence/tempdirs.txt` | `collect_openshift_evidence.sh` | Section 4 temp/session paths |
| `openshift_evidence/serviceaccount.txt` | `collect_openshift_evidence.sh` | Section 8 SA |
| `openshift_evidence/scc.txt` | `collect_openshift_evidence.sh` | Section 8 SCC + securityContext |
| `openshift_evidence/rolebindings.txt` | `collect_openshift_evidence.sh` | Section 8 SA permissions |
| `openshift_evidence/networkpolicy.yaml` | `collect_openshift_evidence.sh` | Section 8 / 9 egress |
| `openshift_evidence/secrets_mounted.txt` | `collect_openshift_evidence.sh` | Section 8 secrets (names only) |
| `openshift_evidence/recent_logs_sanitized.txt` | `collect_openshift_evidence.sh` | Section 10 log evidence |
| `runtime_settings_report.md` | `report_runtime_settings.py` | Sections 2, 3, 4 limits + Section 7 debug surface |
| `openshift_controls_report.md` | `summarize_openshift_controls.py` | Section 8 OpenShift admin checklist |
| `url_check_summary.md` | `check_direct_url_protections.py` | Section 6 direct URL protections |
| `url_check_results.json` | `check_direct_url_protections.py` | Raw probe data for the audit trail |

## Prerequisites

- `oc` CLI on PATH, with `oc login` completed for the OpenShift cluster hosting
  Map-Kinase.
- `oc project <namespace>` set to the Map-Kinase project before running the
  collector.
- Python 3.10+ with `pyyaml` (only required for `summarize_openshift_controls.py`):
  ```bash
  pip install pyyaml
  ```

## 10-minute runbook

Run from the repository root:

```bash
# 1. Log in and select the project
oc login <cluster-url>
oc project <mapkinase-namespace>

# 2. Collect evidence (writes to deployment_validation/results/openshift_evidence/)
bash deployment_validation/openshift/collect_openshift_evidence.sh

# 3. Generate the runtime-settings PASS/REVIEW/FAIL report
python deployment_validation/openshift/report_runtime_settings.py

# 4. Generate the OpenShift controls PASS/REVIEW/FAIL report
python deployment_validation/openshift/summarize_openshift_controls.py

# 5. Probe direct URL protections against the deployed site
#    (URL is printed at the end of step 2; or `oc get route` to find it)
python deployment_validation/check_direct_url_protections.py https://<deployed-url>
```

That's it — every OpenShift row in the OIT-Facing Results Template at the bottom
of `DEPLOYMENT_VALIDATION_PLAN.md` is now backed by a file you can attach or
paste from.

## If the collector can't find the pod

By default the collector tries a list of common Map-Kinase labels:

```
app=mapkinase
app=map-kinase
app=mapkinase-app
deployment=mapkinase
deployment=map-kinase
deploymentconfig=mapkinase
deploymentconfig=map-kinase
name=mapkinase
name=map-kinase
```

If none of those match, it prints the list and exits. To find the real label:

```bash
oc get pods --show-labels
```

Then re-run with an override:

```bash
# Pin to a specific label
MAPKINASE_POD_LABEL=app=map-kinase \
    bash deployment_validation/openshift/collect_openshift_evidence.sh

# Or pin to a specific pod
MAPKINASE_POD_NAME=mapkinase-7c9f88b4d-xyz12 \
    bash deployment_validation/openshift/collect_openshift_evidence.sh

# Specify the Deployment if ownerReferences resolution fails
MAPKINASE_DEPLOYMENT=mapkinase \
    bash deployment_validation/openshift/collect_openshift_evidence.sh

# Widen the log capture window (default 30m)
MAPKINASE_LOG_WINDOW=2h \
    bash deployment_validation/openshift/collect_openshift_evidence.sh
```

**Tip:** when you pass the deployed URL to the URL-protection probe, put it in
single quotes — bare `<...>` placeholders are interpreted by bash as redirects:

```bash
python deployment_validation/check_direct_url_protections.py 'https://mapkinase.example.byu.edu'
```

## Mapping reports to the OIT-Facing Results Template

The template at the bottom of `DEPLOYMENT_VALIDATION_PLAN.md` has these rows.
Use the generated artifacts as follows:

| Template row | Drive PASS/FAIL from | Evidence to attach/cite |
| --- | --- | --- |
| Upload validation | App-level smoke tests + `recent_logs_sanitized.txt` | Browser screenshots, log file |
| Upload size limits | App-level smoke tests + `runtime_settings_report.md` (rows for `MAPKINASE_MAX_UPLOAD_SIZE_MB`, `MAPKINASE_MAX_TOTAL_UPLOAD_SIZE_MB`) | Browser screenshots, runtime report |
| Rate limiting/throttling | App-level smoke tests + `runtime_settings_report.md` (rows for `MAPKINASE_MIN_SECONDS_BETWEEN_RUNS`, `MAPKINASE_MAX_RUNS_PER_MINUTE`) + `recent_logs_sanitized.txt` | Browser screenshots, runtime report, log file |
| Temporary file cleanup | `tempdirs.txt` + `runtime_settings_report.md` (TTL row) + `recent_logs_sanitized.txt` | Tempdir listing, runtime report, log file |
| Session isolation | App-level browser test (two profiles) + `tempdirs.txt` | Browser observations, tempdir listing |
| Direct URL protections | `url_check_summary.md` | URL check summary |
| Debug/admin/developer surface | `runtime_settings_report.md` (debug flags + production env) | Runtime report |
| OpenShift resource/isolation controls | `openshift_controls_report.md` | Controls report + raw `pod.yaml` |
| External service integration | `recent_logs_sanitized.txt` + `networkpolicy.yaml` | Log file, network policy |
| Logging safety | `recent_logs_sanitized.txt` | Log file |

Pre-Validation Setup block at the top of the plan: copy directly from
`openshift_evidence/summary.txt`.

## Privacy and safety notes

- `env_mapkinase.txt` only contains variables matching `^(MAPKINASE_|M5_)`. The
  collector intentionally drops every other variable so that secrets pulled in
  via `envFrom` or `valueFrom.secretKeyRef` do not appear in the file.
- `secrets_mounted.txt` lists secret **names** only. Values are never read.
- `recent_logs_sanitized.txt` is filtered to security-relevant log channels
  (`mapkinase.upload_validation`, `mapkinase.run_guard`,
  `mapkinase.session_workspace`, `RateLimit`, `Cleanup`). Arbitrary stdout
  outside those channels is excluded.
- `tempdirs.txt` lists session directory **names** and modification times only
  — no file contents are read.
- `deployment_validation/results/` is already in `.gitignore`. Do not check in
  the generated `openshift_evidence/` directory unless you have manually
  reviewed each file for sensitive content.

## Re-running

All scripts are idempotent. Re-running the collector overwrites
`openshift_evidence/`. Re-running either summarizer overwrites its corresponding
`*_report.md`. Re-running the URL probe overwrites both `url_check_results.json`
and `url_check_summary.md`.

## Exit codes

The Python scripts (`report_runtime_settings.py`,
`summarize_openshift_controls.py`, `check_direct_url_protections.py`) follow the
same convention:

- `0` — all checks PASS.
- `1` — at least one FAIL.
- `2` — input files missing or bad arguments.
- `3` — no FAILs but at least one REVIEW.

This makes it easy to wire any of them into CI later if desired.
