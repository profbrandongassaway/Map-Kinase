# Map-Kinase Dependency Security Review

Review date: 2026-06-11

Production entry point: `app.py`

Python used for review: Python 3.12.10

Production install command:

```powershell
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install pip-audit
python -m pip_audit -r requirements.txt
```

Development/maintenance audit command:

```powershell
python -m pip_audit -r requirements-dev.txt
```

## Baseline Production Audit

The baseline audit was run from a clean Python 3.12 virtual environment created for this review, before dependency remediation. Exact `pip-audit` output:

```text
Name             Version ID             Fix Versions
---------------- ------- -------------- ------------
requests         2.32.5  CVE-2026-25645 2.33.0
pillow           12.0.0  PYSEC-2026-165 12.2.0
pillow           12.0.0  PYSEC-2026-165 12.2.0
pillow           12.0.0  CVE-2026-25990 12.1.1
pillow           12.0.0  CVE-2026-40192 12.2.0
pillow           12.0.0  CVE-2026-42309 12.2.0
pillow           12.0.0  CVE-2026-42310 12.2.0
pillow           12.0.0  CVE-2026-42311 12.2.0
python-multipart 0.0.21  CVE-2026-24486 0.0.22
python-multipart 0.0.21  CVE-2026-40347 0.0.26
python-multipart 0.0.21  CVE-2026-42561 0.0.27
Found 11 known vulnerabilities in 3 packages
```

## Remediation Table

| Package | Previous version | Updated version | Direct or transitive dependency | Used in production runtime? | Directly imported by Map-Kinase? | pip-audit status after update | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| requests | 2.32.5 | 2.33.0 | Direct runtime pin | Yes | Yes | Clean | Used for KEGG, WikiPathways, UniProt, CST, and related HTTP fetch workflows. Updated to the fixed version reported by `pip-audit`. |
| Pillow | 12.0.0 | 12.2.0 | Direct runtime pin | Yes | Yes, as `PIL` | Clean | Used for pathway/image handling. Updated to cover all Pillow advisories reported by `pip-audit`. |
| python-multipart | 0.0.21 | 0.0.27 | Direct runtime pin; also required by Shiny | Yes | No | Clean | Used by the Shiny/web upload stack. Updated to the highest fixed version required by the reported advisories. |
| lxml | Runtime resolved 6.1.1; dev pin was 6.0.2 | 6.1.1 | Runtime transitive via `pywikipathways`; dev direct optional utility | Yes, transitive | No | Clean | Production already resolved to 6.1.1. Dev pin was separately updated from 6.0.2 to 6.1.1 after `requirements-dev.txt` audit reported PYSEC-2026-87. |
| idna | 3.18 | 3.18 | Transitive via `requests` and `anyio` | Yes, transitive | No | Clean | No vulnerability reported in the clean production audit. |
| urllib3 | 2.7.0 | 2.7.0 | Transitive via `requests` | Yes, transitive | No | Clean | No vulnerability reported in the clean production audit. |
| starlette | 1.3.0 | 1.3.0 | Transitive via `shiny` | Yes, transitive | No | Clean | No vulnerability reported in the clean production audit. |
| orjson | 3.11.9 | 3.11.9 | Transitive via `shiny` | Yes, transitive | No | Clean | No vulnerability reported in the clean production audit. |
| bokeh | Not installed from production requirements | Not installed | Not in production dependency graph | No | No | Not present | Kept out of production runtime. |
| tornado | Not installed from production requirements | Not installed | Not in production dependency graph | No | No | Not present | Kept out of production runtime. |

## Post-Remediation Audit Status

Production `requirements.txt` now passes `pip-audit`:

```text
No known vulnerabilities found
```

Development/maintenance `requirements-dev.txt` also passes `pip-audit` after updating the dev-only `lxml` pin:

```text
No known vulnerabilities found
```

## Verification

The following checks passed after reinstalling `requirements.txt` in a clean Python 3.12 virtual environment:

```powershell
python MapKinase_WebApp/test_dependency_imports_smoke.py
python MapKinase_WebApp/test_upload_validation_smoke.py
python MapKinase_WebApp/test_limits_and_throttle_smoke.py
python MapKinase_WebApp/test_session_workspace_cleanup_smoke.py
python MapKinase_WebApp/test_no_code_execution_smoke.py
python MapKinase_WebApp/test_deployment_debug_surface_smoke.py
python -m compileall MapKinase_WebApp
```

## Exceptions

No production runtime vulnerabilities remain from `pip-audit -r requirements.txt` as of 2026-06-11.

No development/maintenance vulnerabilities remain from `pip-audit -r requirements-dev.txt` as of 2026-06-11.

Existing application-level compensating controls verified by smoke tests include strict upload type validation, upload size limits, per-session throttling, session-scoped temporary storage, no direct URL access to uploaded artifacts, no user-controlled code execution, and production-mode gating of debug surfaces.
