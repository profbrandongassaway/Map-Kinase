#!/usr/bin/env bash
# Collect Map-Kinase OpenShift deployment evidence for the OIT validation plan.
#
# Prerequisites:
#   - `oc` is installed and on PATH.
#   - `oc login` has been run and `oc project <namespace>` points at the
#     Map-Kinase deployment.
#
# Optional environment variables:
#   MAPKINASE_POD_LABEL    Label selector used to find the pod. If unset, the
#                          script tries a list of common Map-Kinase labels
#                          (app=mapkinase, app=map-kinase, deployment=mapkinase,
#                          deploymentconfig=mapkinase, name=mapkinase, and the
#                          hyphenated variants) and uses the first match.
#                          Override examples:
#                            MAPKINASE_POD_LABEL=app=map-kinase
#                            MAPKINASE_POD_LABEL=deployment=mapkinase
#   MAPKINASE_POD_NAME     Exact pod name. If set, MAPKINASE_POD_LABEL is ignored.
#   MAPKINASE_DEPLOYMENT   Deployment name. If unset, the script tries to infer
#                          it from the pod's ownerReferences.
#   MAPKINASE_LOG_WINDOW   Window for `oc logs --since=`. Default: 30m
#
# Outputs (under deployment_validation/results/openshift_evidence/):
#   summary.txt
#   deployment.yaml
#   pod.yaml
#   pod_describe.txt
#   env_mapkinase.txt          (only MAPKINASE_* and M5_* vars)
#   tempdirs.txt
#   serviceaccount.txt
#   scc.txt
#   rolebindings.txt
#   networkpolicy.yaml
#   secrets_mounted.txt        (names only)
#   recent_logs_sanitized.txt  (filtered to security-relevant log channels)
#
# The script never echoes secret values, never reads uploaded file contents,
# and only dumps env vars that match the MAPKINASE_/M5_ allowlist.

set -euo pipefail

# When MAPKINASE_POD_LABEL is unset we try this list in order and use the first
# label that matches a pod. This avoids having to know the deployment's label
# scheme in advance. Override with MAPKINASE_POD_LABEL=<selector> if your pod
# uses something else.
POD_LABEL_DEFAULTS=(
    "app=mapkinase"
    "app=map-kinase"
    "app=mapkinase-app"
    "deployment=mapkinase"
    "deployment=map-kinase"
    "deploymentconfig=mapkinase"
    "deploymentconfig=map-kinase"
    "name=mapkinase"
    "name=map-kinase"
)
POD_LABEL="${MAPKINASE_POD_LABEL:-}"
POD_NAME_OVERRIDE="${MAPKINASE_POD_NAME:-}"
DEPLOY_OVERRIDE="${MAPKINASE_DEPLOYMENT:-}"
LOG_WINDOW="${MAPKINASE_LOG_WINDOW:-30m}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUT_DIR="${REPO_ROOT}/deployment_validation/results/openshift_evidence"
mkdir -p "${OUT_DIR}"

log() { printf '[collect] %s\n' "$*" >&2; }
die() { printf '[collect] ERROR: %s\n' "$*" >&2; exit 1; }

command -v oc >/dev/null 2>&1 || die "'oc' not found on PATH."

# Confirm logged-in context.
if ! oc whoami >/dev/null 2>&1; then
    die "Not logged in to OpenShift. Run 'oc login ...' first."
fi
NAMESPACE="$(oc project -q 2>/dev/null || true)"
if [[ -z "${NAMESPACE}" ]]; then
    die "Could not determine current project. Run 'oc project <namespace>'."
fi
OC_USER="$(oc whoami 2>/dev/null || echo unknown)"
log "Namespace: ${NAMESPACE}"
log "User:      ${OC_USER}"

# Try a single label selector. Echoes the matched pod name on stdout, or empty.
_pod_for_label() {
    local sel="$1"
    local name=""
    name="$(oc get pods -l "${sel}" \
            --field-selector=status.phase=Running \
            -o jsonpath='{.items[0].metadata.name}' 2>/dev/null || true)"
    if [[ -z "${name}" ]]; then
        name="$(oc get pods -l "${sel}" \
                -o jsonpath='{.items[0].metadata.name}' 2>/dev/null || true)"
    fi
    printf '%s' "${name}"
}

# Resolve pod.
POD=""
MATCHED_LABEL=""
if [[ -n "${POD_NAME_OVERRIDE}" ]]; then
    POD="${POD_NAME_OVERRIDE}"
    log "Pod (from MAPKINASE_POD_NAME): ${POD}"
elif [[ -n "${POD_LABEL}" ]]; then
    log "Resolving pod with label selector: ${POD_LABEL}"
    POD="$(_pod_for_label "${POD_LABEL}")"
    [[ -n "${POD}" ]] && MATCHED_LABEL="${POD_LABEL}"
else
    log "Resolving pod by trying common label selectors..."
    for sel in "${POD_LABEL_DEFAULTS[@]}"; do
        candidate="$(_pod_for_label "${sel}")"
        if [[ -n "${candidate}" ]]; then
            POD="${candidate}"
            MATCHED_LABEL="${sel}"
            log "  matched ${sel} -> ${POD}"
            break
        fi
    done
fi
if [[ -z "${POD}" ]]; then
    cat >&2 <<EOF
[collect] ERROR: Could not resolve a Map-Kinase pod in namespace '${NAMESPACE}'.

Tried these label selectors:
$(printf '  - %s\n' "${POD_LABEL_DEFAULTS[@]}")

To find the real label, run:
  oc get pods --show-labels

Then re-run with either:
  MAPKINASE_POD_LABEL=<key>=<value> bash deployment_validation/openshift/collect_openshift_evidence.sh
  # or
  MAPKINASE_POD_NAME=<exact-pod-name> bash deployment_validation/openshift/collect_openshift_evidence.sh
EOF
    exit 1
fi
log "Pod:        ${POD}"
[[ -n "${MATCHED_LABEL}" ]] && log "Matched via: ${MATCHED_LABEL}"

# Resolve deployment name.
if [[ -n "${DEPLOY_OVERRIDE}" ]]; then
    DEPLOY="${DEPLOY_OVERRIDE}"
else
    # ownerReferences -> ReplicaSet -> Deployment
    RS_NAME="$(oc get pod "${POD}" \
        -o jsonpath='{.metadata.ownerReferences[?(@.kind=="ReplicaSet")].name}' \
        2>/dev/null || true)"
    if [[ -n "${RS_NAME}" ]]; then
        DEPLOY="$(oc get rs "${RS_NAME}" \
            -o jsonpath='{.metadata.ownerReferences[?(@.kind=="Deployment")].name}' \
            2>/dev/null || true)"
    fi
    DEPLOY="${DEPLOY:-}"
fi
if [[ -z "${DEPLOY}" ]]; then
    log "Could not auto-resolve Deployment. Set MAPKINASE_DEPLOYMENT to capture deployment.yaml."
fi
log "Deployment: ${DEPLOY:-<unknown>}"

# ---- summary.txt: paste-ready identifiers ------------------------------------
IMAGE_REF="$(oc get pod "${POD}" -o jsonpath='{.status.containerStatuses[0].image}' 2>/dev/null || true)"
IMAGE_ID="$(oc get pod "${POD}" -o jsonpath='{.status.containerStatuses[0].imageID}' 2>/dev/null || true)"
PY_VERSION="$(oc exec "${POD}" -- python --version 2>&1 | head -1 || true)"
GIT_COMMIT="$(
    oc exec "${POD}" -- sh -c '
        if [ -f /app/GIT_COMMIT ]; then
            head -1 /app/GIT_COMMIT
        elif [ -f /app/.git/HEAD ]; then
            (cd /app && git rev-parse HEAD 2>/dev/null) || echo "<git rev-parse failed>"
        elif [ -n "${MAPKINASE_GIT_COMMIT:-}" ]; then
            echo "${MAPKINASE_GIT_COMMIT}"
        else
            echo "<unknown: no /app/GIT_COMMIT, no .git, no MAPKINASE_GIT_COMMIT>"
        fi
    ' 2>/dev/null | tr -d '\r' || true
)"
ROUTE_HOST="$(oc get route -o jsonpath='{.items[?(@.spec.to.name=="'"${DEPLOY}"'")].spec.host}' 2>/dev/null || true)"
if [[ -z "${ROUTE_HOST}" ]]; then
    # Fall back to first route in the project.
    ROUTE_HOST="$(oc get route -o jsonpath='{.items[0].spec.host}' 2>/dev/null || true)"
fi
DEPLOYED_URL=""
if [[ -n "${ROUTE_HOST}" ]]; then
    DEPLOYED_URL="https://${ROUTE_HOST}"
fi

SUMMARY="${OUT_DIR}/summary.txt"
{
    echo "Map-Kinase OpenShift evidence summary"
    echo "Collected: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo "Collected by: ${OC_USER}"
    echo
    echo "Namespace:    ${NAMESPACE}"
    echo "Deployment:   ${DEPLOY:-<unknown>}"
    echo "Pod:          ${POD}"
    echo "Deployed URL: ${DEPLOYED_URL:-<no route detected; check with: oc get route>}"
    echo "Image ref:    ${IMAGE_REF:-<unknown>}"
    echo "Image ID:     ${IMAGE_ID:-<unknown>}"
    echo "Git commit:   ${GIT_COMMIT:-<unknown>}"
    echo "Python:       ${PY_VERSION:-<unknown>}"
    echo
    echo "Paste-ready block for DEPLOYMENT_VALIDATION_PLAN.md Pre-Validation Setup:"
    echo "  - Deployed URL:               ${DEPLOYED_URL:-https://DOMAIN_HERE}"
    echo "  - Python version:             ${PY_VERSION:-<fill in>}"
    echo "  - Image/tag:                  ${IMAGE_REF:-<fill in>}"
    echo "  - Image digest:               ${IMAGE_ID:-<fill in>}"
    echo "  - Git commit:                 ${GIT_COMMIT:-<fill in>}"
    echo "  - OpenShift project:          ${NAMESPACE}"
    echo "  - Pod name:                   ${POD}"
    echo "  - Validation date (UTC):      $(date -u +'%Y-%m-%d')"
} > "${SUMMARY}"
log "Wrote ${SUMMARY}"

# ---- YAML / describe --------------------------------------------------------
if [[ -n "${DEPLOY}" ]]; then
    oc get deploy "${DEPLOY}" -o yaml > "${OUT_DIR}/deployment.yaml" 2>/dev/null \
        || echo "# Failed to get deployment ${DEPLOY}" > "${OUT_DIR}/deployment.yaml"
else
    echo "# Deployment name not resolved. Set MAPKINASE_DEPLOYMENT and re-run." \
        > "${OUT_DIR}/deployment.yaml"
fi
log "Wrote deployment.yaml"

oc get pod "${POD}" -o yaml > "${OUT_DIR}/pod.yaml"
log "Wrote pod.yaml"

oc describe pod "${POD}" > "${OUT_DIR}/pod_describe.txt" 2>&1 || true
log "Wrote pod_describe.txt"

# ---- env (allowlist only) ---------------------------------------------------
# Only MAPKINASE_* and M5_* are exported. Never dump full env to avoid secrets.
ENV_FILE="${OUT_DIR}/env_mapkinase.txt"
{
    echo "# Filtered to MAPKINASE_* and M5_* environment variables only."
    echo "# Other env vars (incl. secrets) intentionally excluded."
    echo "# Source: oc exec ${POD} -- printenv"
    echo "# Collected: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo
    oc exec "${POD}" -- printenv 2>/dev/null \
        | grep -E '^(MAPKINASE_|M5_)' \
        | LC_ALL=C sort \
        || echo "# (no MAPKINASE_* / M5_* env vars found)"
} > "${ENV_FILE}"
log "Wrote env_mapkinase.txt"

# ---- temp/session workspace listing -----------------------------------------
TEMP_FILE="${OUT_DIR}/tempdirs.txt"
{
    echo "# Temp/session workspace listing collected via oc exec python."
    echo "# Names + mtimes only. No file contents are read."
    echo "# Collected: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo
    oc exec "${POD}" -- python - <<'PY' 2>&1 || echo "# (python probe failed)"
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

tmp_root = Path(tempfile.gettempdir()).resolve()
print(f"tempfile.gettempdir() = {tmp_root}")

env_root = os.environ.get("MAPKINASE_SESSION_TMP_ROOT")
print(f"MAPKINASE_SESSION_TMP_ROOT = {env_root or '<unset>'}")

if env_root:
    session_root = Path(env_root).resolve()
else:
    session_root = (tmp_root / "mapkinase_sessions").resolve()
print(f"resolved session root = {session_root}")
print(f"exists = {session_root.exists()}")
print()

if session_root.exists() and session_root.is_dir():
    entries = sorted(session_root.iterdir(), key=lambda p: p.name)[:50]
    print(f"first {len(entries)} entries (name, mtime UTC):")
    for entry in entries:
        try:
            mtime = datetime.fromtimestamp(entry.stat().st_mtime, tz=timezone.utc)
            print(f"  {entry.name}\t{mtime.isoformat()}")
        except OSError as exc:
            print(f"  {entry.name}\t<stat failed: {exc}>")
else:
    print("(session root does not exist yet; created lazily on first session)")
PY
} > "${TEMP_FILE}"
log "Wrote tempdirs.txt"

# ---- ServiceAccount / SCC / RoleBindings ------------------------------------
SA_NAME="$(oc get pod "${POD}" -o jsonpath='{.spec.serviceAccountName}' 2>/dev/null || true)"
SA_NAME="${SA_NAME:-default}"
{
    echo "# ServiceAccount used by pod: ${SA_NAME}"
    echo
    oc get sa "${SA_NAME}" -o yaml 2>&1 || true
} > "${OUT_DIR}/serviceaccount.txt"
log "Wrote serviceaccount.txt (sa=${SA_NAME})"

{
    echo "# SCC assigned to pod (openshift.io/scc annotation):"
    oc get pod "${POD}" -o jsonpath='{.metadata.annotations.openshift\.io/scc}' 2>/dev/null
    echo
    echo
    echo "# Pod securityContext:"
    oc get pod "${POD}" -o jsonpath='{.spec.securityContext}' 2>/dev/null
    echo
    echo
    echo "# Container securityContexts:"
    oc get pod "${POD}" -o jsonpath='{range .spec.containers[*]}{.name}: {.securityContext}{"\n"}{end}' 2>/dev/null
} > "${OUT_DIR}/scc.txt"
log "Wrote scc.txt"

{
    echo "# RoleBindings in namespace ${NAMESPACE}:"
    oc get rolebindings -o wide 2>&1 || true
    echo
    echo "# RoleBindings referencing serviceaccount ${SA_NAME}:"
    oc get rolebindings -o wide 2>/dev/null \
        | grep -E "(NAME|${SA_NAME})" || echo "(none found)"
} > "${OUT_DIR}/rolebindings.txt"
log "Wrote rolebindings.txt"

# ---- NetworkPolicy ----------------------------------------------------------
{
    echo "# NetworkPolicy objects in namespace ${NAMESPACE}"
    echo
    oc get networkpolicy -o yaml 2>&1 || true
} > "${OUT_DIR}/networkpolicy.yaml"
log "Wrote networkpolicy.yaml"

# ---- Mounted secret names (names only) --------------------------------------
{
    echo "# Secrets referenced by pod (names only, no values)."
    echo
    echo "## Mounted as volumes:"
    oc get pod "${POD}" -o jsonpath='{range .spec.volumes[?(@.secret)]}{.name} -> secret/{.secret.secretName}{"\n"}{end}' 2>/dev/null || true
    echo
    echo "## Exposed via envFrom secretRef:"
    oc get pod "${POD}" -o jsonpath='{range .spec.containers[*]}{.name}: {range .envFrom[?(@.secretRef)]}{.secretRef.name} {end}{"\n"}{end}' 2>/dev/null || true
    echo
    echo "## Referenced via env valueFrom.secretKeyRef:"
    oc get pod "${POD}" -o jsonpath='{range .spec.containers[*]}{.name}: {range .env[?(@.valueFrom.secretKeyRef)]}{.valueFrom.secretKeyRef.name}/{.valueFrom.secretKeyRef.key} {end}{"\n"}{end}' 2>/dev/null || true
    echo
    echo "## ImagePullSecrets:"
    oc get pod "${POD}" -o jsonpath='{range .spec.imagePullSecrets[*]}{.name}{"\n"}{end}' 2>/dev/null || true
} > "${OUT_DIR}/secrets_mounted.txt"
log "Wrote secrets_mounted.txt"

# ---- Sanitized recent logs --------------------------------------------------
# Capture only lines from the security-relevant log channels. Anything outside
# this allowlist (e.g. arbitrary stdout that might include user data) is
# excluded. If no matches, write a clear note.
LOG_FILE="${OUT_DIR}/recent_logs_sanitized.txt"
{
    echo "# Log lines from last ${LOG_WINDOW}, filtered to security-relevant channels:"
    echo "#   mapkinase.upload_validation"
    echo "#   mapkinase.run_guard"
    echo "#   mapkinase.session_workspace"
    echo "#   RateLimit / Cleanup markers"
    echo "# Source: oc logs ${POD} --since=${LOG_WINDOW}"
    echo "# Collected: $(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo
    oc logs "${POD}" --since="${LOG_WINDOW}" 2>/dev/null \
        | grep -E 'mapkinase\.(upload_validation|run_guard|session_workspace)|RateLimit|Cleanup' \
        || echo "(no matching log lines in window; run a few uploads / Run clicks and re-collect)"
} > "${LOG_FILE}"
log "Wrote recent_logs_sanitized.txt"

# ---- Final summary ----------------------------------------------------------
log ""
log "All evidence written under: ${OUT_DIR}"
log ""
log "Paste-ready identifier block:"
sed -n '/Paste-ready block/,$p' "${SUMMARY}" | sed 's/^/    /' >&2

cat <<EOF

Next steps:
  1. Review files in ${OUT_DIR}
  2. python deployment_validation/openshift/report_runtime_settings.py
  3. python deployment_validation/openshift/summarize_openshift_controls.py
  4. python deployment_validation/check_direct_url_protections.py ${DEPLOYED_URL:-https://<deployed-url>}
EOF
