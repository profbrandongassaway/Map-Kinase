"""Summarize OpenShift admin controls from collected pod/deployment YAML.

Reads:
    deployment_validation/results/openshift_evidence/deployment.yaml
    deployment_validation/results/openshift_evidence/pod.yaml
    deployment_validation/results/openshift_evidence/scc.txt
    deployment_validation/results/openshift_evidence/networkpolicy.yaml
    deployment_validation/results/openshift_evidence/secrets_mounted.txt

Writes:
    deployment_validation/results/openshift_controls_report.md

Covers Section 8 of DEPLOYMENT_VALIDATION_PLAN.md:
    - CPU/memory requests + limits (per container)
    - Least-privilege container settings (runAsNonRoot, runAsUser,
      allowPrivilegeEscalation, readOnlyRootFilesystem, capabilities.drop)
    - SCC annotation
    - ServiceAccount name
    - Mounted secret names (names only)
    - NetworkPolicy names

PyYAML is required. If missing, the script prints an actionable install hint.

This script does not contact the cluster; it only parses files written by
collect_openshift_evidence.sh.
"""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

try:
    import yaml  # type: ignore
except ImportError:
    print(
        "ERROR: PyYAML is required for this script.\n"
        "Install with:  pip install pyyaml",
        file=sys.stderr,
    )
    raise SystemExit(2)


REPO_ROOT = Path(__file__).resolve().parents[2]
EVIDENCE_DIR = REPO_ROOT / "deployment_validation" / "results" / "openshift_evidence"
DEPLOY_FILE = EVIDENCE_DIR / "deployment.yaml"
POD_FILE = EVIDENCE_DIR / "pod.yaml"
SCC_FILE = EVIDENCE_DIR / "scc.txt"
NETPOL_FILE = EVIDENCE_DIR / "networkpolicy.yaml"
SECRETS_FILE = EVIDENCE_DIR / "secrets_mounted.txt"
OUT_FILE = REPO_ROOT / "deployment_validation" / "results" / "openshift_controls_report.md"


@dataclass
class Row:
    control: str
    expected: str
    actual: str
    status: str  # PASS / REVIEW / FAIL
    notes: str


def _load_yaml(path: Path) -> dict | None:
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8", errors="replace")
    # Filter out leading comment-only documents.
    if not text.strip() or text.strip().startswith("#"):
        # `oc get networkpolicy -o yaml` returns an empty `items:` list when
        # none exist; that still parses fine. Pure comment files (placeholder)
        # would parse to None.
        try:
            data = yaml.safe_load(text)
        except yaml.YAMLError:
            return None
        return data if isinstance(data, dict) else None
    try:
        return yaml.safe_load(text)
    except yaml.YAMLError as exc:
        print(f"WARNING: failed to parse {path}: {exc}", file=sys.stderr)
        return None


def _get_pod_spec(pod_doc: dict) -> dict:
    return (pod_doc or {}).get("spec", {}) or {}


def _get_containers(spec: dict) -> list[dict]:
    return list(spec.get("containers") or [])


def _esc(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def check_resources(containers: list[dict]) -> list[Row]:
    rows: list[Row] = []
    if not containers:
        rows.append(
            Row(
                control="CPU/memory limits",
                expected="requests + limits on every container",
                actual="<no containers found>",
                status="FAIL",
                notes="Pod spec contained no containers. Re-run evidence collection.",
            )
        )
        return rows

    for c in containers:
        name = c.get("name", "<unnamed>")
        resources = c.get("resources") or {}
        requests = resources.get("requests") or {}
        limits = resources.get("limits") or {}

        cpu_req = requests.get("cpu", "<unset>")
        mem_req = requests.get("memory", "<unset>")
        cpu_lim = limits.get("cpu", "<unset>")
        mem_lim = limits.get("memory", "<unset>")

        missing = [
            label
            for label, val in (
                ("cpu request", cpu_req),
                ("memory request", mem_req),
                ("cpu limit", cpu_lim),
                ("memory limit", mem_lim),
            )
            if val == "<unset>"
        ]

        status = "PASS" if not missing else ("FAIL" if "cpu limit" in missing or "memory limit" in missing else "REVIEW")
        notes = (
            "All requests and limits set."
            if not missing
            else f"Missing: {', '.join(missing)}. OIT expects both CPU and memory limits per container."
        )
        rows.append(
            Row(
                control=f"Resources / container `{name}`",
                expected="requests + limits for cpu and memory",
                actual=f"req: cpu={cpu_req}, mem={mem_req}; lim: cpu={cpu_lim}, mem={mem_lim}",
                status=status,
                notes=notes,
            )
        )
    return rows


def check_security_context(spec: dict, containers: list[dict]) -> list[Row]:
    rows: list[Row] = []
    pod_sc = spec.get("securityContext") or {}

    pod_run_as_non_root = pod_sc.get("runAsNonRoot")
    pod_run_as_user = pod_sc.get("runAsUser")

    rows.append(
        Row(
            control="pod.securityContext.runAsNonRoot",
            expected="true",
            actual=str(pod_run_as_non_root) if pod_run_as_non_root is not None else "<unset>",
            status="PASS" if pod_run_as_non_root is True else ("REVIEW" if pod_run_as_non_root is None else "FAIL"),
            notes=(
                "Pod-level non-root enforced."
                if pod_run_as_non_root is True
                else "Check container-level runAsNonRoot below; OpenShift SCC may also enforce this."
            ),
        )
    )
    rows.append(
        Row(
            control="pod.securityContext.runAsUser",
            expected="non-root UID (typically OpenShift-assigned in restricted SCC)",
            actual=str(pod_run_as_user) if pod_run_as_user is not None else "<unset, SCC-assigned>",
            status="PASS" if pod_run_as_user is None or (isinstance(pod_run_as_user, int) and pod_run_as_user != 0) else "FAIL",
            notes=(
                "OpenShift restricted SCC will assign a non-root UID at runtime if unset."
                if pod_run_as_user is None
                else ("Non-root UID." if pod_run_as_user != 0 else "Pod is configured to run as UID 0.")
            ),
        )
    )

    for c in containers:
        name = c.get("name", "<unnamed>")
        sc = c.get("securityContext") or {}
        run_as_non_root = sc.get("runAsNonRoot")
        if run_as_non_root is None:
            run_as_non_root = pod_run_as_non_root
        allow_priv_esc = sc.get("allowPrivilegeEscalation")
        read_only_fs = sc.get("readOnlyRootFilesystem")
        privileged = sc.get("privileged")
        caps = (sc.get("capabilities") or {})
        caps_drop = caps.get("drop") or []
        caps_add = caps.get("add") or []

        # runAsNonRoot
        rows.append(
            Row(
                control=f"container `{name}`.runAsNonRoot",
                expected="true (or enforced by SCC)",
                actual=str(run_as_non_root) if run_as_non_root is not None else "<unset>",
                status="PASS" if run_as_non_root is True else ("REVIEW" if run_as_non_root is None else "FAIL"),
                notes="Restricted SCC enforces non-root at admission time even if unset.",
            )
        )
        # allowPrivilegeEscalation
        rows.append(
            Row(
                control=f"container `{name}`.allowPrivilegeEscalation",
                expected="false",
                actual=str(allow_priv_esc) if allow_priv_esc is not None else "<unset>",
                status="PASS" if allow_priv_esc is False else ("REVIEW" if allow_priv_esc is None else "FAIL"),
                notes="Should be explicitly false for restricted-v2 SCC.",
            )
        )
        # readOnlyRootFilesystem
        rows.append(
            Row(
                control=f"container `{name}`.readOnlyRootFilesystem",
                expected="true (preferred); writable mounts only for temp/session dirs",
                actual=str(read_only_fs) if read_only_fs is not None else "<unset>",
                status="PASS" if read_only_fs is True else "REVIEW",
                notes=(
                    "Read-only root with explicit writable volumes for temp/session is the safer config; "
                    "if unset/false, document why."
                ),
            )
        )
        # privileged
        rows.append(
            Row(
                control=f"container `{name}`.privileged",
                expected="false",
                actual=str(privileged) if privileged is not None else "<unset>",
                status="PASS" if not privileged else "FAIL",
                notes="Privileged containers are not appropriate for this workload.",
            )
        )
        # capabilities
        drop_str = ",".join(caps_drop) if caps_drop else "<none listed>"
        add_str = ",".join(caps_add) if caps_add else "<none>"
        drops_all = any(str(d).upper() == "ALL" for d in caps_drop)
        rows.append(
            Row(
                control=f"container `{name}`.capabilities",
                expected='drop: ["ALL"], add: <none>',
                actual=f"drop=[{drop_str}], add=[{add_str}]",
                status="PASS" if drops_all and not caps_add else ("FAIL" if caps_add else "REVIEW"),
                notes=(
                    "Drops ALL capabilities, none added."
                    if drops_all and not caps_add
                    else (
                        f"Capabilities added: {add_str}. Verify each is required."
                        if caps_add
                        else "Capabilities should explicitly drop ALL."
                    )
                ),
            )
        )

    return rows


def check_scc() -> list[Row]:
    if not SCC_FILE.exists():
        return [
            Row(
                control="SCC annotation",
                expected="restricted-v2 (or equivalent restricted SCC)",
                actual="<scc.txt missing>",
                status="REVIEW",
                notes="Re-run evidence collection.",
            )
        ]
    text = SCC_FILE.read_text(encoding="utf-8", errors="replace")
    # The annotation value is whatever comes after the first heading and before
    # the next blank line.
    match = re.search(r"openshift\.io/scc annotation\):\s*\n([^\n]+)", text)
    scc_value = (match.group(1).strip() if match else "<not found>").strip()

    if scc_value in {"restricted", "restricted-v2"}:
        status = "PASS"
        notes = f"Standard restricted SCC ({scc_value})."
    elif scc_value in {"<not found>", ""}:
        status = "REVIEW"
        notes = "Could not parse SCC annotation. Inspect scc.txt manually."
    elif scc_value in {"privileged", "anyuid", "hostaccess", "hostmount-anyuid"}:
        status = "FAIL"
        notes = f"Elevated SCC '{scc_value}' assigned. Confirm this is intentional and documented; restricted-v2 is preferred."
    else:
        status = "REVIEW"
        notes = f"Non-default SCC '{scc_value}'. Confirm it grants no more than required."

    return [
        Row(
            control="SCC annotation",
            expected="restricted-v2 (or equivalent restricted SCC)",
            actual=scc_value or "<empty>",
            status=status,
            notes=notes,
        )
    ]


def check_service_account(spec: dict) -> list[Row]:
    sa = spec.get("serviceAccountName") or spec.get("serviceAccount") or "default"
    if sa == "default":
        status = "REVIEW"
        notes = (
            "Pod uses the namespace 'default' service account. A dedicated SA with no bound roles "
            "is preferred; confirm rolebindings.txt shows no extra permissions."
        )
    else:
        status = "PASS"
        notes = "Pod uses a dedicated service account. Verify rolebindings.txt shows minimal permissions."
    return [
        Row(
            control="ServiceAccount",
            expected="dedicated SA with minimal permissions",
            actual=str(sa),
            status=status,
            notes=notes,
        )
    ]


def check_secrets() -> list[Row]:
    if not SECRETS_FILE.exists():
        return [
            Row(
                control="Mounted secrets",
                expected="names only; minimum required",
                actual="<secrets_mounted.txt missing>",
                status="REVIEW",
                notes="Re-run evidence collection.",
            )
        ]
    text = SECRETS_FILE.read_text(encoding="utf-8", errors="replace")

    def _names_under(heading: str) -> list[str]:
        # Section starts at "## <heading>:" and ends at next "## " or EOF.
        pattern = rf"## {re.escape(heading)}:\s*\n(.*?)(?=\n## |\Z)"
        m = re.search(pattern, text, flags=re.DOTALL)
        if not m:
            return []
        body = m.group(1)
        names: list[str] = []
        for line in body.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # Lines like "vol-name -> secret/secret-name" or just "secret-name"
            if "->" in line:
                _, _, right = line.partition("->")
                right = right.strip()
                if right.startswith("secret/"):
                    right = right[len("secret/") :]
                names.append(right)
            else:
                names.append(line)
        return [n for n in names if n]

    mounted = _names_under("Mounted as volumes")
    envfrom = _names_under("Exposed via envFrom secretRef")
    valueref = _names_under("Referenced via env valueFrom.secretKeyRef")
    pull = _names_under("ImagePullSecrets")

    rows: list[Row] = []
    rows.append(
        Row(
            control="Secrets mounted as volumes",
            expected="only required secrets, no extraneous mounts",
            actual=", ".join(mounted) if mounted else "<none>",
            status="REVIEW" if mounted else "PASS",
            notes="Confirm each mounted secret is required by the application.",
        )
    )
    rows.append(
        Row(
            control="Secrets via envFrom",
            expected="only required secrets",
            actual=", ".join(envfrom) if envfrom else "<none>",
            status="REVIEW" if envfrom else "PASS",
            notes="Confirm each secret's keys are required by the application.",
        )
    )
    rows.append(
        Row(
            control="Secrets via env valueFrom.secretKeyRef",
            expected="only required keys",
            actual=", ".join(valueref) if valueref else "<none>",
            status="REVIEW" if valueref else "PASS",
            notes="Each entry exposes a single key; verify it is needed.",
        )
    )
    rows.append(
        Row(
            control="ImagePullSecrets",
            expected="zero or one (registry credentials)",
            actual=", ".join(pull) if pull else "<none>",
            status="PASS",
            notes="Typical for a private registry pull credential.",
        )
    )
    return rows


def check_networkpolicy() -> list[Row]:
    if not NETPOL_FILE.exists():
        return [
            Row(
                control="NetworkPolicy",
                expected="ingress/egress restricted per OIT policy",
                actual="<networkpolicy.yaml missing>",
                status="REVIEW",
                notes="Re-run evidence collection.",
            )
        ]
    doc = _load_yaml(NETPOL_FILE)
    items: list[dict] = []
    if isinstance(doc, dict):
        items = list(doc.get("items") or [])
    names = [str((item.get("metadata") or {}).get("name", "<unnamed>")) for item in items]
    if not items:
        return [
            Row(
                control="NetworkPolicy",
                expected="ingress/egress restricted per OIT policy",
                actual="<no NetworkPolicy in namespace>",
                status="REVIEW",
                notes="No NetworkPolicy objects in this namespace. If egress is unrestricted by cluster default, document the reason.",
            )
        ]
    return [
        Row(
            control="NetworkPolicy objects",
            expected="ingress/egress restricted per OIT policy",
            actual=", ".join(names),
            status="REVIEW",
            notes="Confirm each policy's selectors and rules match the documented egress/ingress requirements.",
        )
    ]


def check_image_provenance(spec: dict, pod_doc: dict) -> list[Row]:
    statuses = ((pod_doc or {}).get("status") or {}).get("containerStatuses") or []
    image_id = ""
    image_ref = ""
    if statuses:
        image_id = str(statuses[0].get("imageID", ""))
        image_ref = str(statuses[0].get("image", ""))
    if "@sha256:" in image_id or "@sha256:" in image_ref:
        status = "PASS"
        notes = "Image is pinned by digest."
    elif image_ref and ":" in image_ref and not image_ref.endswith(":latest"):
        status = "REVIEW"
        notes = "Image is referenced by tag, not digest. Digest pinning is preferred for OIT provenance."
    else:
        status = "FAIL"
        notes = "Image reference does not include a tag or digest, or is :latest. Pin to a digest."
    return [
        Row(
            control="Image provenance",
            expected="image pinned by digest (@sha256:...)",
            actual=f"image={image_ref or '<unknown>'}; imageID={image_id or '<unknown>'}",
            status=status,
            notes=notes,
        )
    ]


def render(rows: list[Row]) -> str:
    counts = {"PASS": 0, "REVIEW": 0, "FAIL": 0}
    for r in rows:
        counts[r.status] = counts.get(r.status, 0) + 1

    lines: list[str] = []
    lines.append("# Map-Kinase OpenShift Controls Report")
    lines.append("")
    lines.append("Generated by `deployment_validation/openshift/summarize_openshift_controls.py`.")
    lines.append("")
    lines.append("Sources:")
    for path in (DEPLOY_FILE, POD_FILE, SCC_FILE, NETPOL_FILE, SECRETS_FILE):
        lines.append(f"- `{path.relative_to(REPO_ROOT)}`")
    lines.append("")
    lines.append(
        f"**Totals:** PASS = {counts.get('PASS', 0)}, "
        f"REVIEW = {counts.get('REVIEW', 0)}, "
        f"FAIL = {counts.get('FAIL', 0)}"
    )
    lines.append("")
    lines.append("| Control | Expected | Actual | Status | Notes |")
    lines.append("| --- | --- | --- | --- | --- |")
    for r in rows:
        lines.append(
            f"| {_esc(r.control)} | {_esc(r.expected)} | {_esc(r.actual)} | **{r.status}** | {_esc(r.notes)} |"
        )
    lines.append("")
    lines.append("## Narrative")
    lines.append("")
    if counts.get("FAIL", 0) > 0:
        lines.append(
            "One or more OpenShift controls did not match the expected production configuration. "
            "Address the FAIL rows above before sending the OIT response."
        )
    elif counts.get("REVIEW", 0) > 0:
        lines.append(
            "Core controls (security context, restricted SCC, resource limits) look correct. "
            "The REVIEW rows are items where the OpenShift defaults are typically sufficient but "
            "OIT may want explicit confirmation (e.g. NetworkPolicy selectors, readOnlyRootFilesystem, "
            "dedicated ServiceAccount with no extra rolebindings)."
        )
    else:
        lines.append(
            "All checked OpenShift controls match the expected production configuration: "
            "restricted SCC, non-root container, dropped capabilities, CPU/memory limits, "
            "dedicated service account, and image pinned by digest."
        )
    lines.append("")
    return "\n".join(lines)


def main(argv: list[str]) -> int:
    if not EVIDENCE_DIR.exists():
        print(
            f"ERROR: {EVIDENCE_DIR} not found.\n"
            "Run `bash deployment_validation/openshift/collect_openshift_evidence.sh` first.",
            file=sys.stderr,
        )
        return 2
    if not POD_FILE.exists():
        print(f"ERROR: {POD_FILE} not found. Re-run evidence collection.", file=sys.stderr)
        return 2

    pod_doc = _load_yaml(POD_FILE) or {}
    spec = _get_pod_spec(pod_doc)
    containers = _get_containers(spec)

    rows: list[Row] = []
    rows.extend(check_resources(containers))
    rows.extend(check_security_context(spec, containers))
    rows.extend(check_scc())
    rows.extend(check_service_account(spec))
    rows.extend(check_secrets())
    rows.extend(check_networkpolicy())
    rows.extend(check_image_provenance(spec, pod_doc))

    report = render(rows)

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text(report, encoding="utf-8")

    print(report)
    print(f"\nWrote {OUT_FILE.relative_to(REPO_ROOT)}", file=sys.stderr)

    if any(r.status == "FAIL" for r in rows):
        return 1
    if any(r.status == "REVIEW" for r in rows):
        return 3
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
