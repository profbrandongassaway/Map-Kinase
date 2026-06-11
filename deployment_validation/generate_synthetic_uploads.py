"""Generate safe synthetic files for deployed Map-Kinase validation.

The files are intentionally simple and contain no real user data. They are meant
for manual upload testing in the deployed web UI.
"""

from __future__ import annotations

import json
import zipfile
from pathlib import Path


OUT_DIR = Path(__file__).resolve().parent / "generated_uploads"
MB = 1024 * 1024


def _write_text(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8", newline="\n")


def _write_bytes(path: Path, payload: bytes) -> None:
    path.write_bytes(payload)


def _valid_protein_rows(row_count: int) -> str:
    rows = ["Uniprot Id,Gene Symbol,C: A"]
    for idx in range(row_count):
        rows.append(f"P{10000 + idx},GENE{idx},1.0")
    return "\n".join(rows) + "\n"


def _valid_ptm_rows(row_count: int) -> str:
    rows = ["Uniprot Id\tSite Position\tC: A"]
    for idx in range(row_count):
        rows.append(f"P{10000 + idx}\t{idx + 1}\t0.5")
    return "\n".join(rows) + "\n"


def _valid_metabolite_rows(row_count: int) -> str:
    rows = ["HMDB_ID\tC: A"]
    for idx in range(row_count):
        rows.append(f"HMDB{idx:05d}\t2.0")
    return "\n".join(rows) + "\n"


def _write_padded_valid_csv(path: Path, target_size: int) -> None:
    header = "Uniprot Id,Gene Symbol,C: A,T: Padding\n"
    base_row = "P12345,GENE1,1.0,"
    padding_len = max(1, target_size - len(header.encode("utf-8")) - len(base_row.encode("utf-8")) - 1)
    payload = header + base_row + ("A" * padding_len) + "\n"
    path.write_text(payload, encoding="utf-8", newline="\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    files: list[tuple[str, str]] = []

    fixtures: dict[str, str | bytes] = {
        "valid_protein.csv": _valid_protein_rows(2),
        "valid_ptm.tsv": _valid_ptm_rows(2),
        "valid_metabolite.txt": _valid_metabolite_rows(2),
        "analysis_payload.json": json.dumps({"schema_version": 1, "note": "should be rejected for analysis uploads"}) + "\n",
        "custom_pathway_valid.json": json.dumps(
            {
                "schema_version": 1,
                "pathway_id": "SyntheticCustom",
                "pathway_source": "custom",
                "general_data": {"settings": {"mode": "custom"}},
            },
            indent=2,
        )
        + "\n",
        "payload.exe": b"MZ synthetic executable extension test\n",
        "payload.py": "print('synthetic script extension test')\n",
        "payload.js": "console.log('synthetic script extension test')\n",
        "payload.html": "<html><body>synthetic html extension test</body></html>\n",
        "payload.sh": "#!/bin/sh\necho synthetic script extension test\n",
        "payload.bat": "@echo off\necho synthetic script extension test\r\n",
        "data.csv.exe": _valid_protein_rows(1),
        "binary_renamed.csv": bytes([120, 3, 255, 0, 76, 12, 150, 44, 33]),
        "malformed_non_tabular.txt": "This is plain text, not a tab-delimited table.\nAnother sentence without tabs.\n",
    }

    for name, payload in fixtures.items():
        path = OUT_DIR / name
        if isinstance(payload, bytes):
            _write_bytes(path, payload)
        else:
            _write_text(path, payload)
        files.append((name, "generated"))

    zip_path = OUT_DIR / "payload.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("payload.txt", "synthetic zip extension test\n")
    files.append(("payload.zip", "generated"))

    xlsx_path = OUT_DIR / "payload.xlsx"
    with zipfile.ZipFile(xlsx_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        zf.writestr("[Content_Types].xml", "<Types></Types>")
        zf.writestr("xl/workbook.xml", "<workbook></workbook>")
    files.append(("payload.xlsx", "generated"))

    _write_padded_valid_csv(OUT_DIR / "under_10mb_valid_protein.csv", (10 * MB) - 2048)
    files.append(("under_10mb_valid_protein.csv", "valid CSV under default 10 MB per-file limit"))

    _write_padded_valid_csv(OUT_DIR / "over_10mb_valid_shape.csv", (10 * MB) + 1)
    files.append(("over_10mb_valid_shape.csv", "valid-looking CSV over default 10 MB per-file limit"))

    for idx in range(1, 5):
        _write_padded_valid_csv(OUT_DIR / f"combined_under_10mb_{idx}.csv", (8 * MB) + 128)
        files.append((f"combined_under_10mb_{idx}.csv", "four files together exceed default 30 MB total limit"))

    manifest = {
        "output_dir": str(OUT_DIR),
        "notes": [
            "Use these files only for validation.",
            "Browser upload tests cannot usually preserve path traversal filenames.",
            "Path traversal filename rejection is covered by local smoke tests unless OIT uses a controlled multipart tool.",
        ],
        "files": [{"name": name, "purpose": purpose} for name, purpose in files],
    }
    (OUT_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Generated {len(files)} files in {OUT_DIR}")
    print(f"Manifest: {OUT_DIR / 'manifest.json'}")


if __name__ == "__main__":
    main()
