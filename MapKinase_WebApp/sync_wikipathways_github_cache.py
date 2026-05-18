#!/usr/bin/env python3
"""
Download WikiPathways pathway content from the upstream GitHub repositories and
sync it into the local species-organized cache used by Map-Kinase.

Sources:
- wikipathways/wikipathways-database: GPML + pathway metadata files
- wikipathways/wikipathways-assets: generated JSON/PNG/SVG assets
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import shutil
import subprocess
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple


LOGGER = logging.getLogger("sync_wikipathways_github_cache")
CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parent
DEFAULT_REPO_CACHE_DIR = PROJECT_ROOT / "cache" / "wikipathways_github"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "stored_pathways" / "wikipathways"
DEFAULT_MANIFEST_PATH = DEFAULT_OUTPUT_DIR / "_github_sync_manifest.json"
SPECIES_REF_FILE = CURRENT_DIR / "species_ref_list.csv"
DATABASE_REPO_URL = "https://github.com/wikipathways/wikipathways-database.git"
ASSETS_REPO_URL = "https://github.com/wikipathways/wikipathways-assets.git"


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sync WikiPathways GitHub pathway files into stored_pathways/wikipathways."
    )
    parser.add_argument(
        "--repo-cache-dir",
        default=str(DEFAULT_REPO_CACHE_DIR),
        help=f"Directory used to clone the upstream repositories (default: {DEFAULT_REPO_CACHE_DIR}).",
    )
    parser.add_argument(
        "--out-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help=f"Destination directory for species-organized pathway files (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--manifest",
        default=str(DEFAULT_MANIFEST_PATH),
        help=f"Output manifest path (default: {DEFAULT_MANIFEST_PATH}).",
    )
    parser.add_argument(
        "--org",
        action="append",
        default=[],
        help="Restrict sync to one or more organism codes from species_ref_list.csv, such as hsa or mmu.",
    )
    parser.add_argument(
        "--max-pathways",
        type=int,
        default=None,
        help="For testing, sync at most N pathways after filtering.",
    )
    parser.add_argument(
        "--skip-fetch",
        action="store_true",
        help="Use the already-cloned repositories in --repo-cache-dir without pulling updates.",
    )
    parser.add_argument(
        "--pretty",
        action="store_true",
        help="Pretty-print the sync manifest JSON.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return parser.parse_args(list(argv) if argv is not None else None)


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_json_atomic(path: Path, obj: object, pretty: bool = False) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as fh:
        if pretty:
            json.dump(obj, fh, indent=2, ensure_ascii=False)
        else:
            json.dump(obj, fh, separators=(",", ":"), ensure_ascii=False)
    tmp.replace(path)


def normalize_species_folder(species: str) -> str:
    text = "".join(ch.lower() if ch.isalnum() else "_" for ch in str(species or "").strip())
    while "__" in text:
        text = text.replace("__", "_")
    return text.strip("_") or "unknown"


def load_species_filters(selected_orgs: Sequence[str]) -> Optional[Set[str]]:
    selected = {
        str(item or "").strip().lower()
        for item in list(selected_orgs or [])
        if str(item or "").strip()
    }
    if not selected:
        return None
    if not SPECIES_REF_FILE.exists():
        raise FileNotFoundError(f"species_ref_list.csv not found: {SPECIES_REF_FILE}")
    allowed_species: Set[str] = set()
    with SPECIES_REF_FILE.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for raw in reader:
            row = {
                str(key or "").lstrip("\ufeff").strip(): str(value or "").strip()
                for key, value in dict(raw or {}).items()
            }
            org = row.get("Kegg Gene ID", "").strip().lower()
            species = row.get("Species", "").strip()
            if org in selected and species:
                allowed_species.add(normalize_species_folder(species))
    missing = selected - {
        str(item or "").strip().lower()
        for item in list(selected_orgs or [])
        if normalize_species_folder(str(item or "").strip()) in allowed_species
    }
    if not allowed_species:
        raise ValueError(f"No species rows matched requested orgs: {sorted(selected)}")
    return allowed_species


def run_git(args: Sequence[str], cwd: Optional[Path] = None) -> str:
    completed = subprocess.run(
        list(args),
        cwd=str(cwd) if cwd else None,
        check=True,
        capture_output=True,
        text=True,
    )
    return (completed.stdout or "").strip()


def ensure_sparse_repo(repo_dir: Path, repo_url: str, skip_fetch: bool) -> None:
    if (repo_dir / ".git").exists():
        run_git(["git", "-C", str(repo_dir), "remote", "set-url", "origin", repo_url])
        run_git(["git", "-C", str(repo_dir), "sparse-checkout", "set", "pathways"])
        if not skip_fetch:
            LOGGER.info("Updating %s", repo_dir.name)
            run_git(["git", "-C", str(repo_dir), "pull", "--ff-only", "origin", "main"])
        return
    ensure_dir(repo_dir.parent)
    LOGGER.info("Cloning %s", repo_url)
    run_git(["git", "clone", "--depth", "1", "--filter=blob:none", "--sparse", repo_url, str(repo_dir)])
    run_git(["git", "-C", str(repo_dir), "sparse-checkout", "set", "pathways"])


def get_repo_head(repo_dir: Path) -> str:
    return run_git(["git", "-C", str(repo_dir), "rev-parse", "HEAD"])


def iter_pathway_dirs(repo_dir: Path) -> Iterable[Path]:
    pathways_dir = repo_dir / "pathways"
    if not pathways_dir.exists():
        return []
    return sorted(path for path in pathways_dir.iterdir() if path.is_dir())


def read_json_file(path: Path) -> Dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return payload if isinstance(payload, dict) else {}


def extract_species_from_gpml(path: Path) -> str:
    if not path.exists():
        return ""
    try:
        root = ET.fromstring(path.read_text(encoding="utf-8", errors="replace"))
    except Exception:
        return ""
    for key in ("Organism", "organism", "Species", "species"):
        value = root.attrib.get(key)
        if str(value or "").strip():
            return str(value).strip()
    return ""


def infer_species_label(pathway_id: str, pathway_dir: Path) -> str:
    info_path = pathway_dir / f"{pathway_id}-info.json"
    info = read_json_file(info_path)
    organisms = info.get("organisms")
    if isinstance(organisms, list):
        for raw in organisms:
            if str(raw or "").strip():
                return str(raw).strip()
    gpml_path = pathway_dir / f"{pathway_id}.gpml"
    return extract_species_from_gpml(gpml_path)


def remove_stale_pathway_files(output_root: Path, pathway_id: str, keep_dir: Path) -> int:
    removed = 0
    matches = set(output_root.rglob(f"{pathway_id}.*")) | set(output_root.rglob(f"{pathway_id}-*"))
    for path in matches:
        if not path.is_file():
            continue
        if path.parent.resolve() == keep_dir.resolve():
            continue
        path.unlink(missing_ok=True)
        removed += 1
    return removed


def copy_files(src_dir: Path, dest_dir: Path, overwrite: bool) -> Tuple[int, int]:
    copied = 0
    skipped = 0
    if not src_dir.exists():
        return copied, skipped
    ensure_dir(dest_dir)
    for src in sorted(src_dir.iterdir()):
        if not src.is_file():
            continue
        dest = dest_dir / src.name
        if dest.exists() and not overwrite:
            skipped += 1
            continue
        shutil.copy2(src, dest)
        copied += 1
    return copied, skipped


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, str(args.log_level).upper()),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )

    repo_cache_dir = Path(args.repo_cache_dir)
    output_dir = Path(args.out_dir)
    manifest_path = Path(args.manifest)
    allowed_species = load_species_filters(args.org)

    database_repo_dir = repo_cache_dir / "wikipathways-database"
    assets_repo_dir = repo_cache_dir / "wikipathways-assets"
    ensure_sparse_repo(database_repo_dir, DATABASE_REPO_URL, skip_fetch=bool(args.skip_fetch))
    ensure_sparse_repo(assets_repo_dir, ASSETS_REPO_URL, skip_fetch=bool(args.skip_fetch))

    db_head = get_repo_head(database_repo_dir)
    assets_head = get_repo_head(assets_repo_dir)

    ensure_dir(output_dir)
    stats: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    processed = 0
    skipped_unknown_species = 0

    for pathway_dir in iter_pathway_dirs(database_repo_dir):
        pathway_id = pathway_dir.name.strip()
        if not pathway_id.upper().startswith("WP"):
            continue
        species_label = infer_species_label(pathway_id, pathway_dir)
        if not species_label:
            skipped_unknown_species += 1
            LOGGER.warning("Skipping %s: unable to infer species", pathway_id)
            continue
        species_folder = normalize_species_folder(species_label)
        if allowed_species is not None and species_folder not in allowed_species:
            continue
        target_dir = output_dir / species_folder
        asset_dir = assets_repo_dir / "pathways" / pathway_id

        stale_removed = remove_stale_pathway_files(output_dir, pathway_id, target_dir)
        db_copied, db_skipped = copy_files(pathway_dir, target_dir, overwrite=True)
        asset_copied, asset_skipped = copy_files(asset_dir, target_dir, overwrite=False)

        stats[species_folder]["pathways"] += 1
        stats[species_folder]["db_files_copied"] += db_copied
        stats[species_folder]["db_files_skipped"] += db_skipped
        stats[species_folder]["asset_files_copied"] += asset_copied
        stats[species_folder]["asset_files_skipped"] += asset_skipped
        stats[species_folder]["stale_files_removed"] += stale_removed
        processed += 1

        if args.max_pathways is not None and processed >= int(args.max_pathways):
            break
        if processed % 250 == 0:
            LOGGER.info("Synced %s pathways so far", processed)

    summary = {
        "pathways_processed": processed,
        "species_count": len(stats),
        "skipped_unknown_species": skipped_unknown_species,
    }
    manifest = {
        "meta": {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "database_repo_url": DATABASE_REPO_URL,
            "database_repo_head": db_head,
            "assets_repo_url": ASSETS_REPO_URL,
            "assets_repo_head": assets_head,
            "repo_cache_dir": str(repo_cache_dir),
            "out_dir": str(output_dir),
            "org_filter": list(args.org or []),
            "max_pathways": args.max_pathways,
            "summary": summary,
        },
        "species": {key: dict(value) for key, value in sorted(stats.items())},
    }
    write_json_atomic(manifest_path, manifest, pretty=bool(args.pretty))
    LOGGER.info(
        "Synced %s WikiPathways pathways into %s across %s species folders",
        processed,
        output_dir,
        len(stats),
    )
    LOGGER.info("Wrote sync manifest to %s", manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
