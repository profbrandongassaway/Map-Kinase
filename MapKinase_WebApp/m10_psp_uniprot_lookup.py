from __future__ import annotations

import argparse
import csv
from datetime import UTC, datetime
import json
import re
from pathlib import Path
import time
from typing import Any, Dict, List, Optional
from urllib.parse import parse_qs, urlparse

import requests
from bs4 import BeautifulSoup, Tag


DEFAULT_USER_AGENT = "MapKinase-PSPUniProtLookup/1.0"
DEFAULT_TIMEOUT = 30
DEFAULT_REQUEST_DELAY = 1.0
DEFAULT_CACHE_DIR = Path("MapKinase_WebApp/cache/PSP_UniProt_Lookups")
DEFAULT_CST_PATHWAY_DIR = Path("MapKinase_WebApp/cache/CST_Pathways")


def _normalize_text(value: str) -> str:
    text = str(value or "")
    text = text.replace("\xa0", " ")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _extract_psp_id(value: str | int) -> str:
    text = str(value).strip()
    if not text:
        return ""
    if text.isdigit():
        return text
    parsed = urlparse(text)
    if parsed.query:
        params = parse_qs(parsed.query)
        ids = params.get("id") or params.get("ID") or []
        if ids:
            candidate = str(ids[0]).strip()
            if candidate:
                return candidate
    match = re.search(r"(?:^|[?&])id=(\d+)", text)
    if match:
        return match.group(1)
    return text if text.isdigit() else ""


def _overview_url_for_psp_id(psp_id: str) -> str:
    return f"https://www.phosphosite.org/overviewExecuteAction?id={psp_id}&showAllSites=true"


def _cache_file_for_psp_id(psp_id: str, cache_dir: Path) -> Path:
    return cache_dir / f"{psp_id}.json"


def fetch_overview_html(psp_id_or_url: str | int, timeout: int = DEFAULT_TIMEOUT) -> tuple[str, str]:
    psp_id = _extract_psp_id(psp_id_or_url)
    if not psp_id:
        raise ValueError(f"Could not extract a PSP protein id from: {psp_id_or_url}")
    url = _overview_url_for_psp_id(psp_id)
    response = requests.get(url, headers={"User-Agent": DEFAULT_USER_AGENT}, timeout=timeout)
    response.raise_for_status()
    response.encoding = response.apparent_encoding or "utf-8"
    return psp_id, response.text


def _parse_label_value_rows(soup: BeautifulSoup) -> Dict[str, str]:
    values: Dict[str, str] = {}
    for td in soup.find_all("td"):
        label_span = td.find("span", class_=lambda c: c and "bold02" in c)
        if label_span is None:
            continue
        label = _normalize_text(label_span.get_text(" ", strip=True)).rstrip(":")
        if not label:
            continue
        clone = BeautifulSoup(str(td), "html.parser")
        label_clone = clone.find("span", class_=lambda c: c and "bold02" in c)
        if label_clone is not None:
            label_clone.extract()
        value = _normalize_text(clone.get_text(" ", strip=True))
        if value:
            values[label] = value
    return values


def _find_reference_uniprot(soup: BeautifulSoup) -> tuple[str, str]:
    for td in soup.find_all("td"):
        label_span = td.find("span", class_=lambda c: c and "bold02" in c)
        if label_span is None:
            continue
        label = _normalize_text(label_span.get_text(" ", strip=True)).rstrip(":").lower()
        if label != "reference #":
            continue
        link = td.find("a", href=True)
        accession = _normalize_text(link.get_text(" ", strip=True)) if link else ""
        href = str(link["href"]).strip() if link else ""
        if not accession and href:
            match = re.search(r"/uniprot/([A-Z0-9]+)", href, re.I)
            if match:
                accession = match.group(1).upper()
        return accession, href
    return "", ""


def _find_page_title_gene_symbol(soup: BeautifulSoup) -> str:
    title = soup.find("title")
    if not title:
        return ""
    text = _normalize_text(title.get_text(" ", strip=True))
    match = re.match(r"([A-Za-z0-9._/\-]+)\s+\(", text)
    return match.group(1) if match else text


def lookup_psp_uniprot(
    psp_id_or_url: str | int,
    timeout: int = DEFAULT_TIMEOUT,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
    request_delay: float = DEFAULT_REQUEST_DELAY,
) -> Dict[str, Any]:
    cache_root = Path(cache_dir) if cache_dir is not None else DEFAULT_CACHE_DIR
    psp_id = _extract_psp_id(psp_id_or_url)
    if not psp_id:
        raise ValueError(f"Could not extract a PSP protein id from: {psp_id_or_url}")
    cache_file = _cache_file_for_psp_id(psp_id, cache_root)
    if use_cache and cache_file.exists():
        try:
            return json.loads(cache_file.read_text(encoding="utf-8"))
        except Exception:
            pass
    if request_delay and request_delay > 0:
        time.sleep(float(request_delay))
    _, html = fetch_overview_html(psp_id, timeout=timeout)
    soup = BeautifulSoup(html, "html.parser")
    field_map = _parse_label_value_rows(soup)
    accession, accession_url = _find_reference_uniprot(soup)
    gene_symbol = field_map.get("Gene Symbols", "") or _find_page_title_gene_symbol(soup)
    protein_name = ""
    first_bold = soup.find("span", class_=lambda c: c and "bold01" in c)
    if first_bold and isinstance(first_bold.parent, Tag):
        parent_text = _normalize_text(first_bold.parent.get_text(" ", strip=True))
        head = _normalize_text(first_bold.get_text(" ", strip=True))
        protein_name = parent_text[len(head):].strip(" -,:") if parent_text.startswith(head) else parent_text
    result = {
        "psp_id": psp_id,
        "overview_url": _overview_url_for_psp_id(psp_id),
        "gene_symbol": gene_symbol,
        "uniprot_id": accession,
        "uniprot_url": accession_url,
        "protein_name": protein_name,
        "organism": field_map.get("Organism", ""),
        "reference_label": field_map.get("Reference #", ""),
        "raw_fields": field_map,
    }
    if use_cache:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        cache_file.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    return result


def lookup_many_psp_uniprot(
    items: List[str],
    timeout: int = DEFAULT_TIMEOUT,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
    request_delay: float = DEFAULT_REQUEST_DELAY,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for item in items:
        text = str(item or "").strip()
        if not text:
            continue
        try:
            rows.append(
                lookup_psp_uniprot(
                    text,
                    timeout=timeout,
                    cache_dir=cache_dir,
                    use_cache=use_cache,
                    request_delay=request_delay,
                )
            )
        except Exception as exc:
            rows.append(
                {
                    "psp_id": _extract_psp_id(text),
                    "input": text,
                    "error": str(exc),
                }
            )
    return rows


def _collect_psp_ids_from_module(module: Dict[str, Any]) -> List[str]:
    ids = list(module.get("phosphositeplus_ids") or [])
    if not ids:
        raw = str(module.get("phosphositeplus_id") or "").strip()
        if raw:
            ids = [item for item in re.split(r"[;,]\s*", raw) if item]
    cleaned: List[str] = []
    seen: set[str] = set()
    for item in ids:
        psp_id = _extract_psp_id(item)
        if not psp_id or psp_id in seen:
            continue
        seen.add(psp_id)
        cleaned.append(psp_id)
    return cleaned


def enrich_protein_module_with_uniprot(
    module: Dict[str, Any],
    pathway_name: str,
    timeout: int = DEFAULT_TIMEOUT,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
    request_delay: float = DEFAULT_REQUEST_DELAY,
) -> Dict[str, Any]:
    if module.get("has_uniprot_match") or list(module.get("uniprot_ids") or []):
        return module
    psp_ids = _collect_psp_ids_from_module(module)
    lookup_rows: List[Dict[str, Any]] = []
    uniprot_ids: List[str] = []
    uniprot_urls: List[str] = []
    gene_symbols: List[str] = []
    seen_uniprots: set[str] = set()
    seen_urls: set[str] = set()
    seen_genes: set[str] = set()

    for psp_id in psp_ids:
        search_url = _overview_url_for_psp_id(psp_id)
        print(f"Pathway: {pathway_name}")
        print(f"Searching PSP link: {search_url}")
        try:
            row = lookup_psp_uniprot(
                psp_id,
                timeout=timeout,
                cache_dir=cache_dir,
                use_cache=use_cache,
                request_delay=request_delay,
            )
        except Exception as exc:
            print(f"PSP lookup error for {psp_id}: {exc}")
            row = {
                "psp_id": psp_id,
                "overview_url": search_url,
                "error": str(exc),
            }
        lookup_rows.append(row)
        uniprot_id = str(row.get("uniprot_id") or "").strip()
        uniprot_url = str(row.get("uniprot_url") or "").strip()
        gene_symbol = str(row.get("gene_symbol") or "").strip()
        if uniprot_id and uniprot_id not in seen_uniprots:
            seen_uniprots.add(uniprot_id)
            uniprot_ids.append(uniprot_id)
        if uniprot_url and uniprot_url not in seen_urls:
            seen_urls.add(uniprot_url)
            uniprot_urls.append(uniprot_url)
        if gene_symbol and gene_symbol not in seen_genes:
            seen_genes.add(gene_symbol)
            gene_symbols.append(gene_symbol)

    module["psp_uniprot_lookups"] = lookup_rows
    module["uniprot_ids"] = uniprot_ids
    module["uniprot_urls"] = uniprot_urls
    module["psp_gene_symbols"] = gene_symbols
    module["has_uniprot_match"] = bool(uniprot_ids)
    return module


def enrich_cst_pathway_file(
    pathway_file: Path,
    timeout: int = DEFAULT_TIMEOUT,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
    request_delay: float = DEFAULT_REQUEST_DELAY,
) -> Dict[str, Any]:
    obj = json.loads(pathway_file.read_text(encoding="utf-8"))
    pathway_name = _normalize_text(str(obj.get("pathway_name") or pathway_file.stem))
    modules = list(obj.get("protein_modules") or [])
    for module in modules:
        enrich_protein_module_with_uniprot(
            module,
            pathway_name=pathway_name,
            timeout=timeout,
            cache_dir=cache_dir,
            use_cache=use_cache,
            request_delay=request_delay,
        )
    obj["protein_modules"] = modules
    obj["uniprot_enriched_at"] = datetime.now(UTC).isoformat()
    pathway_file.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    return {
        "pathway_name": pathway_name,
        "pathway_file": str(pathway_file.resolve()),
        "protein_module_count": len(modules),
        "matched_modules": sum(1 for module in modules if module.get("has_uniprot_match")),
    }


def enrich_cst_pathway_directory(
    input_dir: Path,
    timeout: int = DEFAULT_TIMEOUT,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
    limit: int = 0,
    request_delay: float = DEFAULT_REQUEST_DELAY,
) -> Dict[str, Any]:
    files = sorted(
        [
            path
            for path in input_dir.glob("*.json")
            if path.is_file() and path.name.lower() != "_manifest.json"
        ],
        key=lambda p: p.name.lower(),
    )
    if limit and limit > 0:
        files = files[:limit]
    rows: List[Dict[str, Any]] = []
    for idx, pathway_file in enumerate(files, start=1):
        print(f"=== [{idx}/{len(files)}] {pathway_file.name} ===")
        rows.append(
            enrich_cst_pathway_file(
                pathway_file,
                timeout=timeout,
                cache_dir=cache_dir,
                use_cache=use_cache,
                request_delay=request_delay,
            )
        )
    return {
        "input_dir": str(input_dir.resolve()),
        "pathway_files_processed": len(rows),
        "pathways": rows,
    }


def export_results_to_json(result: Any, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")


def export_results_to_csv(results: List[Dict[str, Any]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["input", "psp_id", "gene_symbol", "uniprot_id", "uniprot_url", "protein_name", "organism", "overview_url", "error"]
    with output_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Look up UniProtKB accessions from PhosphoSitePlus protein pages.")
    parser.add_argument(
        "--mode",
        choices=["lookup", "enrich-cst"],
        default="lookup",
        help="Look up one or more PSP ids/URLs, or enrich CST pathway JSON files in place.",
    )
    parser.add_argument("items", nargs="*", help="PSP protein ids or proteinAction/overview URLs.")
    parser.add_argument("--output", default="", help="Optional output file (.json or .csv).")
    parser.add_argument("--timeout", type=int, default=DEFAULT_TIMEOUT, help="HTTP timeout in seconds.")
    parser.add_argument("--delay", type=float, default=DEFAULT_REQUEST_DELAY, help="Delay in seconds before each uncached PSP request.")
    parser.add_argument("--input-dir", default=str(DEFAULT_CST_PATHWAY_DIR), help="CST pathway JSON directory for enrich-cst mode.")
    parser.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR), help="Cache directory for PSP UniProt lookups.")
    parser.add_argument("--limit", type=int, default=0, help="Optional limit for enrich-cst mode.")
    parser.add_argument("--no-cache", action="store_true", help="Disable disk cache for PSP UniProt lookups.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "enrich-cst":
        result = enrich_cst_pathway_directory(
            input_dir=Path(args.input_dir),
            timeout=int(args.timeout),
            cache_dir=Path(args.cache_dir),
            use_cache=not bool(args.no_cache),
            limit=int(args.limit),
            request_delay=float(args.delay),
        )
        if args.output:
            export_results_to_json(result, Path(args.output))
        else:
            print(json.dumps(result, indent=2, ensure_ascii=False))
    else:
        if not args.items:
            raise SystemExit("Provide at least one PSP id or PSP protein URL.")
        results = lookup_many_psp_uniprot(
            [str(item) for item in args.items],
            timeout=int(args.timeout),
            cache_dir=Path(args.cache_dir),
            use_cache=not bool(args.no_cache),
            request_delay=float(args.delay),
        )
        if args.output:
            output_path = Path(args.output)
            if output_path.suffix.lower() == ".csv":
                export_results_to_csv(results, output_path)
            else:
                export_results_to_json(results, output_path)
        else:
            print(json.dumps(results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
