from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, UTC
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional
from urllib.parse import urljoin, urlparse

import requests
from bs4 import BeautifulSoup, Tag


DEFAULT_URL = "https://www.cellsignal.com/pathways/regulation-of-apoptosis-pathway"
DEFAULT_RESEARCH_AREA_URL = "https://www.cellsignal.com/pathways/by-research-area"
DEFAULT_USER_AGENT = "MapKinase-CSTWebpageScraper/1.0"
PHOSPHOSITE_PROTEIN_URL = "https://www.phosphosite.org/proteinAction.action?id={psp_id}"
DEFAULT_PATHWAY_LINKS_FILE = Path("MapKinase_WebApp/cache/CST_pathway_links.json")
DEFAULT_PATHWAY_OUTPUT_DIR = Path("MapKinase_WebApp/cache/CST_Pathways")


def _normalize_text(value: str) -> str:
    text = str(value or "")
    text = text.replace("\xa0", " ")
    text = text.replace("\u2013", "-").replace("\u2014", "-").replace("\u2212", "-")
    text = re.sub(r"\s+", " ", text).strip()
    return text


def _dedupe_preserve_order(values: Iterable[str]) -> List[str]:
    output: List[str] = []
    seen: set[str] = set()
    for value in values:
        text = _normalize_text(value)
        if not text or text in seen:
            continue
        seen.add(text)
        output.append(text)
    return output


def _slugify_filename(value: str) -> str:
    text = _normalize_text(value).lower()
    text = re.sub(r"[^a-z0-9]+", "_", text).strip("_")
    return text or "cst_pathway"


def fetch_html(url: str, timeout: int = 30) -> str:
    headers = {"User-Agent": DEFAULT_USER_AGENT}
    response = requests.get(url, headers=headers, timeout=timeout)
    response.raise_for_status()
    if not response.encoding or response.encoding.lower() == "iso-8859-1":
        response.encoding = response.apparent_encoding or "utf-8"
    return response.text


def _find_text_parent(soup: BeautifulSoup, snippet: str) -> Optional[Tag]:
    found = soup.find(string=lambda s: isinstance(s, str) and snippet in s)
    return found.parent if found and isinstance(found.parent, Tag) else None


def extract_breadcrumbs(soup: BeautifulSoup) -> List[str]:
    crumbs: List[str] = []
    for nav in soup.find_all(["nav", "ol", "ul"]):
        text = " ".join(nav.stripped_strings)
        if "Home" in text and "Pathways" in text:
            crumbs = [item for item in _dedupe_preserve_order(nav.stripped_strings) if item != "/"]
            if crumbs:
                return crumbs
    return crumbs


def extract_media_links(soup: BeautifulSoup, page_url: str) -> Dict[str, str]:
    output: Dict[str, str] = {}
    for label in ["DOWNLOAD PATHWAY", "View PDF", "Request Permission for Pathway"]:
        node = soup.find("a", string=lambda s: isinstance(s, str) and _normalize_text(s) == label)
        if node and node.get("href"):
            output[label] = urljoin(page_url, str(node["href"]))
    return output


def _extract_pathway_diagram_container(soup: BeautifulSoup) -> Optional[Tag]:
    trigger = soup.find("button", string=lambda s: isinstance(s, str) and "VIEW INTERACTIVE PATHWAY" in s)
    if not trigger:
        return None
    current: Optional[Tag] = trigger
    while current is not None:
        if (
            current.name == "div"
            and "tw-relative" in (current.get("class") or [])
            and "tw-block" in (current.get("class") or [])
            and "tw-overflow-scroll" in (current.get("class") or [])
        ):
            return current
        current = current.parent if isinstance(current.parent, Tag) else None
    return None


def extract_pathway_diagram_data(soup: BeautifulSoup, page_url: str) -> Dict[str, Any]:
    container = _extract_pathway_diagram_container(soup)
    if not container:
        return {"diagram_labels": [], "diagram_links": [], "diagram_proteins": []}

    link_rows: List[Dict[str, str]] = []
    for anchor in container.find_all("a", href=True):
        label = _normalize_text(anchor.get_text(" ", strip=True))
        href = urljoin(page_url, str(anchor["href"]))
        if not label:
            continue
        link_rows.append({"label": label, "url": href})

    protein_rows: List[Dict[str, str]] = []
    seen_protein_keys: set[tuple[str, str, str]] = set()
    for anchor in container.find_all("a", href=True):
        text_node = next(
            (
                child
                for child in anchor.children
                if isinstance(child, Tag) and str(getattr(child, "name", "")).lower() == "text"
            ),
            None,
        )
        if text_node is None:
            continue
        label = _normalize_text(text_node.get_text(" ", strip=True))
        cst_url = urljoin(page_url, str(anchor["href"]))
        psp_id = _normalize_text(str(text_node.get("data-psp") or ""))
        if not label:
            continue
        psp_ids = [
            item
            for item in re.split(r"[;,]\s*", psp_id)
            if item and item.lower() not in {"undefined", "null", "none", "nan"}
        ]
        psp_urls = [PHOSPHOSITE_PROTEIN_URL.format(psp_id=item) for item in psp_ids]
        key = (label, cst_url, "|".join(psp_urls))
        if key in seen_protein_keys:
            continue
        seen_protein_keys.add(key)
        protein_rows.append(
            {
                "label": label,
                "cst_url": cst_url,
                "phosphositeplus_id": psp_id if psp_ids else "",
                "phosphositeplus_ids": psp_ids,
                "phosphositeplus_url": psp_urls[0] if len(psp_urls) == 1 else "",
                "phosphositeplus_urls": psp_urls,
            }
        )

    raw_strings = [
        _normalize_text(text)
        for text in container.stripped_strings
        if _normalize_text(text) and _normalize_text(text) != "VIEW INTERACTIVE PATHWAY"
    ]
    labels = _dedupe_preserve_order(raw_strings)
    return {
        "diagram_labels": labels,
        "diagram_links": link_rows,
        "diagram_proteins": protein_rows,
    }


def extract_description_block(soup: BeautifulSoup, page_url: str) -> Dict[str, Any]:
    output: Dict[str, Any] = {
        "description_paragraphs": [],
        "selected_reviews": [],
        "acknowledgements": [],
        "created": "",
        "revised": "",
    }

    desc_button = _find_text_parent(soup, "Pathway description")
    if not desc_button:
        return output

    root = desc_button
    for _ in range(6):
        if root is None:
            break
        if isinstance(root, Tag) and root.find("h2", string=lambda s: isinstance(s, str) and "Selected Reviews" in s):
            break
        root = root.parent if isinstance(root.parent, Tag) else root

    search_root = root if isinstance(root, Tag) else soup
    review_h2 = search_root.find("h2", string=lambda s: isinstance(s, str) and "Selected Reviews" in s)
    if not review_h2:
        review_h2 = soup.find("h2", string=lambda s: isinstance(s, str) and "Selected Reviews" in s)

    paragraphs: List[str] = []
    if review_h2:
        current = review_h2.find_previous_sibling()
        collected: List[Tag] = []
        while current is not None:
            if isinstance(current, Tag) and current.name == "p":
                collected.append(current)
            current = current.find_previous_sibling()
        paragraphs = [_normalize_text(node.get_text(" ", strip=True)) for node in reversed(collected)]
    output["description_paragraphs"] = [p for p in paragraphs if p]

    if review_h2:
        ul = review_h2.find_next("ul")
        if ul:
            reviews: List[Dict[str, str]] = []
            for li in ul.find_all("li", recursive=False):
                citation = _normalize_text(li.get_text(" ", strip=True))
                link = li.find("a", href=True)
                reviews.append(
                    {
                        "citation": citation,
                        "url": urljoin(page_url, str(link["href"])) if link else "",
                    }
                )
            output["selected_reviews"] = reviews

        tail_strings: List[str] = []
        cursor = ul.find_next_sibling() if ul else review_h2.find_next_sibling()
        while cursor is not None:
            if isinstance(cursor, Tag) and cursor.name in {"h2", "footer"}:
                break
            text = _normalize_text(cursor.get_text(" ", strip=True)) if isinstance(cursor, Tag) else ""
            if text:
                tail_strings.append(text)
            cursor = cursor.find_next_sibling() if isinstance(cursor, Tag) else None

        acknowledgements: List[str] = []
        for item in tail_strings:
            lower = item.lower()
            if lower.startswith("created "):
                output["created"] = item.removeprefix("created ").strip()
            elif lower.startswith("revised "):
                output["revised"] = item.removeprefix("revised ").strip()
            elif "thank" in lower or "reviewing this diagram" in lower:
                acknowledgements.append(item)
        output["acknowledgements"] = acknowledgements

    return output


def scrape_cst_pathway_page(url: str) -> Dict[str, Any]:
    html = fetch_html(url)
    soup = BeautifulSoup(html, "html.parser")

    title_tag = soup.find("h1")
    page_title = _normalize_text(title_tag.get_text(" ", strip=True)) if title_tag else ""
    meta_description = ""
    meta = soup.find("meta", attrs={"name": "description"}) or soup.find("meta", attrs={"property": "og:description"})
    if meta and meta.get("content"):
        meta_description = _normalize_text(str(meta["content"]))

    subtitle = ""
    access_line = soup.find(string=lambda s: isinstance(s, str) and "Access the full library of downloadable pathway diagrams" in s)
    if access_line:
        subtitle = _normalize_text(str(access_line))

    result: Dict[str, Any] = {
        "url": url,
        "scraped_at_utc": datetime.now(UTC).isoformat(),
        "page_title": page_title,
        "meta_title": _normalize_text(soup.title.get_text(" ", strip=True)) if soup.title else "",
        "meta_description": meta_description,
        "subtitle": subtitle,
        "breadcrumbs": extract_breadcrumbs(soup),
        "media_links": extract_media_links(soup, url),
    }
    result.update(extract_pathway_diagram_data(soup, url))
    result.update(extract_description_block(soup, url))
    result["raw_page_text"] = _normalize_text(soup.get_text(" ", strip=True))
    return result


def extract_pathway_protein_modules(result: Dict[str, Any], pathway_meta: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    meta = dict(pathway_meta or {})
    pathway_name = _normalize_text(
        str(result.get("page_title") or meta.get("label") or meta.get("name") or "")
    )
    source_url = str(result.get("url") or meta.get("url") or "").strip()
    return {
        "pathway_name": pathway_name,
        "pathway_url": source_url,
        "source_sections": list(meta.get("sections") or []),
        "protein_modules": list(result.get("diagram_proteins") or []),
        "protein_module_count": len(result.get("diagram_proteins") or []),
        "scraped_at_utc": str(result.get("scraped_at_utc") or datetime.now(UTC).isoformat()),
    }


def scrape_cst_research_area_catalog(url: str) -> Dict[str, Any]:
    html = fetch_html(url)
    soup = BeautifulSoup(html, "html.parser")

    page_title = ""
    title_tag = soup.find("h1")
    if title_tag:
        page_title = _normalize_text(title_tag.get_text(" ", strip=True))

    section_rows: List[Dict[str, str]] = []
    pathway_map: Dict[str, Dict[str, Any]] = {}

    for accordion in soup.find_all("div", class_="context-accordion"):
        heading = accordion.find("h3")
        body = accordion.find("div", class_="context-body")
        if heading is None or body is None:
            continue
        current_section = _normalize_text(heading.get_text(" ", strip=True))
        if not current_section:
            continue
        for node in body.find_all("a", href=True):
            href = str(node["href"])
            label = _normalize_text(node.get_text(" ", strip=True))
            full_url = urljoin(url, href)
            if not label or "/pathways/" not in href or "/pathways/by-" in href:
                continue
            section_rows.append(
                {
                    "section": current_section,
                    "label": label,
                    "url": full_url,
                }
            )
            existing = pathway_map.get(full_url)
            if existing is None:
                pathway_map[full_url] = {
                    "label": label,
                    "url": full_url,
                    "sections": [current_section],
                }
            elif current_section not in existing["sections"]:
                existing["sections"].append(current_section)

    return {
        "url": url,
        "scraped_at_utc": datetime.now(UTC).isoformat(),
        "page_title": page_title,
        "breadcrumbs": extract_breadcrumbs(soup),
        "pathway_links": list(pathway_map.values()),
        "pathway_links_by_section": section_rows,
        "total_unique_pathway_links": len(pathway_map),
        "total_section_rows": len(section_rows),
    }


def scrape_cst_pathway_links_file(
    input_file: Path,
    output_dir: Path,
    limit: Optional[int] = None,
) -> Dict[str, Any]:
    obj = json.loads(input_file.read_text(encoding="utf-8"))
    pathway_links = list(obj.get("pathway_links") or [])
    if limit is not None and limit > 0:
        pathway_links = pathway_links[:limit]

    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows: List[Dict[str, Any]] = []
    used_names: set[str] = set()

    for idx, pathway in enumerate(pathway_links, start=1):
        url = str(pathway.get("url") or "").strip()
        if not url:
            continue
        result = scrape_cst_pathway_page(url)
        payload = extract_pathway_protein_modules(result, pathway)
        pathway_name = str(payload.get("pathway_name") or pathway.get("label") or f"pathway_{idx}")
        stem = _slugify_filename(pathway_name)
        if stem in used_names:
            suffix = _slugify_filename(Path(urlparse(url).path).name.replace(".html", "")) or f"pathway_{idx}"
            stem = f"{stem}__{suffix}"
        used_names.add(stem)
        output_path = output_dir / f"{stem}.json"
        export_results_to_json(payload, output_path)
        manifest_rows.append(
            {
                "pathway_name": payload.get("pathway_name"),
                "pathway_url": payload.get("pathway_url"),
                "output_file": str(output_path.resolve()),
                "protein_module_count": payload.get("protein_module_count", 0),
                "source_sections": payload.get("source_sections", []),
            }
        )
        print(f"[{idx}/{len(pathway_links)}] {payload.get('pathway_name')} -> {output_path.name} ({payload.get('protein_module_count', 0)} proteins)")

    summary = {
        "input_file": str(input_file.resolve()),
        "output_dir": str(output_dir.resolve()),
        "scraped_at_utc": datetime.now(UTC).isoformat(),
        "total_pathways_requested": len(pathway_links),
        "total_pathways_written": len(manifest_rows),
        "pathways": manifest_rows,
    }
    export_results_to_json(summary, output_dir / "_manifest.json")
    return summary


def export_results_to_json(result: Dict[str, Any], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Scrape a Cell Signaling Technology pathway webpage into structured JSON.")
    parser.add_argument("--url", default=DEFAULT_URL, help="CST page URL to scrape.")
    parser.add_argument(
        "--mode",
        choices=["pathway", "research-area", "pathway-batch"],
        default="pathway",
        help="Scrape a single pathway page, the research-area catalog page, or a batch of pathway pages from a links file.",
    )
    parser.add_argument(
        "--output",
        default="MapKinase_WebApp/cache/cst_scrape_regulation_of_apoptosis.json",
        help="Output JSON path.",
    )
    parser.add_argument(
        "--print-summary",
        action="store_true",
        help="Print a small summary after scraping.",
    )
    parser.add_argument(
        "--input-file",
        default=str(DEFAULT_PATHWAY_LINKS_FILE),
        help="Input JSON file containing pathway_links for batch mode.",
    )
    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_PATHWAY_OUTPUT_DIR),
        help="Output directory for per-pathway files in batch mode.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional number of pathway links to scrape in batch mode.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "research-area":
        result = scrape_cst_research_area_catalog(str(args.url or DEFAULT_RESEARCH_AREA_URL))
        export_results_to_json(result, Path(args.output))
    elif args.mode == "pathway-batch":
        result = scrape_cst_pathway_links_file(
            input_file=Path(args.input_file),
            output_dir=Path(args.output_dir),
            limit=(args.limit if args.limit and args.limit > 0 else None),
        )
    else:
        result = scrape_cst_pathway_page(str(args.url))
        export_results_to_json(result, Path(args.output))
    if args.print_summary:
        if args.mode == "research-area":
            print(f"Title: {result.get('page_title')}")
            print(f"Unique pathway links: {result.get('total_unique_pathway_links')}")
            print(f"Section rows: {result.get('total_section_rows')}")
            print(f"Output: {Path(args.output).resolve()}")
        elif args.mode == "pathway-batch":
            print(f"Pathways written: {result.get('total_pathways_written')}")
            print(f"Output dir: {result.get('output_dir')}")
        else:
            print(f"Title: {result.get('page_title')}")
            print(f"Diagram labels: {len(result.get('diagram_labels') or [])}")
            print(f"Diagram links: {len(result.get('diagram_links') or [])}")
            print(f"Reviews: {len(result.get('selected_reviews') or [])}")
            print(f"Output: {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()
