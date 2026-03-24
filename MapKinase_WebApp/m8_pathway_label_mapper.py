from __future__ import annotations

import csv
import json
import logging
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

try:
    import requests  # type: ignore
except ImportError:  # pragma: no cover
    requests = None


LOGGER = logging.getLogger(__name__)

BASE_DIR = Path(__file__).resolve().parent
ANNOTATION_DIR = BASE_DIR / "annotation_files"
CACHE_DIR = BASE_DIR / "cache" / "pathway_label_mapper"
LOCAL_UNIPROT_ANNOTATIONS = ANNOTATION_DIR / "uniprot_human_function_annotations.tsv"
DEFAULT_ORGANISM_ID = 9606
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


SIGNALING_LABEL_RULES: Dict[str, Dict[str, Any]] = {
    "P38": {
        "mapping_type": "alias",
        "gene_symbols": ["MAPK14"],
        "notes": "Common shorthand for p38 MAPK in human signaling diagrams; mapped conservatively to MAPK14 (p38alpha).",
    },
    "P38MAPK": {
        "mapping_type": "alias",
        "gene_symbols": ["MAPK14"],
        "notes": "Common shorthand for p38 MAPK in human signaling diagrams; mapped conservatively to MAPK14 (p38alpha).",
    },
    "ERK": {
        "mapping_type": "family",
        "gene_symbols": ["MAPK1", "MAPK3"],
        "notes": "ERK usually refers to the ERK1/2 family in human signaling pathways.",
    },
    "ERK1/2": {
        "mapping_type": "family",
        "gene_symbols": ["MAPK3", "MAPK1"],
        "notes": "ERK1/2 expands to MAPK3 (ERK1) and MAPK1 (ERK2).",
    },
    "P44/42MAPK": {
        "mapping_type": "family",
        "gene_symbols": ["MAPK3", "MAPK1"],
        "notes": "p44/42 MAPK is common pathway shorthand for ERK1/2.",
    },
    "JNK": {
        "mapping_type": "family",
        "gene_symbols": ["MAPK8", "MAPK9", "MAPK10"],
        "notes": "JNK expands to the major human JNK family members.",
    },
    "SAPK/JNK": {
        "mapping_type": "family",
        "gene_symbols": ["MAPK8", "MAPK9", "MAPK10"],
        "notes": "SAPK/JNK is common shorthand for the JNK family.",
    },
    "SAPK": {
        "mapping_type": "family",
        "gene_symbols": ["MAPK8", "MAPK9", "MAPK10"],
        "notes": "SAPK is commonly used interchangeably with JNK in stress signaling.",
    },
    "MEK": {
        "mapping_type": "family",
        "gene_symbols": ["MAP2K1", "MAP2K2"],
        "notes": "MEK without an isoform usually refers to MEK1/2.",
    },
    "MEK1": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP2K1"],
        "notes": "MEK1 maps to MAP2K1.",
    },
    "MEK2": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP2K2"],
        "notes": "MEK2 maps to MAP2K2.",
    },
    "MEK1/2": {
        "mapping_type": "family",
        "gene_symbols": ["MAP2K1", "MAP2K2"],
        "notes": "MEK1/2 expands to MAP2K1 and MAP2K2.",
    },
    "SEK1": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP2K4"],
        "notes": "SEK1 is a common alias for MAP2K4.",
    },
    "MKK4": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP2K4"],
        "notes": "MKK4 maps to MAP2K4.",
    },
    "MKK3": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP2K3"],
        "notes": "MKK3 maps to MAP2K3.",
    },
    "MKK6": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP2K6"],
        "notes": "MKK6 maps to MAP2K6.",
    },
    "MKK3/6": {
        "mapping_type": "family",
        "gene_symbols": ["MAP2K3", "MAP2K6"],
        "notes": "MKK3/6 expands to MAP2K3 and MAP2K6.",
    },
    "ASK1": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP3K5"],
        "notes": "ASK1 maps to MAP3K5.",
    },
    "ASK1/2": {
        "mapping_type": "family",
        "gene_symbols": ["MAP3K5", "MAP3K6"],
        "notes": "ASK1/2 expands to MAP3K5 and MAP3K6.",
    },
    "TAK1": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP3K7"],
        "notes": "TAK1 maps to MAP3K7.",
    },
    "DLK": {
        "mapping_type": "alias",
        "gene_symbols": ["MAP3K12"],
        "notes": "DLK maps to MAP3K12.",
    },
    "MLK2/3": {
        "mapping_type": "family",
        "gene_symbols": ["MAP3K10", "MAP3K11"],
        "notes": "MLK2/3 expands to MAP3K10 and MAP3K11.",
    },
    "ELK-1": {
        "mapping_type": "alias",
        "gene_symbols": ["ELK1"],
        "notes": "Elk-1 maps to ELK1.",
    },
    "C-JUN": {
        "mapping_type": "alias",
        "gene_symbols": ["JUN"],
        "notes": "c-Jun maps to JUN.",
    },
    "P53": {
        "mapping_type": "alias",
        "gene_symbols": ["TP53"],
        "notes": "p53 maps to TP53.",
    },
    "AKT": {
        "mapping_type": "family",
        "gene_symbols": ["AKT1", "AKT2", "AKT3"],
        "notes": "AKT usually refers to the AKT1/2/3 family.",
    },
    "PKB": {
        "mapping_type": "family",
        "gene_symbols": ["AKT1", "AKT2", "AKT3"],
        "notes": "PKB is a common alias for the AKT family.",
    },
    "RAS": {
        "mapping_type": "family",
        "gene_symbols": ["HRAS", "KRAS", "NRAS"],
        "notes": "RAS expands to the major human HRAS/KRAS/NRAS family members.",
    },
    "RAF": {
        "mapping_type": "family",
        "gene_symbols": ["ARAF", "BRAF", "RAF1"],
        "notes": "RAF refers to the ARAF/BRAF/RAF1 kinase family.",
    },
    "PI3K": {
        "mapping_type": "family",
        "gene_symbols": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3"],
        "notes": "PI3K is a family-level label; expanded to common human class I catalytic and regulatory subunits.",
    },
    "MTOR": {
        "mapping_type": "alias",
        "gene_symbols": ["MTOR"],
        "notes": "mTOR maps to MTOR.",
    },
    "P70S6K": {
        "mapping_type": "family",
        "gene_symbols": ["RPS6KB1", "RPS6KB2"],
        "notes": "p70S6K can refer to RPS6KB1 and RPS6KB2.",
    },
    "S6K": {
        "mapping_type": "family",
        "gene_symbols": ["RPS6KB1", "RPS6KB2"],
        "notes": "S6K expands to the common p70 S6 kinase family members.",
    },
    "RSK": {
        "mapping_type": "family",
        "gene_symbols": ["RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA6"],
        "notes": "RSK is a family-level label covering major p90RSK genes.",
    },
    "P90RSK": {
        "mapping_type": "family",
        "gene_symbols": ["RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA6"],
        "notes": "p90RSK is a family-level label covering major p90RSK genes.",
    },
    "PKC": {
        "mapping_type": "family",
        "gene_symbols": ["PRKCA", "PRKCB", "PRKCG", "PRKCD", "PRKCE", "PRKCQ", "PRKCI", "PRKCZ"],
        "notes": "PKC is a family-level label and should not be collapsed to one protein.",
    },
    "PRAK": {
        "mapping_type": "alias",
        "gene_symbols": ["MAPKAPK5"],
        "notes": "PRAK maps to MAPKAPK5.",
    },
    "MAPKAPK-3": {
        "mapping_type": "alias",
        "gene_symbols": ["MAPKAPK3"],
        "notes": "MAPKAPK-3 maps to MAPKAPK3.",
    },
    "MAPKAPK-2": {
        "mapping_type": "alias",
        "gene_symbols": ["MAPKAPK2"],
        "notes": "MAPKAPK-2 maps to MAPKAPK2.",
    },
    "MK2": {
        "mapping_type": "alias",
        "gene_symbols": ["MAPKAPK2"],
        "notes": "MK2 maps to MAPKAPK2.",
    },
    "MK3": {
        "mapping_type": "alias",
        "gene_symbols": ["MAPKAPK3"],
        "notes": "MK3 maps to MAPKAPK3.",
    },
    "MSK1/2": {
        "mapping_type": "family",
        "gene_symbols": ["RPS6KA5", "RPS6KA4"],
        "notes": "MSK1/2 expands to RPS6KA5 and RPS6KA4.",
    },
    "MNK1/2": {
        "mapping_type": "family",
        "gene_symbols": ["MKNK1", "MKNK2"],
        "notes": "MNK1/2 expands to MKNK1 and MKNK2.",
    },
    "HSP27": {
        "mapping_type": "alias",
        "gene_symbols": ["HSPB1"],
        "notes": "HSP27 maps to HSPB1.",
    },
    "CPLA2": {
        "mapping_type": "alias",
        "gene_symbols": ["PLA2G4A"],
        "notes": "cPLA2 maps to PLA2G4A.",
    },
    "RHO-GDI": {
        "mapping_type": "family",
        "gene_symbols": ["ARHGDIA", "ARHGDIB", "ARHGDIG"],
        "notes": "Rho-GDI is a family-level label; expanded to the major human Rho GDP dissociation inhibitors.",
    },
    "RHOGDI": {
        "mapping_type": "family",
        "gene_symbols": ["ARHGDIA", "ARHGDIB", "ARHGDIG"],
        "notes": "Rho-GDI is a family-level label; expanded to the major human Rho GDP dissociation inhibitors.",
    },
    "HISTONEH3": {
        "mapping_type": "family",
        "gene_symbols": ["H3C1", "H3C2", "H3C3", "H3C4", "H3C6", "H3C7", "H3C8", "H3C10", "H3C11", "H3C12", "H3C13", "H3C14", "H3C15"],
        "notes": "Histone H3 is a family-level label and should not be collapsed to one human H3 gene.",
    },
    "HISTONE 3": {
        "mapping_type": "family",
        "gene_symbols": ["H3C1", "H3C2", "H3C3", "H3C4", "H3C6", "H3C7", "H3C8", "H3C10", "H3C11", "H3C12", "H3C13", "H3C14", "H3C15"],
        "notes": "Histone H3 is a family-level label and should not be collapsed to one human H3 gene.",
    },
    "ATF-2": {
        "mapping_type": "alias",
        "gene_symbols": ["ATF2"],
        "notes": "ATF-2 maps to ATF2.",
    },
    "ATF2": {
        "mapping_type": "alias",
        "gene_symbols": ["ATF2"],
        "notes": "ATF2 maps directly to ATF2.",
    },
    "CREB": {
        "mapping_type": "alias",
        "gene_symbols": ["CREB1"],
        "notes": "CREB in signaling diagrams usually refers to CREB1.",
    },
    "NFKB": {
        "mapping_type": "family",
        "gene_symbols": ["NFKB1", "NFKB2", "RELA", "RELB", "REL"],
        "notes": "NF-kB without a subunit expands to the major human NF-kB family members.",
    },
    "P65": {
        "mapping_type": "alias",
        "gene_symbols": ["RELA"],
        "notes": "NF-kB p65 maps to RELA.",
    },
    "NF-KBP65": {
        "mapping_type": "alias",
        "gene_symbols": ["RELA"],
        "notes": "NF-kB p65 maps to RELA.",
    },
    "MEKK1-4": {
        "mapping_type": "family",
        "gene_symbols": ["MAP3K1", "MAP3K2", "MAP3K3", "MAP3K4"],
        "notes": "MEKK1-4 expands to the major MEKK family members MAP3K1/2/3/4.",
    },
    "TAO1/2/3": {
        "mapping_type": "family",
        "gene_symbols": ["TAOK1", "TAOK2", "TAOK3"],
        "notes": "TAO1/2/3 expands to TAOK1, TAOK2, and TAOK3.",
    },
    "RAR": {
        "mapping_type": "family",
        "gene_symbols": ["RARA", "RARB", "RARG"],
        "notes": "RAR expands to the major human retinoic acid receptor family members.",
    },
    "RARA": {
        "mapping_type": "alias",
        "gene_symbols": ["RARA"],
        "notes": "RARalpha maps to RARA.",
    },
    "RARALPHA": {
        "mapping_type": "alias",
        "gene_symbols": ["RARA"],
        "notes": "RARalpha maps to RARA.",
    },
    "DUSPS": {
        "mapping_type": "family",
        "gene_symbols": ["DUSP1", "DUSP2", "DUSP4", "DUSP5", "DUSP6", "DUSP7", "DUSP8", "DUSP9", "DUSP10", "DUSP16"],
        "notes": "DUSPs is a family-level label and should remain expanded.",
    },
    "IRS-1": {
        "mapping_type": "alias",
        "gene_symbols": ["IRS1"],
        "notes": "IRS-1 maps to IRS1.",
    },
    "LIPIN1": {
        "mapping_type": "alias",
        "gene_symbols": ["LPIN1"],
        "notes": "Lipin 1 maps to LPIN1.",
    },
    "LIPIN 1": {
        "mapping_type": "alias",
        "gene_symbols": ["LPIN1"],
        "notes": "Lipin 1 maps to LPIN1.",
    },
    "ULK": {
        "mapping_type": "family",
        "gene_symbols": ["ULK1", "ULK2", "ULK3"],
        "notes": "ULK in signaling diagrams usually refers to the ULK kinase family.",
    },
    "RAGA/B": {
        "mapping_type": "family",
        "gene_symbols": ["RRAGA", "RRAGB"],
        "notes": "Rag A/B expands to RRAGA and RRAGB.",
    },
    "RAGC/D": {
        "mapping_type": "family",
        "gene_symbols": ["RRAGC", "RRAGD"],
        "notes": "Rag C/D expands to RRAGC and RRAGD.",
    },
    "LAMTOR": {
        "mapping_type": "family",
        "gene_symbols": ["LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5"],
        "notes": "LAMTOR expands to the Ragulator/LAMTOR complex components.",
    },
    "LAMTOR1/2/3/4/5": {
        "mapping_type": "family",
        "gene_symbols": ["LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5"],
        "notes": "LAMTOR1/2/3/4/5 expands to the Ragulator/LAMTOR complex components.",
    },
    "REDD1/2": {
        "mapping_type": "family",
        "gene_symbols": ["DDIT4", "DDIT4L"],
        "notes": "REDD1/2 expands to DDIT4 and DDIT4L.",
    },
    "FNIP1/2": {
        "mapping_type": "family",
        "gene_symbols": ["FNIP1", "FNIP2"],
        "notes": "FNIP1/2 expands to FNIP1 and FNIP2.",
    },
    "SESTRIN-1/2": {
        "mapping_type": "family",
        "gene_symbols": ["SESN1", "SESN2"],
        "notes": "Sestrin-1/2 expands to SESN1 and SESN2.",
    },
    "SESTRIN1/2": {
        "mapping_type": "family",
        "gene_symbols": ["SESN1", "SESN2"],
        "notes": "Sestrin-1/2 expands to SESN1 and SESN2.",
    },
    "SREBP-1": {
        "mapping_type": "alias",
        "gene_symbols": ["SREBF1"],
        "notes": "SREBP-1 maps to SREBF1.",
    },
    "GSK-3": {
        "mapping_type": "family",
        "gene_symbols": ["GSK3A", "GSK3B"],
        "notes": "GSK-3 expands to GSK3A and GSK3B.",
    },
    "4E-BP1/2": {
        "mapping_type": "family",
        "gene_symbols": ["EIF4EBP1", "EIF4EBP2"],
        "notes": "4E-BP1/2 expands to EIF4EBP1 and EIF4EBP2.",
    },
    "HIF-1": {
        "mapping_type": "family",
        "gene_symbols": ["HIF1A", "ARNT"],
        "notes": "HIF-1 is a transcription-factor complex; expanded to HIF1A and ARNT.",
    },
    "PGC-1": {
        "mapping_type": "family",
        "gene_symbols": ["PPARGC1A", "PPARGC1B"],
        "notes": "PGC-1 can refer to the PGC-1 coactivator family; expanded conservatively to PPARGC1A and PPARGC1B.",
    },
    "PGC1": {
        "mapping_type": "family",
        "gene_symbols": ["PPARGC1A", "PPARGC1B"],
        "notes": "PGC-1 can refer to the PGC-1 coactivator family; expanded conservatively to PPARGC1A and PPARGC1B.",
    },
    "GALPHAQ/O": {
        "mapping_type": "family",
        "gene_symbols": ["GNAQ", "GNAO1"],
        "notes": "Galphaq/o expands conservatively to GNAQ and GNAO1.",
    },
    "GBETAL": {
        "mapping_type": "alias",
        "gene_symbols": ["GNB2L1"],
        "notes": "GbetaL maps to GNB2L1 (RACK1).",
    },
}


def _is_reviewed_accession(accession: str) -> bool:
    return bool(accession) and not accession.startswith("A0A")


def _split_gene_names(raw: str) -> List[str]:
    if not raw:
        return []
    return [part.strip() for part in re.split(r"\s+", raw.strip()) if part.strip()]


def _normalize_whitespace(text: str) -> str:
    text = re.sub(r"[\u2010\u2011\u2012\u2013\u2014\u2212]", "-", text)
    text = text.replace("\u03ba", "k").replace("\u039a", "K")
    text = re.sub(r"\s+", " ", text.strip())
    return text


def _canonicalize_label(text: str) -> str:
    value = _normalize_whitespace(text).upper()
    value = value.replace("BETA", "B").replace("ALPHA", "A")
    value = value.replace(" ", "")
    value = value.replace(".", "")
    value = value.replace("(", "").replace(")", "")
    value = value.replace("[", "").replace("]", "")
    value = value.replace("{", "").replace("}", "")
    value = value.replace("_", "")
    return value


@dataclass
class UniProtCandidate:
    uniprot_id: str
    gene_symbol: str
    gene_symbols: List[str]
    protein_name: Optional[str] = None
    reviewed: Optional[bool] = None
    source: str = "local"


@dataclass
class LabelMappingResult:
    original_label: str
    normalized_label: str
    mapping_type: str
    suggested_uniprot_ids: List[str]
    suggested_gene_symbols: List[str]
    notes: str
    confidence: str
    candidate_details: List[Dict[str, Any]] = field(default_factory=list)


class PathwayLabelMapper:
    def __init__(
        self,
        organism_id: int = DEFAULT_ORGANISM_ID,
        annotation_file: Optional[Path] = None,
        cache_dir: Optional[Path] = None,
        use_uniprot_rest: bool = True,
        prefer_reviewed: bool = True,
    ) -> None:
        self.organism_id = int(organism_id)
        self.annotation_file = annotation_file or LOCAL_UNIPROT_ANNOTATIONS
        self.cache_dir = cache_dir or CACHE_DIR
        self.use_uniprot_rest = bool(use_uniprot_rest)
        self.prefer_reviewed = bool(prefer_reviewed)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._symbol_index: Dict[str, List[UniProtCandidate]] = {}
        self._canonical_symbol_index: Dict[str, List[UniProtCandidate]] = {}
        self._loaded_local_index = False

    def normalize_label(self, label: str) -> str:
        return _normalize_whitespace(label).upper()

    def resolve_alias(self, label: str) -> Optional[Dict[str, Any]]:
        normalized = self.normalize_label(label)
        stripped = re.sub(r"[-/]+$", "", normalized).strip()
        tokens = [part for part in re.split(r"[\s(),]+", normalized) if part]
        descriptive_context = {
            "ACTIVATION",
            "SIGNALING",
            "PATHWAY",
            "CASCADE",
            "RESPONSE",
            "RESPONSES",
            "NETWORK",
            "MODULE",
            "TRANSCRIPTION",
            "EXPRESSION",
            "REGULATION",
        }
        contains_context_word = any(token in descriptive_context for token in tokens)
        variants = {
            normalized,
            _canonicalize_label(normalized),
            _canonicalize_label(normalized.replace(" MAPK", "")),
        }
        if stripped:
            variants.add(stripped)
            variants.add(_canonicalize_label(stripped))
        if not contains_context_word:
            for token in tokens:
                variants.add(token)
                variants.add(_canonicalize_label(token))
        canonical = _canonicalize_label(normalized)
        if canonical.startswith("RAR"):
            variants.add("RAR")
        if "P65" in canonical:
            variants.add("P65")
        if canonical.startswith("NFKB"):
            variants.add("NFKB")
        for key in variants:
            rule = SIGNALING_LABEL_RULES.get(key)
            if rule:
                return dict(rule)
        return None

    def _load_local_symbol_index(self) -> None:
        if self._loaded_local_index:
            return
        self._loaded_local_index = True
        if self.organism_id != DEFAULT_ORGANISM_ID:
            LOGGER.info("Skipping local annotation load for non-human organism_id=%s", self.organism_id)
            return
        if not self.annotation_file.exists():
            LOGGER.warning("Local UniProt annotation file not found: %s", self.annotation_file)
            return

        with self.annotation_file.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                accession = str(row.get("Entry") or "").strip()
                if not accession:
                    continue
                gene_name_parts = _split_gene_names(str(row.get("Gene Names") or ""))
                if not gene_name_parts:
                    continue
                primary_symbol = gene_name_parts[0].upper()
                candidate = UniProtCandidate(
                    uniprot_id=accession,
                    gene_symbol=primary_symbol,
                    gene_symbols=[part.upper() for part in gene_name_parts],
                    protein_name=None,
                    reviewed=_is_reviewed_accession(accession),
                    source="local_annotations",
                )
                entry_name = str(row.get("Entry Name") or "").strip().upper()
                index_symbols = set(candidate.gene_symbols)
                if entry_name:
                    index_symbols.add(entry_name)
                    if entry_name.endswith("_HUMAN"):
                        index_symbols.add(entry_name[:-6])
                for symbol in index_symbols:
                    clean_symbol = str(symbol or "").strip().upper()
                    if not clean_symbol:
                        continue
                    self._symbol_index.setdefault(clean_symbol, []).append(candidate)
                    canonical = _canonicalize_label(clean_symbol)
                    if canonical:
                        self._canonical_symbol_index.setdefault(canonical, []).append(candidate)

    def _cache_file_for_symbol(self, symbol: str) -> Path:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", symbol.upper())
        return self.cache_dir / f"uniprot_v2_{self.organism_id}_{safe}.json"

    def _query_uniprot_rest(self, symbol: str) -> List[UniProtCandidate]:
        if not self.use_uniprot_rest or requests is None:
            return []
        cache_file = self._cache_file_for_symbol(symbol)
        if cache_file.exists():
            try:
                cached = json.loads(cache_file.read_text(encoding="utf-8"))
                return [UniProtCandidate(**item) for item in cached if isinstance(item, dict)]
            except Exception:
                pass

        headers = {"User-Agent": "MapKinase-PathwayLabelMapper/1.0"}
        def _request(query: str) -> List[UniProtCandidate]:
            params = {
                "query": query,
                "fields": "accession,gene_primary,gene_names,protein_name,reviewed",
                "format": "tsv",
                "size": "50",
            }
            try:
                response = requests.get(UNIPROT_SEARCH_URL, params=params, headers=headers, timeout=20)
                response.raise_for_status()
            except Exception as exc:
                LOGGER.warning("UniProt lookup failed for symbol '%s': %s", symbol, exc)
                return []

            lines = [line for line in response.text.splitlines() if line.strip()]
            if len(lines) <= 1:
                return []

            output: List[UniProtCandidate] = []
            reader = csv.DictReader(lines, delimiter="\t")
            for row in reader:
                accession = str(row.get("Entry") or row.get("accession") or "").strip()
                if not accession:
                    continue
                primary = str(row.get("Gene Names (primary)") or row.get("gene_primary") or "").strip().upper()
                all_names = _split_gene_names(str(row.get("Gene Names") or row.get("gene_names") or ""))
                if primary and primary not in all_names:
                    all_names = [primary] + all_names
                protein_name = str(row.get("Protein names") or row.get("protein_name") or "").strip() or None
                reviewed_text = str(row.get("Reviewed") or row.get("reviewed") or "").strip().lower()
                reviewed = reviewed_text in {"reviewed", "true", "yes"}
                if not primary and all_names:
                    primary = all_names[0].upper()
                output.append(
                    UniProtCandidate(
                        uniprot_id=accession,
                        gene_symbol=primary,
                        gene_symbols=[part.upper() for part in all_names if part],
                        protein_name=protein_name,
                        reviewed=reviewed if reviewed_text else _is_reviewed_accession(accession),
                        source="uniprot_rest",
                    )
                )
            return output

        candidates = _request(f'(gene_exact:{symbol}) AND (organism_id:{self.organism_id})')
        if not candidates:
            candidates = _request(f'(gene:{symbol}) AND (organism_id:{self.organism_id})')

        try:
            cache_file.write_text(json.dumps([asdict(item) for item in candidates], indent=2), encoding="utf-8")
        except Exception:
            pass
        return candidates

    def _sort_candidates(self, candidates: Iterable[UniProtCandidate], requested_symbol: str) -> List[UniProtCandidate]:
        requested = requested_symbol.upper()

        def _key(item: UniProtCandidate) -> tuple:
            reviewed_rank = 0 if (item.reviewed or not self.prefer_reviewed) else 1
            exact_primary_rank = 0 if item.gene_symbol == requested else 1
            synonym_rank = 0 if requested in item.gene_symbols else 1
            accession_rank = item.uniprot_id
            return (reviewed_rank, exact_primary_rank, synonym_rank, accession_rank)

        return sorted(candidates, key=_key)

    def _best_candidates_for_symbol(self, symbol: str) -> List[UniProtCandidate]:
        self._load_local_symbol_index()
        normalized = self.normalize_label(symbol)
        local = list(self._symbol_index.get(normalized, []))
        if not local:
            local = list(self._canonical_symbol_index.get(_canonicalize_label(normalized), []))
        online = self._query_uniprot_rest(normalized)
        pool = online or local
        if not pool:
            return []
        exact_primary_pool = [item for item in pool if item.gene_symbol == normalized]
        exact_synonym_pool = [item for item in pool if normalized in item.gene_symbols]
        pool = exact_primary_pool or exact_synonym_pool or pool

        grouped: Dict[str, List[UniProtCandidate]] = {}
        for item in pool:
            grouped.setdefault(item.gene_symbol or normalized, []).append(item)

        winners: List[UniProtCandidate] = []
        for gene_symbol, items in grouped.items():
            exact_primary = [item for item in items if item.gene_symbol == normalized]
            exact_synonym = [item for item in items if normalized in item.gene_symbols]
            pool = exact_primary or exact_synonym or items
            if self.prefer_reviewed:
                reviewed_pool = [item for item in pool if item.reviewed]
                if reviewed_pool:
                    pool = reviewed_pool
            ranked = self._sort_candidates(pool, normalized)
            if ranked:
                winners.append(ranked[0])
        return self._sort_candidates(winners, normalized)

    def lookup_uniprot(self, symbol_or_symbols: str | Sequence[str], organism: Optional[int] = None) -> List[UniProtCandidate]:
        if organism is not None and int(organism) != self.organism_id:
            temp = PathwayLabelMapper(
                organism_id=int(organism),
                annotation_file=self.annotation_file,
                cache_dir=self.cache_dir,
                use_uniprot_rest=self.use_uniprot_rest,
                prefer_reviewed=self.prefer_reviewed,
            )
            return temp.lookup_uniprot(symbol_or_symbols)

        if isinstance(symbol_or_symbols, str):
            symbols = [symbol_or_symbols]
        else:
            symbols = list(symbol_or_symbols)

        merged: List[UniProtCandidate] = []
        seen_accessions: set[str] = set()
        for symbol in symbols:
            for item in self._best_candidates_for_symbol(symbol):
                if item.uniprot_id in seen_accessions:
                    continue
                seen_accessions.add(item.uniprot_id)
                merged.append(item)
        return merged

    def _confidence_for_mapping_type(self, mapping_type: str) -> str:
        return {
            "exact_gene": "high",
            "alias": "medium-high",
            "family": "medium",
            "ambiguous": "low",
            "unresolved": "none",
        }.get(mapping_type, "low")

    def _result_from_candidates(
        self,
        *,
        original_label: str,
        normalized_label: str,
        mapping_type: str,
        notes: str,
        candidates: Sequence[UniProtCandidate],
        fallback_gene_symbols: Optional[Sequence[str]] = None,
    ) -> LabelMappingResult:
        gene_symbols = [item.gene_symbol for item in candidates if item.gene_symbol]
        if not gene_symbols and fallback_gene_symbols:
            gene_symbols = [str(item).upper() for item in fallback_gene_symbols if str(item).strip()]
        unique_gene_symbols = list(dict.fromkeys(gene_symbols))
        unique_uniprots = list(dict.fromkeys(item.uniprot_id for item in candidates if item.uniprot_id))
        candidate_details = [asdict(item) for item in candidates]
        return LabelMappingResult(
            original_label=original_label,
            normalized_label=normalized_label,
            mapping_type=mapping_type,
            suggested_uniprot_ids=unique_uniprots,
            suggested_gene_symbols=unique_gene_symbols,
            notes=notes,
            confidence=self._confidence_for_mapping_type(mapping_type),
            candidate_details=candidate_details,
        )

    def map_pathway_label(self, label: str) -> Dict[str, Any]:
        original = str(label or "")
        normalized = self.normalize_label(original)
        if not normalized:
            return asdict(
                LabelMappingResult(
                    original_label=original,
                    normalized_label=normalized,
                    mapping_type="unresolved",
                    suggested_uniprot_ids=[],
                    suggested_gene_symbols=[],
                    notes="Label is empty after normalization.",
                    confidence="none",
                )
            )

        canonical = _canonicalize_label(normalized)
        if canonical in {"PIP", "PIP2", "PIP3", "ATP", "ADP", "AMP", "GTP", "GDP", "DNA", "MRNA"}:
            return asdict(
                LabelMappingResult(
                    original_label=original,
                    normalized_label=normalized,
                    mapping_type="unresolved",
                    suggested_uniprot_ids=[],
                    suggested_gene_symbols=[],
                    notes="Common metabolite or nucleic-acid label, not treated as a protein.",
                    confidence="none",
                )
            )
        if canonical in {"MTORC1", "MTORC2"}:
            return asdict(
                LabelMappingResult(
                    original_label=original,
                    normalized_label=normalized,
                    mapping_type="unresolved",
                    suggested_uniprot_ids=[],
                    suggested_gene_symbols=[],
                    notes="mTORC1 and mTORC2 are complexes, not single proteins.",
                    confidence="none",
                )
            )
        if not re.search(r"[A-Z]", canonical) or not re.search(r"[A-Z]", normalized):
            return asdict(
                LabelMappingResult(
                    original_label=original,
                    normalized_label=normalized,
                    mapping_type="unresolved",
                    suggested_uniprot_ids=[],
                    suggested_gene_symbols=[],
                    notes="Label does not look like a protein or gene symbol.",
                    confidence="none",
                )
            )

        rule = self.resolve_alias(normalized)
        if rule:
            gene_symbols = [str(item).upper() for item in rule.get("gene_symbols", [])]
            candidates = self.lookup_uniprot(gene_symbols)
            return asdict(
                self._result_from_candidates(
                    original_label=original,
                    normalized_label=normalized,
                    mapping_type=str(rule.get("mapping_type") or "alias"),
                    notes=str(rule.get("notes") or "Matched using curated pathway alias/family rules."),
                    candidates=candidates,
                    fallback_gene_symbols=gene_symbols,
                )
            )

        candidates = self.lookup_uniprot(normalized)
        if candidates:
            exact_primary = [item for item in candidates if item.gene_symbol == normalized]
            exact_alias = [item for item in candidates if normalized in item.gene_symbols]
            gene_symbols = list(dict.fromkeys(item.gene_symbol for item in candidates if item.gene_symbol))
            if exact_primary and len(set(item.gene_symbol for item in exact_primary)) == 1:
                return asdict(
                    self._result_from_candidates(
                        original_label=original,
                        normalized_label=normalized,
                        mapping_type="exact_gene",
                        notes=f"Matched directly to the human gene symbol {exact_primary[0].gene_symbol}.",
                        candidates=exact_primary,
                    )
                )
            if exact_alias and len(set(item.gene_symbol for item in exact_alias)) == 1:
                return asdict(
                    self._result_from_candidates(
                        original_label=original,
                        normalized_label=normalized,
                        mapping_type="alias",
                        notes=f"Matched via alias/synonym to human gene symbol {exact_alias[0].gene_symbol}.",
                        candidates=exact_alias,
                    )
                )
            return asdict(
                self._result_from_candidates(
                    original_label=original,
                    normalized_label=normalized,
                    mapping_type="ambiguous",
                    notes=f"Multiple human gene symbols matched this label: {', '.join(gene_symbols)}.",
                    candidates=candidates,
                )
            )

        return asdict(
            LabelMappingResult(
                original_label=original,
                normalized_label=normalized,
                mapping_type="unresolved",
                suggested_uniprot_ids=[],
                suggested_gene_symbols=[],
                notes="No curated alias or UniProt-backed symbol match was found.",
                confidence="none",
            )
        )

    def map_pathway_labels(self, label_list: Sequence[str]) -> List[Dict[str, Any]]:
        return [self.map_pathway_label(label) for label in label_list]

    def build_summary(self, results: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
        summary = {
            "total_labels_processed": len(results),
            "exact_matches": 0,
            "alias_matches": 0,
            "family_expansions": 0,
            "ambiguous_labels": 0,
            "unresolved_labels": 0,
        }
        for item in results:
            mapping_type = str(item.get("mapping_type") or "").strip().lower()
            if mapping_type == "exact_gene":
                summary["exact_matches"] += 1
            elif mapping_type == "alias":
                summary["alias_matches"] += 1
            elif mapping_type == "family":
                summary["family_expansions"] += 1
            elif mapping_type == "ambiguous":
                summary["ambiguous_labels"] += 1
            elif mapping_type == "unresolved":
                summary["unresolved_labels"] += 1
        return summary

    def export_results_to_csv(self, results: Sequence[Dict[str, Any]], output_path: str | Path) -> Path:
        path = Path(output_path)
        fieldnames = [
            "original_label",
            "normalized_label",
            "mapping_type",
            "suggested_gene_symbols",
            "suggested_uniprot_ids",
            "confidence",
            "notes",
        ]
        with path.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            for item in results:
                row = dict(item)
                row["suggested_gene_symbols"] = ";".join(str(x) for x in row.get("suggested_gene_symbols", []) or [])
                row["suggested_uniprot_ids"] = ";".join(str(x) for x in row.get("suggested_uniprot_ids", []) or [])
                writer.writerow({key: row.get(key, "") for key in fieldnames})
        return path

    def export_results_to_json(self, results: Sequence[Dict[str, Any]], output_path: str | Path) -> Path:
        path = Path(output_path)
        payload = {
            "summary": self.build_summary(results),
            "results": list(results),
        }
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        return path


def normalize_label(label: str) -> str:
    return PathwayLabelMapper().normalize_label(label)


def resolve_alias(label: str) -> Optional[Dict[str, Any]]:
    return PathwayLabelMapper().resolve_alias(label)


def lookup_uniprot(symbol_or_symbols: str | Sequence[str], organism: int = DEFAULT_ORGANISM_ID) -> List[Dict[str, Any]]:
    mapper = PathwayLabelMapper(organism_id=organism)
    return [asdict(item) for item in mapper.lookup_uniprot(symbol_or_symbols)]


def map_pathway_label(label: str, organism: int = DEFAULT_ORGANISM_ID) -> Dict[str, Any]:
    mapper = PathwayLabelMapper(organism_id=organism)
    return mapper.map_pathway_label(label)


def map_pathway_labels(label_list: Sequence[str], organism: int = DEFAULT_ORGANISM_ID) -> List[Dict[str, Any]]:
    mapper = PathwayLabelMapper(organism_id=organism)
    return mapper.map_pathway_labels(label_list)


def export_results_to_csv(results: Sequence[Dict[str, Any]], output_path: str | Path) -> Path:
    mapper = PathwayLabelMapper()
    return mapper.export_results_to_csv(results, output_path)


def export_results_to_json(results: Sequence[Dict[str, Any]], output_path: str | Path) -> Path:
    mapper = PathwayLabelMapper()
    return mapper.export_results_to_json(results, output_path)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    example_labels = [
        "MAPK14",
        "p38",
        "ERK",
        "ERK1/2",
        "MEK1",
        "JNK",
        "AKT",
        "RAS",
        "PI3K",
        "c-Jun",
        "Elk-1",
        "PRAK",
        "MAPKAPK-3",
        "MKK3/6",
    ]
    mapper = PathwayLabelMapper()
    results = mapper.map_pathway_labels(example_labels)
    print(json.dumps({"summary": mapper.build_summary(results), "results": results}, indent=2))
