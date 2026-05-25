import csv
import io
import logging
import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from MapKinase_WebApp.mk_security_limits import MAX_UPLOAD_SIZE_BYTES, MAX_UPLOAD_SIZE_MB

# Analysis dataset uploads are intentionally narrow to reduce attack surface.
# Custom pathway JSON imports are handled in a dedicated, separate handler.
ALLOWED_EXTENSIONS = {".csv", ".tsv", ".txt"}
CUSTOM_PATHWAY_ALLOWED_EXTENSIONS = {".json"}
_MAX_VALIDATION_BYTES = 128 * 1024
_MAX_VALIDATION_LINES = 120
_FILENAME_CONTROL_CHARS = re.compile(r"[\x00-\x1f\x7f]")
_SAFE_TEXT_BYTES = set(range(32, 127)) | {9, 10, 13}
_UPLOAD_VALIDATION_LOGGER = logging.getLogger("mapkinase.upload_validation")

if not _UPLOAD_VALIDATION_LOGGER.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(name)s - %(message)s"))
    _UPLOAD_VALIDATION_LOGGER.addHandler(_handler)
_UPLOAD_VALIDATION_LOGGER.setLevel(logging.INFO)


@dataclass
class UploadValidationResult:
    valid: bool
    user_message: str
    sanitized_filename: str
    extension: str
    datapath: str


def _log_upload_rejection(
    filename: str,
    reason: str,
    expected_role: Optional[str],
    file_size_bytes: Optional[int] = None,
    session_id: Optional[str] = None,
) -> None:
    role = (expected_role or "unspecified").strip().lower() or "unspecified"
    safe_name = filename or "<unknown>"
    _UPLOAD_VALIDATION_LOGGER.warning(
        "Rejected upload | session=%s | role=%s | filename=%s | size_bytes=%s | reason=%s",
        (session_id or "unknown"),
        role,
        safe_name,
        str(file_size_bytes) if file_size_bytes is not None else "na",
        reason,
    )


def _header_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value or "").strip().lower())


def _is_text_like(payload: bytes) -> bool:
    if not payload:
        return False
    if b"\x00" in payload:
        return False
    sample = payload[:4096]
    if not sample:
        return False
    non_text = 0
    for b in sample:
        if b not in _SAFE_TEXT_BYTES:
            non_text += 1
    return (non_text / max(len(sample), 1)) <= 0.20


def _extract_upload_fields(file_info: Any) -> Tuple[str, str]:
    if isinstance(file_info, dict):
        raw_name = str(file_info.get("name") or "")
        datapath = str(file_info.get("datapath") or "")
        return raw_name, datapath
    raw_name = str(getattr(file_info, "name", "") or "")
    datapath = str(getattr(file_info, "datapath", "") or "")
    return raw_name, datapath


def _validate_role_headers(headers: Sequence[str], expected_role: Optional[str]) -> Optional[str]:
    role = str(expected_role or "").strip().lower()
    if not headers:
        return "Invalid analysis dataset upload: file could not be parsed as tabular text."

    h1 = _header_key(headers[0]) if len(headers) >= 1 else ""
    h2 = _header_key(headers[1]) if len(headers) >= 2 else ""

    uniprot_aliases = {"uniprot", "uniprotid", "uniprotkb", "uniprotaccession", "uniprotkbaccession", "accession"}
    gene_aliases = {"gene", "genesymbol", "symbol", "genesym"}
    ptm_pos_aliases = {"siteposition", "ptmsite", "modifiedresidueposition", "residueposition", "position"}
    hmdb_aliases = {"hmdbid", "hmdb"}

    if role == "protein":
        if len(headers) < 2:
            return "Invalid analysis dataset upload: protein dataset must include at least two columns."
        if h1 not in uniprot_aliases:
            return "Invalid analysis dataset upload: protein dataset first column must be UniProtKB accession."
        if h2 not in gene_aliases:
            return "Invalid analysis dataset upload: protein dataset second column must be gene symbol."
    elif role == "ptm":
        if len(headers) < 2:
            return "Invalid analysis dataset upload: PTM dataset must include at least two columns."
        if h1 not in uniprot_aliases:
            return "Invalid analysis dataset upload: PTM dataset first column must be UniProtKB accession."
        if h2 not in ptm_pos_aliases:
            return "Invalid analysis dataset upload: PTM dataset second column must be modified residue position."
    elif role == "metabolite":
        if len(headers) < 1:
            return "Invalid analysis dataset upload: metabolite dataset must include at least one column."
        if h1 not in hmdb_aliases:
            return "Invalid analysis dataset upload: metabolite dataset first column must be HMDB ID."
    return None


def validate_uploaded_file(file_info: Any, expected_role: Optional[str] = None, session_id: Optional[str] = None) -> UploadValidationResult:
    """
    Centralized upload gate.

    Roles:
      - protein / ptm / metabolite: analysis dataset uploads (.csv/.tsv/.txt only)
      - custom_pathway: separate pathway import flow (.json only)

    Enforcement order for analysis dataset uploads:
      A) filename/path sanitization
      B) extension allowlist
      C) per-file size limit
      D) text-vs-binary content checks
      E) tabular parse checks
      F) role-specific required columns
    """
    raw_name, datapath = _extract_upload_fields(file_info)
    role = str(expected_role or "").strip().lower()
    filename = raw_name.strip()
    if not filename and datapath:
        filename = Path(datapath).name
    sanitized_filename = Path(filename).name

    def _reject(user_message: str, reason: str, file_size_bytes: Optional[int] = None) -> UploadValidationResult:
        if role == "custom_pathway" and user_message.startswith("Invalid upload:"):
            user_message = user_message.replace("Invalid upload:", "Invalid custom pathway import:", 1)
        _log_upload_rejection(sanitized_filename, reason, expected_role, file_size_bytes=file_size_bytes, session_id=session_id)
        return UploadValidationResult(
            valid=False,
            user_message=user_message,
            sanitized_filename=sanitized_filename,
            extension=Path(sanitized_filename).suffix.lower(),
            datapath=datapath,
        )

    if not filename:
        return _reject("Invalid upload: missing filename.", "missing filename")
    if not datapath or not os.path.isfile(datapath):
        return _reject("Invalid upload: file could not be read.", "missing uploaded temp file")
    if filename != sanitized_filename:
        return _reject("Invalid upload: filename is unsafe.", "filename includes path elements")
    if any(token in filename for token in ("../", "..\\", "/", "\\")) or ".." in filename:
        return _reject("Invalid upload: filename is unsafe.", "path traversal token in filename")
    if _FILENAME_CONTROL_CHARS.search(filename):
        return _reject("Invalid upload: filename contains unsupported characters.", "control character in filename")
    if "\x00" in filename:
        return _reject("Invalid upload: filename contains unsupported characters.", "null byte in filename")
    if sanitized_filename in {"", ".", ".."}:
        return _reject("Invalid upload: filename is unsafe.", "empty or reserved filename")
    if sanitized_filename.startswith("."):
        return _reject("Invalid upload: hidden filenames are not allowed.", "hidden filename")

    extension = Path(sanitized_filename).suffix.lower()
    allowed_extensions = CUSTOM_PATHWAY_ALLOWED_EXTENSIONS if role == "custom_pathway" else ALLOWED_EXTENSIONS
    if extension not in allowed_extensions:
        if role == "custom_pathway":
            return _reject("Invalid custom pathway import: only .json files are accepted.", f"extension '{extension}' not allowed for custom pathway")
        return _reject("Invalid analysis dataset upload: only .csv, .tsv, and tab-delimited .txt files are accepted.", f"extension '{extension}' not allowlisted")

    try:
        file_size_bytes = Path(datapath).stat().st_size
    except OSError:
        return _reject("Invalid upload: file could not be read.", "failed to stat uploaded file")

    if file_size_bytes > MAX_UPLOAD_SIZE_BYTES:
        if role == "custom_pathway":
            return _reject(
                f"Upload rejected: custom pathway import files must be {MAX_UPLOAD_SIZE_MB} MB or smaller.",
                f"file size {file_size_bytes} exceeds max {MAX_UPLOAD_SIZE_BYTES}",
                file_size_bytes=file_size_bytes,
            )
        return _reject(
            f"Upload rejected: each file must be {MAX_UPLOAD_SIZE_MB} MB or smaller.",
            f"file size {file_size_bytes} exceeds max {MAX_UPLOAD_SIZE_BYTES}",
            file_size_bytes=file_size_bytes,
        )

    try:
        with open(datapath, "rb") as fh:
            sample_bytes = fh.read(_MAX_VALIDATION_BYTES)
    except OSError:
        return _reject("Invalid upload: file could not be read.", "failed to read uploaded bytes", file_size_bytes=file_size_bytes)

    if not sample_bytes:
        if role == "custom_pathway":
            return _reject("Invalid custom pathway import: file is empty.", "empty file", file_size_bytes=file_size_bytes)
        return _reject("Invalid analysis dataset upload: file is empty.", "empty file", file_size_bytes=file_size_bytes)
    if b"\x00" in sample_bytes:
        if role == "custom_pathway":
            return _reject("Invalid custom pathway import: file appears binary and is not accepted.", "null byte found in file content", file_size_bytes=file_size_bytes)
        return _reject("Invalid analysis dataset upload: file appears binary and is not accepted.", "null byte found in file content", file_size_bytes=file_size_bytes)
    if not _is_text_like(sample_bytes):
        if role == "custom_pathway":
            return _reject("Invalid custom pathway import: file appears binary and is not accepted.", "byte distribution appears binary", file_size_bytes=file_size_bytes)
        return _reject("Invalid analysis dataset upload: file appears binary and is not accepted.", "byte distribution appears binary", file_size_bytes=file_size_bytes)

    decoded = sample_bytes.decode("utf-8", errors="replace")
    if not decoded.strip():
        if role == "custom_pathway":
            return _reject("Invalid custom pathway import: file could not be decoded as text.", "decoded text is blank", file_size_bytes=file_size_bytes)
        return _reject("Invalid analysis dataset upload: file could not be parsed as tabular text.", "decoded text is blank", file_size_bytes=file_size_bytes)
    replacement_ratio = decoded.count("\ufffd") / max(len(decoded), 1)
    if replacement_ratio > 0.10:
        if role == "custom_pathway":
            return _reject("Invalid custom pathway import: file appears binary and is not accepted.", "decoded text contains too many replacement characters", file_size_bytes=file_size_bytes)
        return _reject("Invalid analysis dataset upload: file appears binary and is not accepted.", "decoded text contains too many replacement characters", file_size_bytes=file_size_bytes)

    if role == "custom_pathway":
        return UploadValidationResult(True, "", sanitized_filename, extension, datapath)

    delimiter = "," if extension == ".csv" else "\t"
    sample_lines = decoded.splitlines()[:_MAX_VALIDATION_LINES]
    sample_text = "\n".join(sample_lines)
    if delimiter not in sample_text:
        return _reject("Invalid analysis dataset upload: file could not be parsed as tabular text.", f"delimiter '{delimiter}' not found in validation sample", file_size_bytes=file_size_bytes)
    try:
        csv.Sniffer().sniff(sample_text, delimiters=[delimiter])
    except Exception:
        return _reject("Invalid analysis dataset upload: file could not be parsed as tabular text.", "csv sniffer could not identify required delimiter", file_size_bytes=file_size_bytes)

    parsed_rows: List[List[str]] = []
    try:
        reader = csv.reader(io.StringIO(sample_text), delimiter=delimiter)
        for row in reader:
            if not row or not any(str(cell or "").strip() for cell in row):
                continue
            parsed_rows.append([str(cell or "").strip() for cell in row])
    except Exception:
        return _reject("Invalid analysis dataset upload: file could not be parsed as tabular text.", "csv parser failed on validation sample", file_size_bytes=file_size_bytes)

    if len(parsed_rows) < 2:
        return _reject("Invalid analysis dataset upload: file must include a header row and at least one data row.", "missing header or data row", file_size_bytes=file_size_bytes)
    headers = parsed_rows[0]
    if not headers or all(not str(h).strip() for h in headers):
        return _reject("Invalid analysis dataset upload: file has no valid column headers.", "empty header row", file_size_bytes=file_size_bytes)
    if len(headers) == 0:
        return _reject("Invalid analysis dataset upload: file has no columns.", "zero column header", file_size_bytes=file_size_bytes)

    valid_data_rows = 0
    malformed_rows = 0
    for row in parsed_rows[1:]:
        if len(row) != len(headers):
            malformed_rows += 1
            continue
        if any(str(cell or "").strip() for cell in row):
            valid_data_rows += 1
    if valid_data_rows == 0:
        return _reject("Invalid analysis dataset upload: file could not be parsed as tabular text.", "no valid data rows with header-aligned columns", file_size_bytes=file_size_bytes)
    if malformed_rows > 0:
        return _reject("Invalid analysis dataset upload: file has malformed rows and is not accepted.", f"{malformed_rows} malformed rows in validation sample", file_size_bytes=file_size_bytes)

    role_header_error = _validate_role_headers(headers, expected_role)
    if role_header_error:
        return _reject(role_header_error, "role-specific required columns not present", file_size_bytes=file_size_bytes)

    return UploadValidationResult(
        valid=True,
        user_message="",
        sanitized_filename=sanitized_filename,
        extension=extension,
        datapath=datapath,
    )


def _clean_header(header: str, idx: int) -> str:
    h = header or ""
    if idx == 0:
        h = h.lstrip("\ufeff")
    return h.strip()

def _normalize_fc_suffix(header: str) -> str:
    cleaned = (header or "").strip()
    if ":" in cleaned:
        cleaned = cleaned.split(":", 1)[1]
    cleaned = re.sub(r"\s+", " ", cleaned.strip())
    return cleaned.lower()


class ProteinValidationResult:
    def __init__(self, valid: bool, errors: List[str], summary: Dict[str, int], comparisons: List[str]):
        self.valid = valid
        self.errors = errors
        self.summary = summary
        self.comparisons = comparisons


class PTMValidationResult:
    def __init__(self, valid: bool, errors: List[str], summary: Dict[str, int], comparisons: List[str]):
        self.valid = valid
        self.errors = errors
        self.summary = summary
        self.comparisons = comparisons


def _detect_delimiter(file_path: str) -> str:
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".csv":
        return ","
    return "\t"


def validate_protein_file(file_path: str) -> ProteinValidationResult:
    """
    Validate a protein data file according to the expected schema.
    Rules:
      - Col0: Uniprot ID (required values)
      - Col1: Gene Symbol (blank values are auto-filled from Uniprot ID)
      - Col2+: Comparison columns must start with "C:" (at least one required)
      - Optional outline comparison columns start with "O:" and must match a "C:" header
      - Optional tooltip columns start with "T:"
      - Comparison cells must be numeric (float/int) and non-empty
      - Outline comparison cells must be numeric (float/int) or "NA"
      - Uniprot/Comparison cells required on every row
      - Only .txt (tab-delimited) or .csv files are allowed
    Returns a ProteinValidationResult with validity, errors, and summary counts.
    """
    if not os.path.isfile(file_path):
        return ProteinValidationResult(False, [f"File not found: {file_path}"], {}, [])
    ext = os.path.splitext(file_path)[1].lower()
    if ext not in {".txt", ".tsv", ".csv"}:
        return ProteinValidationResult(False, [f"Unsupported file type '{ext}'. Use .txt (tab) or .csv."], {}, [])

    delimiter = _detect_delimiter(file_path)
    errors: List[str] = []
    comparison_col_indexes: List[int] = []
    outline_col_indexes: List[int] = []
    tooltip_col_indexes: List[int] = []
    row_count = 0

    try:
        with open(file_path, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter=delimiter)
            try:
                headers = next(reader)
            except StopIteration:
                return ProteinValidationResult(False, ["The file is empty."], {}, [])
            row_count = 0
            if len(headers) < 3:
                errors.append("Header must contain at least 3 columns (Uniprot, Gene Symbol, and one Comparison column).")
                return ProteinValidationResult(False, errors, {"rows": row_count, "comparisons": 0, "tooltips": 0}, [])

            for idx, header in enumerate(headers):
                header_clean = _clean_header(header, idx)
                if idx == 0 or idx == 1:
                    continue
                if header_clean.lower().startswith("c:"):
                    comparison_col_indexes.append(idx)
                elif header_clean.lower().startswith("o:"):
                    outline_col_indexes.append(idx)
                elif header_clean.lower().startswith("t:"):
                    tooltip_col_indexes.append(idx)
                else:
                    errors.append(
                        f"Invalid header in column {idx + 1}: '{header_clean}'. "
                        "Comparison headers must start with 'C:', outline headers with 'O:', and tooltip headers with 'T:'."
                    )

            if not comparison_col_indexes:
                errors.append("At least one Comparison column (header starting with 'C:') is required (third column or later).")
            comparison_headers = [_clean_header(headers[idx], idx) for idx in comparison_col_indexes]
            comparison_suffixes = {_normalize_fc_suffix(h) for h in comparison_headers if _normalize_fc_suffix(h)}
            for idx in outline_col_indexes:
                outline_header = _clean_header(headers[idx], idx)
                outline_suffix = _normalize_fc_suffix(outline_header)
                if outline_suffix and outline_suffix not in comparison_suffixes:
                    errors.append(
                        f"Outline header '{outline_header}' must match a Comparison header (C:) with the same label."
                    )

            for row_idx, row in enumerate(reader, start=2):  # 1-based with header at row 1
                row_count += 1
                # Normalize row length
                if len(row) < len(headers):
                    row.extend([""] * (len(headers) - len(row)))
                uniprot = (row[0] or "").strip()
                gene_symbol = (row[1] or "").strip()
                if not uniprot:
                    errors.append(f"Row {row_idx}, Column 1 (Uniprot ID) is empty.")
                if not gene_symbol and uniprot:
                    # Allow blank gene symbols by falling back to UniProt ID.
                    row[1] = uniprot
                for c_idx in comparison_col_indexes:
                    cell = (row[c_idx] or "").strip()
                    if not cell:
                        errors.append(f"Row {row_idx}, Column {c_idx + 1} (Comparison) is empty.")
                        continue
                    try:
                        float(cell)
                    except ValueError:
                        errors.append(f"Row {row_idx}, Column {c_idx + 1} (Comparison) must be numeric. Found '{cell}'.")
                for o_idx in outline_col_indexes:
                    cell = (row[o_idx] or "").strip()
                    if not cell or cell.strip().lower() == "na":
                        continue
                    try:
                        float(cell)
                    except ValueError:
                        errors.append(f"Row {row_idx}, Column {o_idx + 1} (Outline) must be numeric or NA. Found '{cell}'.")
                # Tooltip columns: values optional, no validation needed

                if len(errors) >= 200:
                    errors.append("Stopped after 200 errors; fix these issues and retry.")
                    break

    except UnicodeDecodeError:
        errors.append("File could not be decoded as UTF-8. Please provide UTF-8 encoded text.")
    except Exception as exc:
        errors.append(f"Unexpected error while reading file: {exc}")

    valid = len(errors) == 0
    summary = {
        "rows": row_count,
        "comparisons": len(comparison_col_indexes),
        "tooltips": len(tooltip_col_indexes),
    }
    comparison_headers = [_clean_header(headers[idx], idx) for idx in comparison_col_indexes]
    return ProteinValidationResult(valid, errors, summary, comparison_headers)


def validate_ptm_file(file_path: str, required_comparisons: List[str]) -> PTMValidationResult:
    """
    Validate a PTM data file.
    Rules:
      - Col0: Uniprot ID (required values)
      - Col1: Site Position (required, positive int)
      - Col2+: Comparison columns must start with "C:" (at least one required)
      - Optional outline comparison columns start with "O:" and must match a "C:" header
      - Optional tooltip columns start with "T:"
      - Comparison cells must be numeric (float/int) and non-empty
      - Outline comparison cells must be numeric (float/int) or "NA"
      - Site position must be an integer > 0
      - Uniprot/Site/Comparison cells required on every row
      - File must be .txt (tab) or .csv
      - All required_comparisons (from protein file) must exist in PTM headers (case-sensitive exact match)
    """
    if not os.path.isfile(file_path):
        return PTMValidationResult(False, [f"File not found: {file_path}"], {}, [])
    ext = os.path.splitext(file_path)[1].lower()
    if ext not in {".txt", ".tsv", ".csv"}:
        return PTMValidationResult(False, [f"Unsupported file type '{ext}'. Use .txt (tab) or .csv."], {}, [])

    delimiter = _detect_delimiter(file_path)
    errors: List[str] = []
    comparison_col_indexes: List[int] = []
    outline_col_indexes: List[int] = []
    tooltip_col_indexes: List[int] = []
    row_count = 0

    try:
        with open(file_path, newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter=delimiter)
            try:
                headers = next(reader)
            except StopIteration:
                return PTMValidationResult(False, ["The file is empty."], {}, [])
            row_count = 0
            if len(headers) < 3:
                errors.append("Header must contain at least 3 columns (Uniprot, Site Position, and one Comparison column).")
                return PTMValidationResult(False, errors, {"rows": row_count, "comparisons": 0, "tooltips": 0}, [])

            for idx, header in enumerate(headers):
                header_clean = _clean_header(header, idx)
                if idx == 0 or idx == 1:
                    continue
                if header_clean.lower().startswith("c:"):
                    comparison_col_indexes.append(idx)
                elif header_clean.lower().startswith("o:"):
                    outline_col_indexes.append(idx)
                elif header_clean.lower().startswith("t:"):
                    tooltip_col_indexes.append(idx)
                else:
                    errors.append(
                        f"Invalid header in column {idx + 1}: '{header_clean}'. "
                        "Comparison headers must start with 'C:', outline headers with 'O:', and tooltip headers with 'T:'."
                    )

            if not comparison_col_indexes:
                errors.append("At least one Comparison column (header starting with 'C:') is required (third column or later).")

            missing = [c for c in required_comparisons if c not in [_clean_header(headers[idx], idx) for idx in comparison_col_indexes]]
            if missing:
                errors.append(f"Missing Comparison columns present in protein file: {', '.join(missing)}.")
            comparison_headers = [_clean_header(headers[idx], idx) for idx in comparison_col_indexes]
            comparison_suffixes = {_normalize_fc_suffix(h) for h in comparison_headers if _normalize_fc_suffix(h)}
            for idx in outline_col_indexes:
                outline_header = _clean_header(headers[idx], idx)
                outline_suffix = _normalize_fc_suffix(outline_header)
                if outline_suffix and outline_suffix not in comparison_suffixes:
                    errors.append(
                        f"Outline header '{outline_header}' must match a Comparison header (C:) with the same label."
                    )

            for row_idx, row in enumerate(reader, start=2):
                row_count += 1
                if len(row) < len(headers):
                    row.extend([""] * (len(headers) - len(row)))
                uniprot = (row[0] or "").strip()
                site = (row[1] or "").strip()
                if not uniprot:
                    errors.append(f"Row {row_idx}, Column 1 (Uniprot ID) is empty.")
                if not site:
                    errors.append(f"Row {row_idx}, Column 2 (Site Position) is empty.")
                else:
                    try:
                        site_val = int(site)
                        if site_val <= 0:
                            errors.append(f"Row {row_idx}, Column 2 (Site Position) must be a positive integer. Found '{site}'.")
                    except ValueError:
                        errors.append(f"Row {row_idx}, Column 2 (Site Position) must be an integer. Found '{site}'.")
                for c_idx in comparison_col_indexes:
                    cell = (row[c_idx] or "").strip()
                    if not cell:
                        errors.append(f"Row {row_idx}, Column {c_idx + 1} (Comparison) is empty.")
                        continue
                    try:
                        float(cell)
                    except ValueError:
                        errors.append(f"Row {row_idx}, Column {c_idx + 1} (Comparison) must be numeric. Found '{cell}'.")
                for o_idx in outline_col_indexes:
                    cell = (row[o_idx] or "").strip()
                    if not cell or cell.strip().lower() == "na":
                        continue
                    try:
                        float(cell)
                    except ValueError:
                        errors.append(f"Row {row_idx}, Column {o_idx + 1} (Outline) must be numeric or NA. Found '{cell}'.")
                if len(errors) >= 200:
                    errors.append("Stopped after 200 errors; fix these issues and retry.")
                    break

    except UnicodeDecodeError:
        errors.append("File could not be decoded as UTF-8. Please provide UTF-8 encoded text.")
    except Exception as exc:
        errors.append(f"Unexpected error while reading file: {exc}")

    valid = len(errors) == 0
    summary = {
        "rows": row_count,
        "comparisons": len(comparison_col_indexes),
        "tooltips": len(tooltip_col_indexes),
    }
    comparison_headers = [_clean_header(headers[idx], idx) for idx in comparison_col_indexes] if comparison_col_indexes else []
    return PTMValidationResult(valid, errors, summary, comparison_headers)
