import csv
import os
from typing import Dict, List, Tuple


def _clean_header(header: str, idx: int) -> str:
    h = header or ""
    if idx == 0:
        h = h.lstrip("\ufeff")
    return h.strip()


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
      - Col1: Gene Symbol (required values)
      - Col2+: Comparison columns must start with "C:" (at least one required)
      - Optional tooltip columns start with "T:"
      - Comparison cells must be numeric (float/int) and non-empty
      - Uniprot/GeneSymbol/Comparison cells required on every row
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
                elif header_clean.lower().startswith("t:"):
                    tooltip_col_indexes.append(idx)
                else:
                    errors.append(f"Invalid header in column {idx + 1}: '{header_clean}'. Comparison headers must start with 'C:' and tooltip headers with 'T:'.")

            if not comparison_col_indexes:
                errors.append("At least one Comparison column (header starting with 'C:') is required (third column or later).")

            for row_idx, row in enumerate(reader, start=2):  # 1-based with header at row 1
                row_count += 1
                # Normalize row length
                if len(row) < len(headers):
                    row.extend([""] * (len(headers) - len(row)))
                uniprot = (row[0] or "").strip()
                gene_symbol = (row[1] or "").strip()
                if not uniprot:
                    errors.append(f"Row {row_idx}, Column 1 (Uniprot ID) is empty.")
                if not gene_symbol:
                    errors.append(f"Row {row_idx}, Column 2 (Gene Symbol) is empty.")
                for c_idx in comparison_col_indexes:
                    cell = (row[c_idx] or "").strip()
                    if not cell:
                        errors.append(f"Row {row_idx}, Column {c_idx + 1} (Comparison) is empty.")
                        continue
                    try:
                        float(cell)
                    except ValueError:
                        errors.append(f"Row {row_idx}, Column {c_idx + 1} (Comparison) must be numeric. Found '{cell}'.")
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
      - Optional tooltip columns start with "T:"
      - Comparison cells must be numeric (float/int) and non-empty
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
                elif header_clean.lower().startswith("t:"):
                    tooltip_col_indexes.append(idx)
                else:
                    errors.append(f"Invalid header in column {idx + 1}: '{header_clean}'. Comparison headers must start with 'C:' and tooltip headers with 'T:'.")

            if not comparison_col_indexes:
                errors.append("At least one Comparison column (header starting with 'C:') is required (third column or later).")

            missing = [c for c in required_comparisons if c not in [_clean_header(headers[idx], idx) for idx in comparison_col_indexes]]
            if missing:
                errors.append(f"Missing Comparison columns present in protein file: {', '.join(missing)}.")

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
