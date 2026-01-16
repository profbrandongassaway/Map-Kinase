
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PhosProtNorm — Embedded Selectors GUI (no popups, auto-detect ID columns)
- Main column pickers live on the main window (no modal dialogs)
- Auto-detect & pre-select: Uniprot ID, Gene Symbol, Site Position (overridable via dropdowns)
- Matching is by Uniprot ID only
- Supports TSV/TXT and XLSX inputs
- ALWAYS appends both Uniprot_site_ID and Gene_site_ID to the output (blanks if missing)
"""
import tkinter as tk
from tkinter import ttk, filedialog, font
import csv
import os

try:
    import pandas as pd
except Exception:
    pd = None  # XLSX support will be unavailable if pandas is not installed

# -------------------- UI helpers --------------------
def show_error_popup(parent, msg: str):
    win = tk.Toplevel(parent)
    win.title("Error")
    win.resizable(False, False)
    window_width, window_height = 420, 140
    screen_width = win.winfo_screenwidth()
    screen_height = win.winfo_screenheight()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    win.geometry(f"{window_width}x{window_height}+{x}+{y}")
    ttk.Label(win, text=msg, wraplength=380, justify="center").pack(pady=16, padx=12)
    ttk.Button(win, text="OK", command=win.destroy).pack(pady=6)

def show_success_popup(parent, msg="File has been successfully saved!"):
    popup = tk.Toplevel(parent)
    popup.title("Success")
    window_width, window_height = 280, 120
    screen_width = popup.winfo_screenwidth()
    screen_height = popup.winfo_screenheight()
    x = (screen_width - window_width) // 2
    y = (screen_height - window_height) // 2
    popup.geometry(f"{window_width}x{window_height}+{x}+{y}")
    ttk.Label(popup, text=msg, padding=10).pack(padx=20, pady=20)
    ttk.Button(popup, text="Close", command=popup.destroy).pack(pady=10)

def create_tooltip(widget, text):
    tooltip = tk.Toplevel(widget)
    tooltip.wm_overrideredirect(True)
    label = tk.Label(tooltip, text=text, background="white", relief="solid", borderwidth=1, font=("Arial", 8))
    label.pack()
    def show(event):
        x, y = widget.winfo_rootx(), widget.winfo_rooty() + widget.winfo_height()
        tooltip.wm_geometry(f"+{x}+{y}")
        tooltip.deiconify()
    def hide(event):
        tooltip.withdraw()
    tooltip.withdraw()
    widget.bind("<Enter>", show)
    widget.bind("<Leave>", hide)

# -------------------- Core helpers --------------------
def align_tab_fields(file_content: str):
    try:
        out = []
        lines = file_content.split('\n')
        if not lines:
            return None
        headers = lines[0].split('\t')
        header_count = len(headers)
        out.append(lines[0])
        for line in lines[1:]:
            row = line.split('\t')
            if len(row) < header_count:
                row.extend([''] * (header_count - len(row)))
            out.append('\t'.join(row))
        return '\n'.join(out)
    except Exception:
        return None

def search_column(file_content: str, column_name: str):
    try:
        headers = file_content.split('\n')[0].split('\t')
        column_name = column_name.lower().replace(' ', '_')
        for i, header in enumerate(headers):
            if column_name in header.lower().replace(' ', '_'):
                return i
        return -1
    except Exception:
        return -1

def phosscaled(file_content: str, column_indices):
    try:
        out = []
        lines = file_content.split('\n')
        if not lines or not lines[0].strip():
            return None
        headers = lines[0].split('\t')
        if max(column_indices, default=-1) >= len(headers):
            return None
        for col in column_indices:
            headers.append(f"{headers[col]}__PhosScaled")
        out.append('\t'.join(headers))
        for line in lines[1:]:
            row = line.split('\t')
            if len(row) <= max(column_indices, default=-1):
                out.append(line); continue
            values = []
            ok = True
            for col in column_indices:
                value = row[col].strip() if col < len(row) else ''
                try:
                    values.append(float(value) if value else 0.0)
                except (ValueError, TypeError):
                    ok = False; break
            if (not ok) or (not any(values)):
                out.append(line); continue
            row_sum = sum(values)
            if row_sum == 0:
                out.append(line); continue
            normalized_values = [(v / row_sum) * 100 for v in values]
            row.extend([str(round(val, 4)) for val in normalized_values])
            out.append('\t'.join(row))
        return '\n'.join(out)
    except Exception:
        return None

def append_protein_values(prot_content: str, phos_content: str, prot_column_indices,
                          prot_match_header=None, phos_match_header=None):
    try:
        prot_lines = prot_content.split('\n'); phos_lines = phos_content.split('\n')
        if not prot_lines or not phos_lines:
            return None
        prot_headers = prot_lines[0].split('\t'); phos_headers = phos_lines[0].split('\t')
        if max(prot_column_indices, default=-1) >= len(prot_headers):
            return None

        if not prot_match_header or not phos_match_header:
            return None
        if prot_match_header not in prot_headers or phos_match_header not in phos_headers:
            return None
        prot_match_col = prot_headers.index(prot_match_header)
        phos_match_col = phos_headers.index(phos_match_header)

        out = []
        out_headers = phos_headers + [f'prot_{prot_headers[i]}' for i in prot_column_indices]
        out.append('\t'.join(out_headers))

        prot_data_dict = {}
        for prot_line in prot_lines[1:]:
            prot_cells = prot_line.split('\t')
            if len(prot_cells) > prot_match_col:
                match_id = prot_cells[prot_match_col]
                if not match_id:
                    continue
                values = []
                for col in prot_column_indices:
                    try:
                        raw = prot_cells[col].strip()
                        if not raw:
                            values.append('')
                        else:
                            values.append(raw)
                    except (ValueError, IndexError):
                        values.append('')
                prot_data_dict[match_id] = values

        for phos_line in phos_lines[1:]:
            phos_cells = phos_line.split('\t')
            if len(phos_cells) > phos_match_col:
                match_id = phos_cells[phos_match_col]
                phos_cells.extend(prot_data_dict.get(match_id, [''] * len(prot_column_indices)))
            out.append('\t'.join(phos_cells))
        return '\n'.join(out)
    except Exception:
        return None

def phosprotnorm(file_content: str, phos_cols, prot_cols):
    try:
        if len(phos_cols) != len(prot_cols):
            return None

        lines = file_content.split('\n')
        headers = lines[0].split('\t')

        # Guard against bad indices
        if max(phos_cols + prot_cols, default=-1) >= len(headers):
            return None

        # Build names from the corresponding phospho headers
        new_norm_headers = []
        for ph_col in phos_cols:
            base = headers[ph_col]
            if base.endswith('__PhosScaled'):
                base = base[: -len('__PhosScaled')]
            elif base.lower().endswith('__phosscaled'):
                base = base[: -len('__phosscaled')]
            new_norm_headers.append(f"{base}_phosprotnorm")

        out = []

        # Append the new headers
        out_headers = headers + new_norm_headers
        out.append('\t'.join(out_headers))

        # Compute row-wise phos/protein normalization, preserving alignment
        for line in lines[1:]:
            row = line.split('\t')
            # If row is shorter than headers, pad to prevent index errors
            if len(row) < len(headers):
                row += [''] * (len(headers) - len(row))

            norm_values = []
            for phos_col, prot_col in zip(phos_cols, prot_cols):
                try:
                    phos_val = float(row[phos_col]) if row[phos_col].strip() else 0.0
                    prot_val = float(row[prot_col]) if row[prot_col].strip() else 0.0
                    if prot_val == 0:
                        norm_values.append('')  # avoid div-by-zero
                    else:
                        norm_values.append(str(round(phos_val / prot_val, 4)))
                except (ValueError, IndexError):
                    norm_values.append('')

            out.append('\t'.join(row + norm_values))

        return '\n'.join(out)
    except Exception:
        return None

def always_add_site_ids(file_content: str, uniprot_header: str, gene_header: str, sitepos_header: str):
    """
    ALWAYS adds both Uniprot_site_ID and Gene_site_ID columns.
    If source columns are missing or values are blank, the new cell is left blank.
    """
    try:
        lines = file_content.split('\n')
        headers = lines[0].split('\t')

        # indices (may be -1 if not found)
        def idx_of(hname):
            try:
                return headers.index(hname) if hname in headers else -1
            except ValueError:
                return -1

        uni_idx = idx_of(uniprot_header) if uniprot_header else -1
        gene_idx = idx_of(gene_header) if gene_header else -1
        site_idx = idx_of(sitepos_header) if sitepos_header else -1

        # Append new headers
        headers = headers + ["Uniprot_site_ID", "Gene_site_ID"]
        out = ['\t'.join(headers)]

        for line in lines[1:]:
            if not line:
                out.append(line); continue
            row = line.split('\t')
            # pad row if needed
            if len(row) < len(headers) - 2:  # minus the two we're adding
                row.extend([''] * ((len(headers) - 2) - len(row)))
            uni_val = row[uni_idx] if 0 <= uni_idx < len(row) else ''
            gene_val = row[gene_idx] if 0 <= gene_idx < len(row) else ''
            site_val = row[site_idx] if 0 <= site_idx < len(row) else ''

            uni_site = f"{uni_val}_{site_val}" if (uni_val and site_val) else ''
            gene_site = f"{gene_val}_{site_val}" if (gene_val and site_val) else ''

            row.extend([uni_site, gene_site])
            out.append('\t'.join(row))
        return '\n'.join(out)
    except Exception:
        return None

# -------------------- File IO --------------------
def _read_file_any(path: str):
    """
    Returns (headers: list[str], rows: list[list[str]]). All returned as strings.
    Supports TSV/TXT and XLSX (first sheet). XLSX requires pandas.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in [".tsv", ".txt"]:
        with open(path, newline='', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)
        if not rows:
            return [], []
        return rows[0], rows[1:]
    elif ext in [".xlsx", ".xls"]:
        if pd is None:
            raise RuntimeError("XLSX/XLS requires pandas. Please install pandas to read Excel files.")
        df = pd.read_excel(path, sheet_name=0, dtype=str)
        df = df.fillna('')
        headers = list(df.columns.map(str))
        rows = df.astype(str).values.tolist()
        return headers, rows
    else:
        with open(path, newline='', encoding='utf-8') as f:
            reader = csv.reader(f, delimiter='\t')
            rows = list(reader)
        if not rows:
            return [], []
        return rows[0], rows[1:]

def _to_tsv(headers, rows):
    return '\n'.join(['\t'.join(map(str, headers))] + ['\t'.join(map(str, r)) for r in rows])

# -------------------- Auto-detect heuristics --------------------
def _detect_candidates(headers):
    L = [h.lower() for h in headers]
    def find_first(patterns):
        for i, h in enumerate(L):
            for p in patterns:
                if p in h:
                    return headers[i]
        return None
    uniprot = find_first(["uniprot", "uniprot_id", "uniprot id", "accession", "protein accession", "protein id"])
    gene = find_first(["gene_symbol", "gene symbol", "gene", "genesymbol"])
    sitepos = find_first(["site position", "position in peptide", "site_position", "siteposition", "position"])
    return uniprot, gene, sitepos

# -------------------- GUI --------------------
class App:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("PhosProtNorm — Embedded Selectors")
        self.root.geometry("1000x640")
        self.root.minsize(900, 600)
        self.style = ttk.Style()
        self.style.configure("Bold.TLabel", font=("Helvetica", 12, "bold"))
        self.italic_font = font.Font(family="Helvetica", size=9, slant="italic")

        # State
        self.phos_path = tk.StringVar()
        self.prot_path = tk.StringVar()
        self.save_path = tk.StringVar()

        # Selected "Main columns"
        self.phos_main_cols = []
        self.prot_main_cols = []

        # Combo selections (overridable; will be auto-filled when files are loaded)
        self.phos_uniprot_col = tk.StringVar()
        self.phos_gene_col = tk.StringVar()
        self.phos_sitepos_col = tk.StringVar()
        self.prot_uniprot_col = tk.StringVar()
        self.prot_gene_col = tk.StringVar()

        self._build_ui()
        self.root.mainloop()

    def _build_ui(self):
        top = ttk.Frame(self.root, padding=10)
        top.pack(fill="both", expand=True)

        # ---- Left = Protein, Right = Phospho ----
        left = ttk.LabelFrame(top, text="Protein Input", padding=10)
        right = ttk.LabelFrame(top, text="Phospho Input", padding=10)
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
        right.grid(row=0, column=1, sticky="nsew", padx=(8, 0))
        top.columnconfigure(0, weight=1)
        top.columnconfigure(1, weight=1)
        top.rowconfigure(0, weight=1)

        # Protein section
        self._build_prot_section(left)

        # Phospho section
        self._build_phos_section(right)

        # Middle controls (matching + save)
        mid = ttk.Frame(self.root, padding=(10,0,10,8))
        mid.pack(fill="x")
        ttk.Label(mid, text="Matching: Uniprot ID", style="Bold.TLabel").grid(row=0, column=0, sticky="w", padx=4, pady=4)

        # Save footer
        footer = ttk.Frame(self.root, padding=10)
        footer.pack(fill="x")
        ttk.Button(footer, text="Choose Destination", command=self._choose_destination, width=18).grid(row=0, column=0, sticky="w")
        self.dest_entry = ttk.Entry(footer, textvariable=self.save_path)
        self.dest_entry.grid(row=0, column=1, sticky="ew", padx=8)
        self.save_btn = ttk.Button(footer, text="Save", command=self._on_save, width=16)
        self.save_btn.grid(row=0, column=2, sticky="e")
        footer.columnconfigure(1, weight=1)
        self._update_save_state()

    def _build_prot_section(self, frame):
        row = 0
        btn = ttk.Button(frame, text="Load Protein File", command=self._load_prot)
        btn.grid(row=row, column=0, sticky="w")
        self.prot_label = ttk.Label(frame, text="No file selected", font=self.italic_font)
        self.prot_label.grid(row=row, column=1, sticky="w", padx=8)
        row += 1

        ttk.Label(frame, text="Uniprot ID column:").grid(row=row, column=0, sticky="w", pady=(10, 2))
        self.prot_uniprot_combo = ttk.Combobox(frame, textvariable=self.prot_uniprot_col, state="readonly", width=40, values=[])
        self.prot_uniprot_combo.grid(row=row, column=1, sticky="w", pady=(10, 2))
        row += 1

        ttk.Label(frame, text="Gene Symbol column:").grid(row=row, column=0, sticky="w", pady=(2, 8))
        self.prot_gene_combo = ttk.Combobox(frame, textvariable=self.prot_gene_col, state="readonly", width=40, values=[])
        self.prot_gene_combo.grid(row=row, column=1, sticky="w", pady=(2, 8))
        row += 1

        ttk.Label(frame, text="Main Columns (Protein):", style="Bold.TLabel").grid(row=row, column=0, sticky="w", pady=(2, 4))
        row += 1
        self.prot_list = tk.Listbox(frame, selectmode="extended", height=14, width=40, exportselection=False)
        self.prot_list.grid(row=row, column=0, columnspan=2, sticky="nsew")
        frame.rowconfigure(row, weight=1)
        row += 1

        ttk.Label(frame, text="Tip: Ctrl/Cmd or Shift to multi-select.\nSelections are used as 'Main' protein columns.", font=self.italic_font).grid(row=row, column=0, columnspan=2, sticky="w", pady=(6,0))

    def _build_phos_section(self, frame):
        row = 0
        btn = ttk.Button(frame, text="Load Phospho File", command=self._load_phos)
        btn.grid(row=row, column=0, sticky="w")
        self.phos_label = ttk.Label(frame, text="No file selected", font=self.italic_font)
        self.phos_label.grid(row=row, column=1, sticky="w", padx=8)
        row += 1

        ttk.Label(frame, text="Uniprot ID column:").grid(row=row, column=0, sticky="w", pady=(10, 2))
        self.phos_uniprot_combo = ttk.Combobox(frame, textvariable=self.phos_uniprot_col, state="readonly", width=40, values=[])
        self.phos_uniprot_combo.grid(row=row, column=1, sticky="w", pady=(10, 2))
        row += 1

        ttk.Label(frame, text="Gene Symbol column:").grid(row=row, column=0, sticky="w", pady=(2, 2))
        self.phos_gene_combo = ttk.Combobox(frame, textvariable=self.phos_gene_col, state="readonly", width=40, values=[])
        self.phos_gene_combo.grid(row=row, column=1, sticky="w", pady=(2, 2))
        row += 1

        ttk.Label(frame, text="Site Position column:").grid(row=row, column=0, sticky="w", pady=(2, 8))
        self.phos_sitepos_combo = ttk.Combobox(frame, textvariable=self.phos_sitepos_col, state="readonly", width=40, values=[])
        self.phos_sitepos_combo.grid(row=row, column=1, sticky="w", pady=(2, 8))
        row += 1

        ttk.Label(frame, text="Main Columns (Phospho):", style="Bold.TLabel").grid(row=row, column=0, sticky="w", pady=(2, 4))
        row += 1
        self.phos_list = tk.Listbox(frame, selectmode="extended", height=14, width=40, exportselection=False)
        self.phos_list.grid(row=row, column=0, columnspan=2, sticky="nsew")
        frame.rowconfigure(row, weight=1)
        row += 1

        ttk.Label(frame, text="Tip: Ctrl/Cmd or Shift to multi-select.\nSelections are used as 'Main' phospho columns.", font=self.italic_font).grid(row=row, column=0, columnspan=2, sticky="w", pady=(6,0))

    # -------------------- Events --------------------
    def _choose_destination(self):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".tsv",
            filetypes=[("TSV", "*.tsv"), ("Tab-Delimited", "*.txt *.tsv"), ("All files", "*.*")]
        )
        if file_path:
            self.save_path.set(file_path)
            self._update_save_state()

    def _load_prot(self):
        path = filedialog.askopenfilename(
            title="Select Protein file (TSV/TXT/XLSX)",
            filetypes=(("TSV/Tab/Excel", "*.tsv *.txt *.xlsx *.xls"), ("All Files", "*.*"))
        )
        if not path:
            return
        try:
            headers, rows = _read_file_any(path)
        except Exception as e:
            show_error_popup(self.root, f"Could not read Protein file:\n{e}")
            return

        self.prot_path.set(path)
        disp = path if len(path) < 70 else path[:70] + "..."
        self.prot_label.configure(text=disp)
        create_tooltip(self.prot_label, path)

        self.prot_uniprot_combo['values'] = headers
        self.prot_gene_combo['values'] = headers
        self.prot_list.delete(0, tk.END)
        for h in headers:
            self.prot_list.insert(tk.END, h)

        uni, gene, _ = _detect_candidates(headers)
        if uni: self.prot_uniprot_col.set(uni)
        if gene: self.prot_gene_col.set(gene)

        self._update_save_state()

    def _load_phos(self):
        path = filedialog.askopenfilename(
            title="Select Phospho file (TSV/TXT/XLSX)",
            filetypes=(("TSV/Tab/Excel", "*.tsv *.txt *.xlsx *.xls"), ("All Files", "*.*"))
        )
        if not path:
            return
        try:
            headers, rows = _read_file_any(path)
        except Exception as e:
            show_error_popup(self.root, f"Could not read Phospho file:\n{e}")
            return

        self.phos_path.set(path)
        disp = path if len(path) < 70 else path[:70] + "..."
        self.phos_label.configure(text=disp)
        create_tooltip(self.phos_label, path)

        self.phos_uniprot_combo['values'] = headers
        self.phos_gene_combo['values'] = headers
        self.phos_sitepos_combo['values'] = headers

        self.phos_list.delete(0, tk.END)
        for h in headers:
            self.phos_list.insert(tk.END, h)

        uni, gene, site = _detect_candidates(headers)
        if uni: self.phos_uniprot_col.set(uni)
        if gene: self.phos_gene_col.set(gene)
        if site: self.phos_sitepos_col.set(site)

        self._update_save_state()

    def _collect_main_cols(self):
        self.prot_main_cols = [self.prot_list.get(i) for i in self.prot_list.curselection()]
        self.phos_main_cols = [self.phos_list.get(i) for i in self.phos_list.curselection()]

    def _update_save_state(self):
        self._collect_main_cols()
        ok = bool(self.phos_path.get() and self.prot_path.get() and self.save_path.get()
                  and self.prot_main_cols and self.phos_main_cols
                  and self.phos_uniprot_col.get()
                  and self.prot_uniprot_col.get()
                 )
        try:
            self.save_btn.state(["!disabled"] if ok else ["disabled"])
        except Exception:
            pass

    def _on_save(self):
        self._collect_main_cols()
        if not (self.phos_path.get() and self.prot_path.get() and self.save_path.get() and self.prot_main_cols and self.phos_main_cols):
            show_error_popup(self.root, "Please load both files, pick 'Main Columns' for each, and choose a destination.")
            return

        try:
            phos_headers, phos_rows = _read_file_any(self.phos_path.get())
            prot_headers, prot_rows = _read_file_any(self.prot_path.get())
        except Exception as e:
            show_error_popup(self.root, f"Failed to re-read inputs:\n{e}")
            return

        phos_tsv = _to_tsv(phos_headers, phos_rows)
        prot_tsv = _to_tsv(prot_headers, prot_rows)

        phos_tsv = align_tab_fields(phos_tsv)
        prot_tsv = align_tab_fields(prot_tsv)
        if phos_tsv is None or prot_tsv is None:
            show_error_popup(self.root, "Failed to align fields. Are the inputs valid?")
            return

        def headers_to_indices(headers, chosen):
            return [headers.index(h) for h in chosen if h in headers]

        phos_indices = headers_to_indices(phos_tsv.split('\n')[0].split('\t'), self.phos_main_cols)
        prot_indices = headers_to_indices(prot_tsv.split('\n')[0].split('\t'), self.prot_main_cols)

        if not phos_indices or not prot_indices:
            show_error_popup(self.root, "No valid 'Main Columns' detected. Make sure your selections exist in the files.")
            return

        if len(phos_indices) != len(prot_indices):
            show_error_popup(self.root, "Phospho and Protein 'Main Columns' must be the same count for Phos/Prot.")
            return

        # step 1: match by Uniprot and append raw protein values
        pre_cols = len(phos_tsv.split('\n')[0].split('\t'))
        phos_tsv = append_protein_values(
            prot_tsv,
            phos_tsv,
            prot_indices,
            prot_match_header=self.prot_uniprot_col.get(),
            phos_match_header=self.phos_uniprot_col.get()
        )
        if phos_tsv is None:
            show_error_popup(self.root, "Matching failed. Ensure Uniprot ID columns exist in BOTH files.")
            return
        post_cols = len(phos_tsv.split('\n')[0].split('\t'))
        prot_cols = list(range(pre_cols, post_cols))

        # step 2: phosprotnorm (raw Phos / raw Prot)
        phos_tsv_final = phosprotnorm(phos_tsv, phos_indices, prot_cols)
        if phos_tsv_final is None:
            show_error_popup(self.root, "phosprotnorm failed. Check that counts of phospho/protein 'Main Columns' match.")
            return

        # step 3: ALWAYS add both site IDs (using PHOSPHO-side headers the user selected/auto-detected)
        phos_header_list = phos_tsv_final.split('\n')[0].split('\t')
        # We pass the dropdown values (may or may not exist; helper will still add columns)
        phos_tsv_final = always_add_site_ids(
            phos_tsv_final,
            uniprot_header=self.phos_uniprot_col.get(),
            gene_header=self.phos_gene_col.get(),
            sitepos_header=self.phos_sitepos_col.get()
        )

        if phos_tsv_final is None:
            show_error_popup(self.root, "Adding Uniprot_site_ID / Gene_site_ID failed.")
            return

        try:
            with open(self.save_path.get(), 'w', encoding='utf-8') as outf:
                outf.write(phos_tsv_final)
        except Exception as e:
            show_error_popup(self.root, f"Failed to write output:\n{e}")
            return

        show_success_popup(self.root, "Saved TSV successfully.")

# -------------------- main --------------------
if __name__ == "__main__":
    App()
