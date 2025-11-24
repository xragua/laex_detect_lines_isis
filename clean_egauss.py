#!/usr/bin/env python3
"""
clean_egauss.py

Usage:
    python clean_egauss.py input.par output.par

What it does:
- Finds egauss(N) parameter triplets (area, center, sigma)
- Drops any triplet where area == 0 (numerically)
- Renumbers the remaining egauss indices in ascending order starting from 1
- Renumbers the leading "idx" column sequentially starting at 1
- Preserves column formatting as closely as practical
"""
import sys
import re
from pathlib import Path

AREA_RE = re.compile(r'(egauss\()\s*(\d+)\s*(\)\.area)')
CENTER_RE = re.compile(r'(egauss\()\s*(\d+)\s*(\)\.center)')
SIGMA_RE = re.compile(r'(egauss\()\s*(\d+)\s*(\)\.sigma)')

def parse_float_safe(s: str):
    try:
        return float(s)
    except Exception:
        return None

def is_egauss_area(line: str):
    return bool(AREA_RE.search(line))

def get_area_value(line: str):
    """
    Extract the 'value' field from a parameter line like:
      "  45  egauss(1).area         0     0     9.613535e-06           0       1e+08  photons/s/cm^2"
    We need the 5th column (value). We'll split by whitespace and read the column positions:
    idx  name  tie  freeze  value  min  max  [unit]
    """
    # Collapse multiple spaces, but we must be careful not to lose units.
    parts = line.strip().split()
    # We expect at least 7 columns before unit; if unit exists it's 8+
    if len(parts) < 7:
        return None
    # parts[0]=idx, parts[1]=name, parts[2]=tie-to, parts[3]=freeze, parts[4]=value
    return parse_float_safe(parts[4])

def renumber_idx_column(lines):
    """Rewrite the first numeric column (idx) to be 1..N, preserving right alignment width (3 digits or more)."""
    out = []
    idx_counter = 1
    idx_field_width = 3  # observed in sample
    for line in lines:
        m = re.match(r'\s*\d+\s+(.*)$', line)
        if not m:
            out.append(line)
            continue
        rest = m.group(1)
        out.append(f"{idx_counter:4d}  {rest}")
        idx_counter += 1
    return out

def main(inp: Path, outp: Path):
    text = inp.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    # Keep header lines (typically first 2 lines). We'll process from line 3 onward.
    if len(lines) < 3:
        outp.write_text(text, encoding="utf-8")
        return

    header0 = lines[0]
    header1 = lines[1]
    body = lines[2:]

    # Walk through body; identify egauss triplets and decide whether to keep them.
    kept_lines = []
    i = 0
    # We'll also track renumber mapping old_index -> new_index
    old_to_new = {}
    next_egauss_index = 1

    while i < len(body):
        line = body[i]
        if is_egauss_area(line):
            # Expect the following two lines to be center and sigma with the same index
            if i + 2 >= len(body):
                # malformed; just keep as-is
                kept_lines.append(line)
                i += 1
                continue
            area_line = body[i]
            center_line = body[i+1]
            sigma_line = body[i+2]

            # Check they belong to the same egauss index
            ma = AREA_RE.search(area_line)
            mc = CENTER_RE.search(center_line) if CENTER_RE.search(center_line) else None
            ms = SIGMA_RE.search(sigma_line) if SIGMA_RE.search(sigma_line) else None

            same_index = False
            if ma and mc and ms:
                ia = int(ma.group(2))
                ic = int(mc.group(2))
                isg = int(ms.group(2))
                same_index = (ia == ic == isg)

            if same_index:
                val = get_area_value(area_line)
                if val is not None and abs(val) == 0.0:
                    # Drop this triplet
                    i += 3
                    continue
                else:
                    # Keep this triplet but rewrite egauss(index) to new contiguous index
                    old_idx = int(ma.group(2))
                    new_idx = next_egauss_index
                    next_egauss_index += 1
                    old_to_new[old_idx] = new_idx

                    def renumber_line(l, pat):
                        return pat.sub(lambda m: f"{m.group(1)}{new_idx}{m.group(3)}", l)

                    kept_lines.append(renumber_line(area_line, AREA_RE))
                    kept_lines.append(renumber_line(center_line, CENTER_RE))
                    kept_lines.append(renumber_line(sigma_line, SIGMA_RE))
                    i += 3
                    continue
            else:
                # Not a well-formed triplet; just keep and move forward
                kept_lines.append(line)
                i += 1
                continue
        else:
            kept_lines.append(line)
            i += 1

    # Now renumber the idx column from 1..N within the body
    body_renumbered = renumber_idx_column(kept_lines)

    # Build new header with all remaining egauss components explicitly
    area_pat = re.compile(r'egauss\((\d+)\)\.area')
    idxs = sorted({int(m.group(1)) for m in map(lambda mm: re.search(area_pat, mm), kept_lines) if m})
    if idxs:
        egauss_terms = "+".join(f"egauss({i})" for i in idxs)
        header0 = f"tbnew(1)*(powerlaw(1)+{egauss_terms})"
    else:
        header0 = "tbnew(1)*(powerlaw(1))"


    # Reassemble
    final_lines = [header0, header1] + body_renumbered
    outp.write_text("\n".join(final_lines) + "\n", encoding="utf-8")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python clean_egauss.py input.par output.par")
        sys.exit(1)
    inp = Path(sys.argv[1])
    outp = Path(sys.argv[2])
    main(inp, outp)
