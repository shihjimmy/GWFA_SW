import re
import os

# Folders (可自行調整)
folder1     = "./GWFA_result_all"
folder2     = "./GWFA_result_all"
out_folder  = "./GWFA_result_compare"

os.makedirs(out_folder, exist_ok=True)

# List all files
files1 = sorted([f for f in os.listdir(folder1) if os.path.isfile(os.path.join(folder1, f))])
files2 = sorted([f for f in os.listdir(folder2) if os.path.isfile(os.path.join(folder2, f))])

# Pair files: file2 contains file1 but with _github added
paired_files = []
for f1 in files1:
    f1_base = os.path.splitext(f1)[0]
    matched = False
    for f2 in files2:
        f2_base = os.path.splitext(f2)[0]
        if re.fullmatch(f"{re.escape(f1_base)}_github.*", f2_base):
            paired_files.append((f1, f2))
            matched = True
            break
    if not matched:
        print(f"Warning: No matching file2 found for {f1}")

def read_nums_file1(path):
    """
    從 file1 讀出:
    - edit distances
    - case ids
    - traceback pass/fail 數量
    - ending pass/fail 數量
    """
    nums = []
    case_ids = []

    traceback_pass = 0
    traceback_fail = 0
    ending_pass = 0
    ending_fail = 0

    pat_num  = re.compile(r"final edit distance\s*:\s*(\d+)", re.IGNORECASE)
    pat_case = re.compile(r"Running test case pair:\s+.*?/([^/]+)\.fa\b", re.IGNORECASE)

    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line_stripped = line.strip()

            m = pat_num.search(line_stripped)
            if m:
                nums.append(float(m.group(1)))

            c = pat_case.search(line_stripped)
            if c:
                case_ids.append(c.group(1))

            low = line_stripped.lower()
            if "traceback pass" in low:
                traceback_pass += 1
            elif "traceback fail" in low:
                traceback_fail += 1

            if "ending pass" in low:
                ending_pass += 1
            elif "ending fail" in low:
                ending_fail += 1

    return nums, case_ids, traceback_pass, traceback_fail, ending_pass, ending_fail

def read_nums_file2(path):
    """從 file2 每行結尾抓數字，回傳 list[int]"""
    nums = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            m = re.search(r"(\d+)$", line.strip())
            if m:
                nums.append(int(m.group(1)))
    return nums

def smape(n1, n2):
    """SMAPE (%): 200 * |n1 - n2| / (|n1| + |n2|)；n1==n2==0 → 0"""
    if n1 == 0 and n2 == 0:
        return 0.0
    return 200.0 * abs(n1 - n2) / (abs(n1) + abs(n2))

def write_aligned_txt(rows, out_txt_path, smapes, f1_name, f2_name, tb_pass, tb_fail, end_pass, end_fail):
    headers = ["File1", "File2", "Case", "File1_Value", "File2_Value", "AbsError", "SMAPE(%)"]
    cols = list(zip(*([headers] + rows)))
    widths = [max(len(str(x)) for x in col) for col in cols]

    def fmt_row(r):
        return "  ".join(str(v).rjust(w) for v, w in zip(r, widths))

    with open(out_txt_path, 'w', encoding='utf-8') as f:
        f.write(fmt_row(headers) + "\n")
        f.write("  ".join("-" * w for w in widths) + "\n")
        for r in rows:
            f.write(fmt_row(r) + "\n")

        # Summary（SMAPE 版本）
        if smapes:
            summary = (
                f"[{f1_name} vs {f2_name}]  "
                f"SMAPE avg={sum(smapes)/len(smapes):.4f}%  "
                f"max={max(smapes):.4f}%  "
                f"min={min(smapes):.4f}%"
            )
        else:
            summary = f"[{f1_name} vs {f2_name}]  No SMAPE rows generated."
        f.write("\n" + summary + "\n")

        f.write(
            f"Traceback: pass={tb_pass}, fail={tb_fail}   |   "
            f"Ending: pass={end_pass}, fail={end_fail}\n"
        )

# 對齊：先算最大標籤寬度
pair_labels = [f"[{f1} vs {f2}]" for (f1, f2) in paired_files]
label_w = max((len(s) for s in pair_labels), default=0)

for f1_name, f2_name in paired_files:
    path1 = os.path.join(folder1, f1_name)
    path2 = os.path.join(folder2, f2_name)

    nums1, case_ids, tb_pass, tb_fail, end_pass, end_fail = read_nums_file1(path1)
    nums2 = read_nums_file2(path2)

    if len(nums1) != len(nums2):
        print(f"Warning: Length mismatch! {f1_name} file1={len(nums1)}, file2={len(nums2)}")
        continue

    n = min(len(nums1), len(nums2))
    rows = []
    smapes = []

    for i in range(n):
        n1 = nums1[i]
        n2 = nums2[i]

        if case_ids and i < len(case_ids) and n1 < n2:
            print(f"Should not be smaller than golden one: case={case_ids[i]}, {f1_name}={n1}, {f2_name}={n2}")

        abs_err = abs(n2 - n1)
        s = smape(n1, n2)
        smapes.append(s)

        rows.append([
            f1_name,
            f2_name,
            case_ids[i] if i < len(case_ids) else str(i),
            int(n1),
            int(n2),
            int(abs_err),
            f"{s:.6f}"
        ])

    f1_base = os.path.splitext(f1_name)[0]
    f2_base = os.path.splitext(f2_name)[0]
    out_txt = os.path.join(out_folder, f"{f1_base}__vs__{f2_base}.txt")

    write_aligned_txt(rows, out_txt, smapes, f1_name, f2_name, tb_pass, tb_fail, end_pass, end_fail)

    # 只印 SMAPE（對齊版）
    label = f"[{f1_name} vs {f2_name}]"
    if smapes:
        avg = sum(smapes) / len(smapes)
        print(
            f"{label:<{label_w}}  "
            f"SMAPE avg={avg:7.4f}%  "
            f"TB PASS={tb_pass:4d} FAIL={tb_fail:4d}  "
            f"END PASS={end_pass:4d} FAIL={end_fail:4d}"
        )
    else:
        print(
            f"{label:<{label_w}}  "
            f"No SMAPE rows generated.  "
            f"TB PASS={tb_pass:4d} FAIL={tb_fail:4d}  "
            f"END PASS={end_pass:4d} FAIL={end_fail:4d}"
        )
