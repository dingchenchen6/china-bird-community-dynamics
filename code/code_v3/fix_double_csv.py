#!/usr/bin/env python3
"""把 *.csv.csv 收敛为 *.csv。
安全规则：目标已存在时，比较两者——内容相同则删冗余；不同则保留较新的一份，
旧的改名加 .superseded 后缀，绝不静默覆盖。
Collapse *.csv.csv to *.csv without silently overwriting an existing file."""
import filecmp
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2] / "results"
renamed = same = superseded = 0

for src in sorted(ROOT.rglob("*.csv.csv")):
    dst = src.with_suffix("")            # 去掉最后一个 .csv
    if not dst.exists():
        src.rename(dst); renamed += 1
        continue
    if filecmp.cmp(src, dst, shallow=False):
        src.unlink(); same += 1          # 内容一致，删冗余
        continue
    # 内容不同：保留较新者
    if src.stat().st_mtime >= dst.stat().st_mtime:
        dst.rename(dst.with_suffix(".csv.superseded"))
        src.rename(dst)
    else:
        src.rename(src.with_suffix(".csv.superseded"))
    superseded += 1

print(f"改名 {renamed} 个；内容一致删冗余 {same} 个；版本冲突保留两份 {superseded} 个")
rest = list(ROOT.rglob("*.csv.csv"))
print("剩余 .csv.csv:", len(rest))
