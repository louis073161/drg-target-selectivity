#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""从 GENCODE GTF 提取基因坐标（GRCh38），供基因层面 GWAS 汇总使用。"""
import gzip
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
GTF = ROOT / "gencode.v44.gtf.gz"
OUT = ROOT / "gene_coords_grch38.tsv"

rows = []
with gzip.open(GTF, "rt", encoding="utf-8") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "gene":
            continue
        attr = {}
        for kv in f[8].split(";"):
            kv = kv.strip()
            if not kv or " " not in kv:
                continue
            k, v = kv.split(" ", 1)
            attr[k] = v.strip('"')
        gtype = attr.get("gene_type", "")
        if gtype not in ("protein_coding", "lncRNA"):
            continue
        rows.append({
            "gene": attr.get("gene_name", ""),
            "ensembl": attr.get("gene_id", "").split(".")[0],
            "chrom": f[0].replace("chr", ""),
            "start": int(f[3]),
            "end": int(f[4]),
            "strand": f[6],
            "gene_type": gtype,
        })

df = pd.DataFrame(rows).drop_duplicates("gene")
df = df[df.chrom.isin([str(i) for i in range(1, 23)] + ["X"])]
df.to_csv(OUT, sep="\t", index=False)
print("基因记录 %d -> %s" % (len(df), OUT))
print(df.head().to_string(index=False))
