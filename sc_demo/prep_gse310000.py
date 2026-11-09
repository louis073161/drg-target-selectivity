#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把 GSE310000_RAW.tar 展开成 scanpy 可读的 h5ad。

GEO 的 RAW.tar 里每个样本是 10x 三件套，但文件名带 GSM 前缀，
read_10x_mtx 要求目录内固定叫 matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz，
所以这里按样本建子目录并改名，再合并。
"""
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

ROOT = Path(__file__).resolve().parent
TAR = ROOT / "GSE310000_RAW.tar"
OUT = ROOT / "GSE310000.h5ad"


def main():
    if not TAR.exists():
        sys.exit("找不到 %s" % TAR)
    raw = ROOT / "raw"
    raw.mkdir(exist_ok=True)
    with tarfile.open(TAR) as tf:
        names = tf.getnames()
        tf.extractall(raw)
    print("[解包] %d 个文件" % len(names))

    groups = {}
    for n in names:
        m = re.match(r"(GSM\d+)_(.+?)_(barcodes|features|matrix)\.(tsv|mtx)\.gz$", Path(n).name)
        if m:
            groups.setdefault((m.group(1), m.group(2)), {})[m.group(3)] = raw / n

    ads, labels = [], []
    for (gsm, samp), files in sorted(groups.items()):
        if len(files) != 3:
            print("[跳过] %s_%s 缺少文件：%s" % (gsm, samp, list(files)))
            continue
        d = raw / ("%s_%s" % (gsm, samp))
        d.mkdir(exist_ok=True)
        for k, src in files.items():
            dst = d / {"barcodes": "barcodes.tsv.gz", "features": "features.tsv.gz",
                       "matrix": "matrix.mtx.gz"}[k]
            if not dst.exists():
                dst.write_bytes(src.read_bytes())
        ad = sc.read_10x_mtx(d, var_names="gene_symbols", make_unique=True, cache=False)
        ad.obs["sample"] = samp
        ad.obs_names = [gsm + "_" + b for b in ad.obs_names]
        print("[读入] %-8s %5d 细胞 × %5d 基因" % (samp, ad.n_obs, ad.n_vars))
        ads.append(ad)

    if not ads:
        sys.exit("没有可用的样本")

    ad = ads[0] if len(ads) == 1 else sc.concat(ads, join="outer", fill_value=0)
    ad.obs["sample"] = ad.obs["sample"].astype("category")
    ad.var_names_make_unique()

    X = ad.X
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    print("[合计] %d 细胞 × %d 基因；每细胞计数中位数 %.0f"
          % (ad.n_obs, ad.n_vars, float(np.median(np.asarray(X.sum(1)).ravel()))))
    print("[样本] %s" % ad.obs["sample"].value_counts().to_dict())

    ad.write_h5ad(OUT, compression="gzip")
    print("[写出] %s" % OUT)


if __name__ == "__main__":
    main()
