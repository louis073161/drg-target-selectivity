# 发布指南：GitHub + Zenodo（30 分钟可完成）

目标：让稿件里能写「代码与派生数据见 GitHub，永久归档于 Zenodo，DOI: 10.5281/zenodo.XXXXXXX」，
而不是「可向作者索取」。这是高分刊的硬性要求，也是审稿人判断可信度的第一条依据。

---

## 0. 三个账号（只做一次）

| 账号 | 用途 | 说明 |
|---|---|---|
| GitHub | 放代码，审稿人点开就能看 | 免费，public 仓库无限私有额度 |
| ORCID | Zenodo 认人，DOI 元数据里绑 ORCID | 你已有，直接用 |
| Zenodo | 发 DOI，CERN 运营，永久免费 | 用 GitHub 账号一键登录，自动关联 |

> Zenodo 不需要单独注册：官网点 Sign in → GitHub → 授权即可。

## 1. 本地已备好的东西

`drg-target-selectivity/` 这个目录就是完整的仓库，已包含：

- 31 个分析脚本（`sc_demo/`，原样保留，未做任何"美化改写"）
- 79 个派生数据文件（`data/01…12/`，17.5 MB），每张表每张图都能回溯到一个 CSV
- 7 张成品图（`figures/`）
- `README.md`、`LICENSE`（MIT）、`CITATION.cff`、`requirements.txt`、`.gitignore`
- `tools/quick_check.py`：审稿人一条命令即可复核表 10、3.13、3.11、图 2，逐字节比对 MD5

合计约 28 MB。GitHub 单文件上限 100 MB、仓库建议 1 GB 以内，完全够。

**尚未完成的两处占位**，发布前必须改：

1. `CITATION.cff` 里的 `repository-code`、`orcid`；
2. Zenodo 生成 DOI 后，把 `10.5281/zenodo.XXXXXXX` 替换成真实 DOI（三处：`CITATION.cff` 两处、`README.md` 一处）。

## 2. 上传到 GitHub

### 方式 A：命令行（最快，本机已装 git 2.55）

```bash
cd drg-target-selectivity
git init -b main
git add -A
git commit -m "Code and derived data for DRG target selectivity benchmark (v1.0.0)"

# 先在网页上建一个空仓库 drg-target-selectivity（不要勾选 README/.gitignore）
git remote add origin https://github.com/<你的账号>/drg-target-selectivity.git
git push -u origin main
```

首次 push 会弹出浏览器授权，登录即可。国内网络若超时，改用下面的方式 B。

### 方式 B：GitHub Desktop（图形界面，无需命令行）

1. 装 GitHub Desktop → File → Add local repository → 选 `drg-target-selectivity`；
   提示 "not a git repository" 时点 **create a repository** → Publish repository；
2. 勾掉 "Keep this code private"，点 Publish。

### 方式 C：网页拖拽（零安装，但不推荐）

新建仓库 → Add file → Upload files → 把文件夹拖进去。缺点：文件夹结构要手动保持，
`.gitignore` 不生效，后续改版麻烦。

## 3. 打一个版本标签（Zenodo 靠 tag 触发归档）

GitHub 网页：右侧 **Releases → Create a new release**

- Tag: `v1.0.0`（必须新建 tag，不要只写标题）
- Title: `v1.0.0 — manuscript submission`
- 描述里写一句：`Code and derived data matching the 2026-11-09 manuscript.`

此时**先不要点 Publish**，先把 Zenodo 开关打开（第 4 步），否则这次 release 不会被归档。

## 4. 打开 Zenodo 的 GitHub 同步

1. 登录 https://zenodo.org → 右上头像 → **Settings** → 左侧 **GitHub**；
2. 点 **Connect**，授权 Zenodo 读取你的公开仓库；
3. 回到列表，找到 `drg-target-selectivity`，把开关拨到 **ON**。

## 5. 发布 Release，自动拿到 DOI

回到 GitHub 那个 draft release，点 **Publish release**。几分钟后：

- Zenodo 收件箱（https://zenodo.org/account/settings/github/ 下方）出现新条目；
- 条目页右上角有两个 DOI：
  - **Version DOI**（`10.5281/zenodo.1234567`）—— 指向 v1.0.0 这一版，**写进稿件用这个**；
  - **Concept DOI**（`10.5281/zenodo.1234566`）—— 永远指向最新版，适合放在个人主页。

把 Version DOI 回填到 `CITATION.cff`、`README.md`，再 commit 一次即可（这一改动不会再触发新 DOI）。

## 6. 3.1 GB 的原始归档怎么办

GitHub 放不下，但 Zenodo 单条记录上限 **50 GB**，正好。做法：

1. 把 `原始数据/` 压成一个 zip（排除 `11_未复制_大体积原始件` 里 282 GB 的 `timsTOF.zip`，
   它本来就只有登记没有实体）：
   ```bash
   cd 原始数据
   zip -r ../drg-archive-2026-11-09.zip . -x "11_*/*" -x "*.h5ad"
   ```
   `.h5ad` 是否排除取决于你：含进去约 2.2 GB，能让 `power_sim.py` 等脚本完全复跑，
   **建议含进去**（审稿人最想跑的就是 3.12 节）。
2. Zenodo → **New upload** → 拖入 zip → 填元数据：
   - Upload type: `Dataset`
   - Title: `Archive of intermediates for the human DRG target selectivity benchmark`
   - Authors: 同上
   - Description: 一段话说明这是稿件配套的中间结果与 `.h5ad`，并注明与 GitHub 仓库的关系
   - License: `CC BY 4.0`
   - 在 **Related identifiers** 里填 GitHub 仓库 URL，关系选 `is supplemented by`
3. Publish → 得第二个 DOI。两个 DOI 都要写进稿件（一个是代码，一个是数据）。

## 7. 之后要改代码怎么办

**不要**删旧记录重发。做法：改代码 → 新 tag `v1.1.0` → 新 release → Zenodo 自动生成新
version DOI，Concept DOI 不变。稿件已发表就引 v1.0.0 的 DOI，不必动。

## 8. 稿件里那句话怎么写

中文（放在"数据可用性"一节末尾）：

> 全部分析脚本、派生数据表与图件见 GitHub（https://github.com/&lt;账号&gt;/drg-target-selectivity，
> MIT 许可），并永久归档于 Zenodo（DOI: 10.5281/zenodo.XXXXXXX）；中间结果与 `.h5ad`
> 对象归档于 Zenodo（DOI: 10.5281/zenodo.YYYYYYY）。仓库内 `tools/quick_check.py`
> 可一键复核表 10、3.11、3.13 与图 2，输出与本稿逐字节一致。原始测序数据来自
> GSE249746、GSE59739、GSE310000、GSE201586，遗传学数据来自 FinnGen R12，
> 蛋白组数据来自 PRIDE PXD061129，均按各自许可公开。

英文：

> All analysis scripts, derived tables and figure source files are available at
> https://github.com/&lt;account&gt;/drg-target-selectivity (MIT) and are permanently archived at
> Zenodo (DOI: 10.5281/zenodo.XXXXXXX); intermediates and `.h5ad` objects are archived
> separately (DOI: 10.5281/zenodo.YYYYYYY). `tools/quick_check.py` reproduces Table 10,
> sections 3.11 and 3.13, and Figure 2 byte-for-byte. Raw sequencing data were obtained
> from GSE249746, GSE59739, GSE310000 and GSE201586, genetic associations from FinnGen
> release 12, and proteomic data from PRIDE PXD061129, each under its own terms.

## 9. 四条容易踩的坑

1. **DOI 不可撤销。** Zenodo 一旦 publish，记录永久存在、不能删除（只能限制访问）。
   所以 publish 前务必确认文件是最终版、没有患者信息、没有未发表的他方数据。
2. **不要在仓库里放患者级数据。** 本文全部为公开数据重算，不存在此问题；但若以后加入
   自采样本，只放汇总表。
3. **`.gitignore` 要先于第一次 commit。** 本仓库已配好（`sc_demo/*.csv` 等运行产物不入库），
   千万别 `git add -f` 把 216 MB 的 `.h5ad` 塞进去 —— GitHub 会拒绝，且历史里删不干净。
4. **匿名评审。** 少数刊要求双盲，此时 GitHub 仓库先设 private、Zenodo 记录设
   "Restricted"，录用后开放；或使用 https://anonymous.4open.science/ 生成匿名镜像。
