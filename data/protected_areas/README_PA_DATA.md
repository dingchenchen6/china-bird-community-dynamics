# 中国自然保护区边界数据（分析用包）

## 数据来源与许可

China Nature Reserve Specimen Resource Sharing Platform (2024).
*List and Vector Boundaries of Nature Reserves in China* [Data set]. Zenodo.
https://doi.org/10.5281/zenodo.14875797

**许可：CC-BY-4.0**。转发与改编均允许，但必须保留上述署名。本目录中的
`china_nature_reserves_utf8.gpkg` 与 `*_list_utf8.csv` 是对原始数据的格式转换
（编码与容器格式），内容未作实质修改。

## 文件说明

| 文件 | 说明 | 推荐用途 |
|---|---|---|
| `china_nature_reserves_utf8.gpkg` | **推荐使用**。UTF-8 编码的 GeoPackage，图层名 `nature_reserves`，1028 个保护区，EPSG:4326，几何已修复（`buffer(0)`，无无效几何） | R/Python 分析主入口 |
| `china_nature_reserves_list_utf8.csv` | 保护区名录（3376 条）转成 UTF-8 CSV | 无 `openpyxl`/`readxl` 时使用 |
| `全国自然保护区名录+矢量边界/` | 原始解压数据（shp/dbf/shx/prj/xlsx/kmz） | 存档与回退 |
| `china_nature_reserves.rar` | 原始下载包 | 存档 |

## ⚠️ 编码问题（重要）

**原始 shapefile 是 GBK 编码，不是 UTF-8。** 在 Linux 服务器（UTF-8 locale）下
用 `st_read()` 直接读取会导致中文字段乱码，属性连接会**静默失败**而不报错——
结果看似跑通，实则错误。

因此：

1. **优先使用 `china_nature_reserves_utf8.gpkg`**（已转 UTF-8，跨平台安全）。
2. 若必须读原始 shapefile，显式声明编码：
   ```r
   st_read("保护区.shp", options = "ENCODING=GBK")
   ```
   本目录已补 `保护区.cpg`（内容为 `GBK`）以帮助 GDAL 自动识别。
3. `31_protected_area_effectiveness.R` 已实现：优先读 GeoPackage → 回退时显式
   声明 GBK → 并对名称字段做乱码自检，检测到 `�` 即停止运行。

## 字段说明（GeoPackage）

中英双份字段并存：

| 中文字段 | 英文字段 | 内容 |
|---|---|---|
| `保护区名称` | `name` | 保护区名称 |
| `级别` | `level` | 国家级 / 省级 / 市级 / 县级 / 州级 |
| `类型` | `type` | 野生动物 / 内陆湿地 / 海洋海岸 / 森林生态 等 |
| `始建时间` | `years` | 建立时间（YYYYMMDD） |
| `面积` | `area` | 面积（公顷） |
| `主管部门` | `dept` | 林业 / 环保 / 海洋 等 |
| `行政区域` | `adminarea` | 所属行政区 |
| `主保护对象` | `protect` | 主要保护对象描述 |

**级别分布**：国家级 544、省级 235、县级 171、市级 69、州级 2（合计 1028）。

## ⚠️ 覆盖度局限（写作时必须声明）

矢量边界仅 **1028** 条，而名录有 **3376** 条——即约 **30%** 的保护区有边界数据，
且明显偏向**国家级与大面积**保护区。因此：

- 保护成效结论偏向大型/高级别保护区，不能外推到全部保护地体系；
- 未被覆盖的中小型保护区会被误判为"非保护区"，这会使保护成效估计**偏保守**
  （即真实效应可能比估计值更强）；
- 论文方法与局限部分必须显式说明这一点。

## 上传到服务器

```bash
rsync -avz "$HOME/Documents/New project/bird_dynamic_occupancy_analysis/data/external/protected_areas/" \
  dingchenchen@<SERVER_IP>:"~/Documents/New project/bird_dynamic_occupancy_analysis/data/external/protected_areas/"
```

或上传打包文件后解压：

```bash
tar -xzf protected_areas_package.tar.gz -C ~/Documents/New\ project/bird_dynamic_occupancy_analysis/data/external/
```
