from importlib import metadata
# from cell_sketcher.utils import (
#     get_package_dir,
#     base64uri,
#     csv2dict,
#     help_collapse,
#     read_mtx,
# )
from space_sketcher.rna.src.plotly_to_report import (
    plot_barcode_ranks,
    saturation_plot,
    distribution_violin,
    spatial_scatter,
    plot_sb_cb_umi_knee
)
from pathlib import Path
from jinja2 import FileSystemLoader, Environment
import pandas as pd
import numpy as np
import plotly.express as px
import scanpy as sc
from space_sketcher.tools.utils import (
    csv2dict,
    get_package_dir,
    base64uri,
    help_collapse,
    read_mtx
)  

def guess_type(filepath):
    """
    类型自动猜测
    """
    import mimetypes

    return mimetypes.guess_type(filepath)[0]

def file_to_base64(filepath):
    """
    文件转 base64
    """
    import base64

    with open(filepath, "rb") as f:
        encoded_str = base64.b64encode(f.read())
    return encoded_str.decode("utf-8")

def base64uri(filepath):
    """
    文件转为 base64uri
    """
    return "data:%s;base64,%s" % (guess_type(filepath), file_to_base64(filepath))


def help_icon():
    """
    帮助图标
    """
    icon = """<svg t="1687338260183" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2757" width="18" height="18" data-spm-anchor-id="a313x.7781069.0.i9"><path d="M512 0C229.23 0 0 229.23 0 512s229.23 512 512 512 512-229.23 512-512S794.77 0 512 0zM512 928c-229.75 0-416-186.25-416-416S282.25 96 512 96s416 186.25 416 416S741.75 928 512 928z" p-id="2758" fill="#707070"></path><path d="M537.64 343.452c47.074 0 83.266-37.528 83.266-78.072 0-32.46-20.832-60.878-62.496-60.878-54.816 0-82.178 44.618-82.178 77.11C475.144 320.132 498.152 343.452 537.64 343.452z" p-id="2759" fill="#707070"></path><path d="M533.162 728.934c-7.648 0-10.914-10.136-3.264-39.55l43.25-166.406c16.386-60.848 10.944-100.398-21.92-100.398-39.456 0-131.458 39.83-211.458 107.798l16.416 27.392c25.246-17.256 67.906-34.762 77.792-34.762 7.648 0 6.56 10.168 0 35.508l-37.746 158.292c-23.008 89.266 1.088 109.538 33.984 109.538 32.864 0 117.808-30.47 195.57-109.632l-18.656-25.34C575.354 716.714 543.05 728.934 533.162 728.934z" p-id="2760" fill="#707070"></path></svg>"""
    return icon


def pct(n):
    """
    百分比 without zero
    """
    return "%g%%" % (n * 100)


class Counter:
    """
    用于在html 中生成唯一id，不使用uuid
    """

    def __init__(self, start=0):
        self.count = start

    def increment(self):
        self.count += 1
        return self.count


_ID = Counter(0)

def help_collapse(help_msg: list | None = None):
    """
    生成折叠的帮助文档
    """
    if help_msg is None:
        return {"icon": "", "text": ""}
    id = f"collapse-{_ID.increment()}"
    if isinstance(help_msg, list):
        help_msg_str = ""
        for i in help_msg:
            if isinstance(i, str):
                help_msg_str += f"<p>{i}</p>"
            else:
                help_msg_str += (
                    f'<h6><strong>{i["name"]}</strong></h6><p>{i["value"]}</p>'
                )
    else:
        help_msg_str = help_msg
    msg = {
        "icon": f' <a data-bs-toggle = "collapse" href="#{id}" aria-expanded="false" aria-controls="help">{help_icon()}</a>',
        "text": f'<div class = "collapse" id={id}>{help_msg_str}</div>',
    }
    return msg


def report_card(
    title: str = None,
    help_msg: list | None = None,
    content: str = "",
    style="margin-top: 1.6rem;",
):
    """生成报告的卡片
    to do: 是否使用 htmltools，目前代码有点难以阅读

    Args:
        title (str): 标题
        help_msg (list | None): 帮助文档
        content (str): html 的内容，可以是表格，图片等
    """
    help_msg = help_collapse(help_msg)
    if title is None:
        id = ""
        title = ""
    else:
        id = f'id="{title}"'
        title = f'<h5 class = "fs-2 card-title fw-bold">{title}{help_msg["icon"]}</h5>'

    html_str = (
        f'<div style = "{style}" class = "card" {id}>'
        '<div style="padding: var(--bs-card-spacer-y) var(--bs-card-spacer-x) 0 var(--bs-card-spacer-x);">'
        f'{title}'
        f'{help_msg['text']}'
        "</div>"
        f"{content}</div>"
    )
    return html_str


def to_table_html(values: dict):
    """
    表格转 html
    """
    table_str = ""
    for k, v in values.items():
        table_str += f'<tr><td class="gt_row gt_left">{k}</td><td class="gt_row gt_right">{v}</td></tr>'

    return f"""
    <table class="gt_table">
    <tbody class="gt_table_body">
        {table_str}
    </tbody>
    </table>
    """


def tbl(df: pd.DataFrame, id=None):
    """
    数据框转 html
    """
    if id is None:
        id = ""
    else:
        id = f'id="{id}"'
    t_body = ""
    for i, rows in df.iterrows():
        td_str = ""
        for value in rows:
            td_str += f'<td class="gt_row gt_left">{value}</td>'
        t_body += f"<tr>{td_str}</tr>"
    t_head = ""
    for column in df.columns:
        t_head += f'<th class="gt_left" scope="col">{column}</th>'
    html_str = f"""
<table {id} class="gt_table">
<thead>
<tr class="header gt_col_headings">
{t_head}
</tr>
</thead>
<tbody class="gt_table_body">
{t_body}
</tbody>
</table>
    """
    return html_str


_version = metadata.version("cell_sketcher")
source_dir = get_package_dir()


def tbl(df: pd.DataFrame, id=None):
    """
    数据框转 html
    """
    if id is None:
        id = ""
    else:
        id = f'id="{id}"'
    t_body = ""
    for i, rows in df.iterrows():
        td_str = ""
        for value in rows:
            td_str += f"<td>{value}</td>"
        t_body += f"<tr>{td_str}</tr>"
    t_head = ""
    for column in df.columns:
        t_head += f"<th>{column}</th>"
    html_str = f"""
<table {id} class="table table-striped">
<thead>
<tr class="header">
{t_head}
</tr>
</thead>
<tbody">
{t_body}
</tbody>
</table>
    """
    return html_str

def to_bootstrap_table(values: dict):
    """
    表格转 html
    """
    table_str = ""
    for k, v in values.items():
        table_str += f'<tr><td class="gt_row gt_left">{k}</td><td class="gt_row gt_right">{v}</td></tr>'

    return f"""
    <table class="gt_table">
    <tbody class="gt_table_body">
        {table_str}
    </tbody>
    </table>
    """

def report_card(
    title: str = None, help_msg: list | None = None, content: str = "", style: str = ""
):
    """
    生成报告的卡片

    Args:
        title (str): 标题
        help_msg (list | None): 帮助文档
        content (str): html 的内容，可以是表格，图片等
    """
    help_msg = help_collapse(help_msg)
    if title is None:
        id = ""
        title = ""
    else:
        id = f'id="{title}"'
        title = f'<div class = "card-header fw-bold">{title}{help_msg["icon"]}</div>'

    html_str = f"""
<div class="card" {id} style="{style}">
    {title}
    <div class="card-body">
        {help_msg["text"]}
        {content}
    </div>
</div>
"""
    return html_str

def get_resource(cdn=False):
    """获取css 以及 js等资源"""
    if cdn:
        links = [
            "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css",
            "https://cdn.jsdelivr.net/npm/@jstable/jstable@1.6.5/dist/jstable.min.css"
        ]
        scripts = [
            "https://cdn.plot.ly/plotly-3.0.0.min.js",
            "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/js/bootstrap.min.js",
            "https://cdn.jsdelivr.net/npm/@jstable/jstable@1.6.5/dist/jstable.min.js"
        ]
    else:
        path = source_dir / "template" / "libs"
        links = [
            base64uri(path / "bootstrap" / "bootstrap.min.css"),
            base64uri(path / "jstable.min.css"),
        ]

        import plotly

        plotly_js = Path(plotly.__file__).parent / "package_data" / "plotly.min.js"
        # 轻量
        scripts = [
            base64uri(plotly_js),
            base64uri(path / "bootstrap" / "bootstrap.min.js"),
            base64uri(path / "jstable.min.js"),
        ]
    return links, scripts


def format_summary(summary: dict):
    """
    格式化 summary，返回一个新的字典
    """
    import re
    formatted_summary = summary.copy()

    int_keys = [
        "Number of Reads",
        "Total Gene Detected",
        "Estimated Number of Cells",
        "Cells with Certain Location",
        "Unique Reads in Cells Mapped to Gene",
        "Total Spatial Reads",
        "Spatial Reads with Valid Cellbarcode",
        "Valid Spatial Reads",
        "Valid Spatial UMIs",
        "Valid Spatial Reads in Cells",
        "Valid Spatial UMIs in Cell",
        "Total Spatial Barcodes with Location in Chip",
        "Unique Valid Spatial Barcodes with Location in Chip",
        "Spatial Reads in Cells with Location in Chip",
        *[key for key in formatted_summary.keys() if re.search(r"per Cell$", key)],
    ]
    for k in int_keys:
        if not isinstance(formatted_summary[k], str):
            formatted_summary[k] = f"{formatted_summary[k]:,}"

    pct_keys = [
        "Reads With Valid Barcodes",
        "Sequencing Saturation",
        *[key for key in formatted_summary.keys() if re.search(r"^Q30|^Reads Mapped|^Fraction of", key)],
    ]
    for k in pct_keys:
        if not isinstance(formatted_summary[k], str):
            formatted_summary[k] = f"{formatted_summary[k]:.2%}"

    return formatted_summary


def generate_report(starsolo_output: Path | str, metric_file, outfile, dev=True,cdn = False):
    starsolo_output = Path(starsolo_output)




    template_dir = source_dir / "template"
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("count_template.html")

    _summary = csv2dict(metric_file)
    # 鉴定到的细胞数
    n_cells = _summary["Estimated Number of Cells"]
    summary = format_summary(_summary)

    ###顶部指标项
    keys = [
        "Estimated Number of Cells",
        "Median Gene per Cell",
        "Mean Reads per Cell",
        "Median UMI per Cell",
    ]
    metric_str = "".join(
        [
            f"""
        <div class="col-6 col-md-3">
            <div class="metric-box d-flex">
                <div class="flex-grow-1">
                    <div class="metric-title">{key}</div>
                    <div class="metric-number">{summary[key]}</div>
                </div>
                <div class="metric-icon">
                    <img src={base64uri(template_dir / "bar-icon.png")} alt="bar" width="60px">
                </div>
            </div>
        </div>
        """
            for key in keys
        ]
    )

    # Cells
    # fig = plot_knee(umi_sorted_file=umi_sorted_file, n_cells=n_cells)
    knee_str = fig.to_html(
        full_html=False,
        default_height="380px",
        include_plotlyjs=False,
    )
    if dev:
        keys = [
            "Estimated Number of Cells",
            "Cells with Certain Location",
            "Fraction of Unique Reads in Cells",
            "Fraction of UMIs in Cells",
            "Mean Reads per Cell",
            "Median Reads per Cell",
            "Mean UMI per Cell",
            "Median UMI per Cell",
            "Mean Gene per Cell",
            "Median Gene per Cell",
            "Total Gene Detected",
        ]
        cell_help_msg = {
            "Estimated Number of Cells": "The number of barcodes associated with at least one cell.",
            "Cells with Certain Location": "The number of cell barcodes that can be assigned to a certain location on a chip.",
            "Fraction of Unique Reads in Cells": "The ratio of mapped reads that are assigned to valid cell barcodes.",
            "Fraction of UMIs in Cells": "The proportion of UMIs assigned to valid cells to the total UMIs.",
            "Mean Reads per Cell": "The average number of mapped reads per cell detected.",
            "Median Reads per Cell": "The median number of unique mapped reads per cell detected.",
            "Mean UMI per Cell": "The average number of UMIs per cell detected.",
            "Median UMI per Cell": "The median number of UMIs per cell detected.",
            "Mean Gene per Cell": "The average number of Genes per cell detected.",
            "Median Gene per Cell": "The median number of Genes per cell detected.",
            "Total Gene Detected": "The number of genes with at least one UMI count in cell detected.",
        }
    else:
        keys = [
            "Estimated Number of Cells",
            "Cells with Certain Location",
            "Mean Reads per Cell",
            "Median Reads per Cell",
            "Mean UMI per Cell",
            "Median UMI per Cell",
            "Mean Gene per Cell",
            "Median Gene per Cell",
            "Total Gene Detected",
        ]
        cell_help_msg = {
            "Estimated Number of Cells": "The number of barcodes associated with at least one cell.",
            "Cells with Certain Location": "The number of cell barcodes that can be assigned to a certain location on a chip.",
            "Mean Reads per Cell": "The average number of mapped reads per cell detected.",
            "Median Reads per Cell": "The median number of unique mapped reads per cell detected.",
            "Mean UMI per Cell": "The average number of UMIs per cell detected.",
            "Median UMI per Cell": "The median number of UMIs per cell detected.",
            "Mean Gene per Cell": "The average number of Genes per cell detected.",
            "Median Gene per Cell": "The median number of Genes per cell detected.",
            "Total Gene Detected": "The number of genes with at least one UMI count in cell detected.",
        }

    help_msg = [{"name": k, "value": v} for (k, v) in cell_help_msg.items()]
    content = f"""
    <div class="row">
        <div class="col-6">
            {to_bootstrap_table({key: summary[key] for key in keys})}
        </div>
        <div class="col-6">
            {knee_str}
        </div>
    </div>
    """
    cell_str = report_card(
        title="Cells", help_msg=help_msg, content=content,
    )

    # 序列饱和度
    fig1, fig2 = plot_sequencing_saturation(sequencing_saturation_file)
    sequencing_saturation_str = report_card(
        title="Sequencing_saturation",
        help_msg=[
            {
                "name": "Left",
                "value": "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth (measured in mean reads per cell), up to the observed sequencing depth. Sequencing Saturation is a measure of the observed library complexity, and approaches 1.0 (100%) when all converted mRNA transcripts have been sequenced. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point. The dotted line is drawn at a value reasonably approximating the saturation point.",
            },
            {
                "name": "Right",
                "value": "This plot shows the Median Genes per Cell as a function of downsampled sequencing depth in mean reads per cell, up to the observed sequencing depth. The slope of the curve near the endpoint can be interpreted as an upper bound to the benefit to be gained from increasing the sequencing depth beyond this point.",
            },
        ],
        content=f"""
<div class="row">
    <div class="col-6">
        {
            fig1.to_html(
                full_html=False,
                default_height="380px",
                include_plotlyjs=False,
            )
        }
    </div>
    <div class="col-6">
        {
            fig2.to_html(
                full_html=False,
                default_height="380px",
                include_plotlyjs=False,
            )
        }
    </div>
</div>
""",
    )

    # mapping
    mapping_help_msg = {
        "Reads Mapped to Genome": "The number and fraction of reads mapped to the genome.",
        "Reads Mapped Confidently to Genome": "The number and fraction of reads uniquely mapped to the genome.",
        "Reads Mapped Confidently to Transcriptome": "The number and fraction of reads uniquely mapped to the transcriptome. The read must be consistent with annotated splice junctions.",
        "Reads Mapped Confidently to Exonic": "The number and fraction of reads that uniquely mapped to an exonic region of the genome.",
        "Reads Mapped Confidently to Intronic": "The number and fraction of reads that uniquely mapped to an intronic region of the genome.",
        "Reads Mapped Confidently to Intergenic": "The number and fraction of reads that uniquely mapped to the intergenic region of the genome.",
        "Reads Mapped to antisense": "The number and fraction of reads uniquely mapped to the antisense strand of gene.",
        "Reads Mapped to mito": "Fraction of reads that mapped uniquely to mitochondria region of the genome.",
    }

    keys = [
        "Reads Mapped to Genome",
        "Reads Mapped Confidently to Genome",
        "Reads Mapped Confidently to Transcriptome",
        "Reads Mapped Confidently to Exonic",
        "Reads Mapped Confidently to Intronic",
        "Reads Mapped Confidently to Intergenic",
        "Reads Mapped to antisense",
        "Reads Mapped to mito",
    ]
    content = (
        f'<div class="row"><div class = "col-6">{to_bootstrap_table({key: summary[key] for key in keys[0:4]})}'
        f'</div><div class = "col-6">'
        f"{to_bootstrap_table({key: summary[key] for key in keys[4:8]})}</div></div>"
    )
    mapping_str = report_card(
        title="Mapping",
        help_msg=[{"name": k, "value": v} for (k, v) in mapping_help_msg.items()],
        content=f"{content}",
    )

    # summary
    summary_help_msg = {
        "Number of Reads": "The number of total raw reads.",
        "Reads With Valid Barcodes": "The number and fraction of reads with valid barcodes.",
        "Sequencing Saturation": "Sequencing saturation refers to the percentage of reads originating from a duplicate UMI (in other words, an mRNA molecule that was sequenced more than 1 time).",
        "Q30 Bases in CB+UMI": "The fraction of cell-barcode and UMI bases with Q-scores >= 30.",
        "Q30 Bases in RNA read": "The fraction of RNA read bases with Q-scores >= 30.",
        "Sequencing Saturation Plot": "This plot shows the Sequencing Saturation metric as a function of downsampled sequencing depth, up to the observed sequencing depth. A high saturation indicates we are detecting the vast majority of mRNA molecules in the samples, and thus don't need to sequence the libraries any deeper.",
    }
    keys = [
        "Sample",
        "Reference",
        "Kit",
        "Pipeline version",
        "Number of Reads",
        "Reads With Valid Barcodes",
        "Sequencing Saturation",
        "Q30 Bases in CB+UMI",
        "Q30 Bases in RNA read",
    ]
    summary["Pipeline version"] = _version

    content = (
        f'<div class="row"><div class = "col-6">{to_bootstrap_table({key: summary[key] for key in keys[0:5]})}'
        f'</div><div class = "col-6">'
        f"{to_bootstrap_table({key: summary[key] for key in keys[5:10]})}</div></div>"
    )
    summary_str = report_card(
        title="Summary",
        help_msg=[{"name": k, "value": v} for (k, v) in summary_help_msg.items()],
        content=content,
    )

    # 人鼠混合
    mixed_str = ""
    fig, mixed_ratio = plot_cross(matrix_path)
    if fig is not None:
        content = fig.to_html(
            full_html=False,
            default_height="450px",
            default_width="500px",
            include_plotlyjs=False,
        )
        mixed_str = report_card(
            title="Mixed",
            help_msg=None,
            content=f'{content}',
        )
    # to do: 优化此代码，修改输入文件不是一种很明确的代码
    if mixed_ratio is not None:
        _summary["Mixed Ratio"] = mixed_ratio
        pd.DataFrame(zip(_summary.keys(), _summary.values())).to_csv(
            metric_file,index = False,header=False
        )

    summary_str = f"""
    <div class="row mb-4">
        <div class="col-12">
            <div class="card">
                <div class="card-body">
                    <div class="row">
                        {metric_str}
                    </div>
                </div>
            </div>
        </div>
    </div>
    <div class="row">
        <div class="col-12 mb-4 col-lg-6">
            {cell_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {sequencing_saturation_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {summary_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {mapping_str}
        </div>
        <div class="col-12 mb-4 col-lg-6">
            {mixed_str}
        </div>
    </div>
    """

    adata = sc.read_h5ad(h5ad_file, backed="r")
    # 小提琴图
    fig = plot_violin(adata, samplename=summary["Sample"])
    content = fig.to_html(
        full_html=False,
        default_height="450px",
        default_width="950px",
        include_plotlyjs=False,
    )
    violin_str = report_card(
        title="Distribution",
        help_msg=[
            {
                "name": "Gene counts",
                "value": "The distribution of effective gene counts detected in each cell.",
            },
            {
                "name": "UMI counts",
                "value": "The distribution of UMI counts detected in each cell.",
            },
            {
                "name": "Mito percentage",
                "value": "The distribution of mitochondrial fraction detected in each cell.",
            },
        ],
        content=f'<div class="d-flex justify-content-end">{content}</div>',
        style="width: 1000px",
    )

    #Position and UMAP by UMI count


    fig1, fig2 = plot_umap(adata)

    def _to_html(fig):
        return fig.to_html(
            full_html=False,
            default_height="420px",
            default_width="450px",
            include_plotlyjs=False,
        )

    content = (
        f'<div class="row"><div class = "col-6">{_to_html(fig1)}'
        f'</div><div class = "col-6">'
        f"{_to_html(fig2)}</div></div>"
    )
    cluster_str = report_card(
        title="Cluster",
        help_msg=[
            "The display is limited to a random subset of cells.",
            {
                "name": "left",
                "value": "This plot shows the UMI counts for each cell barcode. Each dot associated with a cell barcode and is colored by the number of total UMI counts. The coordinate axes represents the 2-dimensional embedding produced by the UMAP(Uniform Manifold Approximation and Projection)algorithm.",
            },
            {
                "name": "right",
                "value": "This plot shows the automated clustering result for each cell-barcode by UMAP algorithm. Each dot associated with a cell barcode and is colored according to different cluster.",
            },
        ],
        content=content,
        style="width: 1000px",
    )

    # marker gene
    marker_genes = sc.get.rank_genes_groups_df(
        adata, group=None, key="rank_genes_groups"
    )
    marker_genes.columns = ["cluster", "gene", "scores", "log2FC", "pvals", "pvals_adj"]
    marker_genes = marker_genes[marker_genes["log2FC"] >= 0.25]
    marker_genes = (
        marker_genes.groupby("cluster", observed=True)
        .apply(
            lambda x: x.sort_values(by="pvals_adj", ascending=True).head(30),
            include_groups=False,
        )
        .reset_index(level=0)
    )
    # 格式化
    marker_genes["scores"] = marker_genes["scores"].map(lambda x: f"{x:.4f}")
    marker_genes["log2FC"] = marker_genes["log2FC"].map(lambda x: f"{x:.4f}")
    marker_genes["pvals"] = marker_genes["pvals"].map(lambda x: f"{x:.4f}")
    marker_genes["pvals_adj"] = marker_genes["pvals_adj"].map(lambda x: f"{x:.4f}")

    marker_str = report_card(
        title="Genes",
        help_msg=[
            "The table shows the top 30 differentially expressed genes for each cluster. Here a differential expression test was performed between each cluster and the rest of the sample for each feature. The avg_log2FC is an estimate of the log2 ratio of expression in a cluster to that in all other cells. The p_val is a measure of the statistical significance of the expression difference and the p_val_adj is adjusted p-value, based on bonferroni correction using all features in the dataset."
        ],
        content=tbl(marker_genes, id="table"),
    )
    js = """<script>
let myTable = new JSTable("#table", {
    sortable: true,
    searchable: true,
    perPage: 10,
});
</script>"""

    analysis_str = f"""
<div class="row">
    <div class="col mb-4">
        {violin_str}
        <div class="mt-4">
        {cluster_str}
        </div>
    </div>
    <div class="col mb-4">
        {marker_str}
    </div>
</div>
"""
    analysis_str += js

    links, scripts = get_resource(cdn=cdn)

    content = template.render(
        title="Count",
        logo_uri=base64uri(template_dir / "logo-high.png"),
        favicon=base64uri(template_dir / "favicon.ico"),
        links=links,
        scripts=scripts,
        summary=summary_str,
        analysis=analysis_str,
    )

    with open(outfile, "w", encoding="utf-8") as f:
        f.write(content)

    return None
