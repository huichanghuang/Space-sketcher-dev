import pandas as pd
import os
from pathlib import Path
from space_sketcher.tools.utils import csv2dict
from space_sketcher.__init__ import __root_dir__, __version__

def base64uri(filepath):
    """
    将文件转换为 base64 URI 格式   
    参数:
        filepath: 文件路径
    返回:
        base64 URI 字符串，格式为: data:[MIME类型];base64,[base64编码数据]
    """
    import mimetypes
    import base64
    # 猜测文件类型
    mime_type = mimetypes.guess_type(filepath)[0] or 'application/octet-stream'
    # 读取文件并编码为base64
    with open(filepath, "rb") as f:
        encoded_data = base64.b64encode(f.read()).decode('utf-8')
    
    return f"data:{mime_type};base64,{encoded_data}"

def help_icon():
    """
    帮助图标
    """
    icon = """<svg t="1687338260183" class="icon" viewBox="0 0 1024 1024" version="1.1" xmlns="http://www.w3.org/2000/svg" p-id="2757" width="18" height="18" data-spm-anchor-id="a313x.7781069.0.i9"><path d="M512 0C229.23 0 0 229.23 0 512s229.23 512 512 512 512-229.23 512-512S794.77 0 512 0zM512 928c-229.75 0-416-186.25-416-416S282.25 96 512 96s416 186.25 416 416S741.75 928 512 928z" p-id="2758" fill="#707070"></path><path d="M537.64 343.452c47.074 0 83.266-37.528 83.266-78.072 0-32.46-20.832-60.878-62.496-60.878-54.816 0-82.178 44.618-82.178 77.11C475.144 320.132 498.152 343.452 537.64 343.452z" p-id="2759" fill="#707070"></path><path d="M533.162 728.934c-7.648 0-10.914-10.136-3.264-39.55l43.25-166.406c16.386-60.848 10.944-100.398-21.92-100.398-39.456 0-131.458 39.83-211.458 107.798l16.416 27.392c25.246-17.256 67.906-34.762 77.792-34.762 7.648 0 6.56 10.168 0 35.508l-37.746 158.292c-23.008 89.266 1.088 109.538 33.984 109.538 32.864 0 117.808-30.47 195.57-109.632l-18.656-25.34C575.354 716.714 543.05 728.934 533.162 728.934z" p-id="2760" fill="#707070"></path></svg>"""
    return icon

def help_collapse(help_msg: list | str | None = None, _counter=[0]):
    """
    生成折叠的帮助文档，自动生成唯一ID
    参数:
        help_msg: 帮助信息，可以是字符串、列表（支持混合str和dict[name,value]）
        _counter: 内部使用的计数器（无需手动设置）
    
    返回:
        {"icon": 可点击的折叠图标HTML, "text": 折叠内容HTML}
    """
    if help_msg is None:
        return {"icon": "", "text": ""}
    # 生成唯一ID（闭包实现计数器）
    _counter[0] += 1
    collapse_id = f"collapse-{_counter[0]}"
    
    # 处理帮助信息（支持字符串或列表）
    if isinstance(help_msg, list):
        help_msg_str = ""
        for item in help_msg:
            if isinstance(item, str):
                help_msg_str += f"<p>{item}</p>"
            else:
                help_msg_str += f'<h6><strong>{item["name"]}</strong></h6><p>{item["value"]}</p>'
    else:
        help_msg_str = help_msg
    
    # 返回HTML结构
    return {
        "icon": f'<a data-bs-toggle="collapse" href="#{collapse_id}" aria-expanded="false" aria-controls="help">{help_icon()}</a>',
        "text": f'<div class="collapse" id="{collapse_id}">{help_msg_str}</div>'
    }

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
        path = __root_dir__ / "data" / "template" / "libs"
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


def extract_qc_from_cellread_stat(cellread_stat_file):
    """
    从 STARsolo 输出文件 CellReads.stats 提取 QC 信息
    """
    stats = {}
    df = pd.read_csv(cellread_stat_file, sep="\t")
    stats["Reads Mapped to Genome"] = df["genomeU"].sum() + df["genomeM"].sum()
    stats["Reads Mapped Confidently to Genome"] = df["genomeU"].sum()
    stats["Reads Mapped Confidently to Transcriptome"] = df["featureU"].sum()
    stats["Reads Mapped Confidently to Exonic"] = df["exonic"].sum()
    stats["Reads Mapped Confidently to Intronic"] = df["intronic"].sum()
    # 比对到反义链
    stats["Reads Mapped to antisense"] = df["exonicAS"].sum() + df["intronicAS"].sum()
    stats["Reads Mapped Confidently to Intergenic"] = (
        stats["Reads Mapped to Genome"]
        - stats["Reads Mapped Confidently to Exonic"]
        - stats["Reads Mapped Confidently to Intronic"]
        - stats["Reads Mapped to antisense"]
    )
    stats["Reads Mapped to mito"] = df["mito"].sum()
    for k, v in stats.items():
        stats[k] = v / df["cbMatch"].sum()
    return stats

def get_metrics(old_rna_summary_file, new_rna_summary_file, spatial_summary_file, out_summary, args):
    """
    Generate Total summary metrix
    """
    summary = {"Sample": args.sample, "Reference": args.reference,
               "Kit": args.kit}

    ###old rna summary generated by starsolo
    old_rna_summary = csv2dict(old_rna_summary_file)
    old_rna_summary["Total Gene Detected"] = old_rna_summary.pop("Total GeneFull_Ex50pAS Detected")
    old_rna_summary = old_rna_summary["Number of Reads", "Reads With Valid Barcodes", 
                              "Sequencing Saturation", "Q30 Bases in CB+UMI", 
                              "Q30 Bases in RNA read", "Estimated Number of Cells",
                              "Total Gene Detected"]
    ###new rna summary after dbscan filtering
    new_rna_summary = csv2dict(new_rna_summary_file)
    ###spatial library summary generated by spatial_barcode_extraction and assign_coordinate
    spatial_summary = csv2dict(spatial_summary_file)

    ###merge
    summary.update(old_rna_summary)
    summary.update(new_rna_summary)
    summary.update(spatial_summary)

    # Mapping
    cellread_stat_file = (
        os.path.dirname(old_rna_summary_file)/ "CellReads.stats"
    )
    stats = extract_qc_from_cellread_stat(cellread_stat_file)
    summary.update(stats)

    with open(out_summary, "wt") as outf:
        for key, value in summary:
            outf.write(key,",", value)

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