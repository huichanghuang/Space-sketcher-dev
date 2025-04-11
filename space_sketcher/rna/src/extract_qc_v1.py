import os
import argparse
import pandas as pd
import base64
from jinja2 import Environment, FileSystemLoader
from Plotly_to_report import plot_barcoderanks_rna, saturation_plot, spatial_scatter, distribution_violin, plot_sb_cb_umi_knee
###细胞数，基因数等信息提取dbscanfilter后重新统计的


def extract_qc_from_starmapping_log(_starlog):
    """
    Extract QC information from STAR mapping log file.(XXXX_Log.final.out)
    """
    mapping_result = {}
    with open(_starlog, 'r') as f:
        for line in f:      
            lineinfo = line.strip().split('|')
            if lineinfo[0].strip() == 'Average input read length':
                mapping_result['average_input_read_length'] = lineinfo[1].strip()
            if lineinfo[0].strip() == 'Number of reads mapped to multiple loci':
                mapping_result['multi_loci_mapped_reads'] = lineinfo[1].strip()
            if lineinfo[0].strip() == '% of reads mapped to multiple loci':
                mapping_result['multi_loci_mapped_rate'] = lineinfo[1].strip()
            if lineinfo[0].strip() == 'Number of reads mapped to too many loci':
                mapping_result['too_many_loci_mapped_reads'] = lineinfo[1].strip()
            if lineinfo[0].strip() == '% of reads mapped to too many loci':
                mapping_result['too_many_loci_mapped_rate'] = lineinfo[1].strip()
            if lineinfo[0].strip() == 'Number of reads unmapped: too many mismatches':
                mapping_result['too_many_mismatches_reads'] = lineinfo[1].strip()
            if lineinfo[0].strip() == '% of reads unmapped: too many mismatches':
                mapping_result['too_many_mismatches_rate'] = lineinfo[1].strip()
            if lineinfo[0].strip() == 'Number of reads unmapped: too short':
                mapping_result['too_short_reads'] = lineinfo[1].strip()
            if lineinfo[0].strip() == '% of reads unmapped: too short':
                mapping_result['too_short_rate'] = lineinfo[1].strip()
            if lineinfo[0].strip() == 'Number of reads unmapped: other':
                mapping_result['other_unmapped_reads'] = lineinfo[1].strip()
            if lineinfo[0].strip() == '% of reads unmapped: other':
                mapping_result['other_unmapped_rate'] = lineinfo[1].strip()
    f.close()

    return mapping_result


def extract_qc_from_cellread_stat(_cellread_stat, outdir):
    """
    Extract QC information from CellReads.stats file.
    """
    cellread_result = {}
    df1 = pd.read_csv(_cellread_stat, sep="\t")
    cellread_result["input_reads"] = df1["cbMatch"].sum()
    ###barcode match rate
    cbMatch = df1[df1["CB"] != "CBnotInPasslist"]
    cellread_result["cbMatchreads"] = cbMatch["cbMatch"].sum()
    cellread_result["cbMatch_rate"] = round(cellread_result["cbMatchreads"]*100/cellread_result["input_reads"], 3) if cellread_result["input_reads"] != 0 else 0
    cellread_result["perfect_cbmatch"] = cbMatch["cbPerfect"].sum()
    cellread_result["perfect_cbmatch_rate"] = round(cellread_result["perfect_cbmatch"]*100/cellread_result["input_reads"], 3) if cellread_result["input_reads"] != 0 else 0
    ###genome mapping
    cellread_result["genome_mapping_reads"] = df1["genomeU"].sum()+df1["genomeM"].sum()
    cellread_result["genome_mapping_rate"] = round(cellread_result["genome_mapping_reads"]*100/cellread_result["input_reads"], 3) if cellread_result["input_reads"] != 0 else 0
    cellread_result["unique_genome_mapping_reads"] = df1["genomeU"].sum()
    cellread_result["unique_genome_mapping_rate"] = round(cellread_result["unique_genome_mapping_reads"]*100/cellread_result["input_reads"], 3) if cellread_result["input_reads"] != 0 else 0
    ###map to transcriptome (confidently mapped reads)
    cellread_result["transcriptome_mapping_reads"] = df1["featureU"].sum() 
    cellread_result["transcriptome_mapping_rate"] = round(cellread_result["transcriptome_mapping_reads"]*100/cellread_result["input_reads"], 3) if cellread_result["input_reads"] != 0 else 0
    ##mito
    cellread_result["mito_mapping_reads"] = df1["mito"].sum()
    cellread_result["mito_mapping_rate"] = round(cellread_result["mito_mapping_reads"]*100/cellread_result["input_reads"], 3) if cellread_result["input_reads"] != 0 else 0
       
    ###exome, intron, antisense mapping rate
    infile3 = os.path.join(os.path.dirname(os.path.dirname(_cellread_stat)) + "/GeneFull_Ex50pAS/CellReads.stats")
    df3 = pd.read_csv(infile3, sep="\t")
    ###exon mapping
    cellread_result["exon_mapping_reads"] = df3["exonic"].sum()
    cellread_result["exon_mapping_rate"] = round(cellread_result["exon_mapping_reads"]*100/cellread_result["genome_mapping_reads"], 3) if cellread_result["genome_mapping_reads"] != 0 else 0
    ###intron mapping
    cellread_result["intron_mapping_reads"] = df3["intronic"].sum()
    cellread_result["intron_mapping_rate"] = round(cellread_result["intron_mapping_reads"]*100/cellread_result["genome_mapping_reads"], 3) if cellread_result["genome_mapping_reads"] != 0 else 0
    ###antisense mapping
    cellread_result["antisense_mapping_reads"] = df3["exonicAS"].sum()+df3["intronicAS"].sum()
    cellread_result["antisense_mapping_rate"] = round(cellread_result["antisense_mapping_reads"]*100/cellread_result["genome_mapping_reads"], 3) if cellread_result["genome_mapping_reads"] != 0 else 0
    ###Intergenic mapping
    INtergenic = cellread_result["genome_mapping_reads"] - cellread_result["exon_mapping_reads"] - cellread_result["intron_mapping_reads"] - cellread_result["antisense_mapping_reads"]
    cellread_result["Intergenic_mapping_reads"] = INtergenic
    cellread_result["Intergenic_mapping_rate"] = round(INtergenic*100/cellread_result["genome_mapping_reads"], 3) if cellread_result["genome_mapping_reads"] != 0 else 0

    outfile = os.path.join(outdir, "cellread_stat_qc.txt")
    with open(outfile, "w") as outf:
        for key, value in cellread_result.items():
            outf.write(f"{key}: {value}\n")
    outf.close()
    
    return cellread_result


def extract_qc_from_star_summary(_starsum):
    """
    Extract QC information from STAR mapping summary file.(Summary.csv)
    """
    summary_result = {}
    with open(_starsum, 'r') as f:
        for line in f:
            lineinfo = line.strip().split(",")
            if lineinfo[0].strip() == "Sequencing Saturation":
                summary_result['sequencing_saturation'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Q30 Bases in CB+UMI":
                summary_result['Q30_bases_in_CB_UMI'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Q30 Bases in RNA read":
                summary_result['Q30_bases_in_RNA_read'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Total GeneFull_Ex50pAS Detected":
                summary_result['total_gene_detected'] = lineinfo[1].strip()
    f.close()

    return summary_result

def extract_qc_from_callcell_summary(_callcellsum):
    callcell_result = {}
    with open(_callcellsum, 'r') as f:
        for line in f:
            lineinfo = line.strip().split(",")
            if lineinfo[0].strip() == "Estimated Number of Cells":
                callcell_result['estimated_number_of_cells'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Unique Reads in Cells Mapped to GeneFull_Ex50pAS":
                callcell_result['unique_reads_in_cells_mapped_to_gene'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Fraction of Unique Reads in Cells":
                callcell_result['fraction_of_unique_reads_in_cells'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Mean Reads per Cell":
                callcell_result['mean_reads_per_cell'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Median Reads per Cell":
                callcell_result['median_reads_per_cell'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "UMIs in Cells":
                callcell_result['UMIs_in_cells'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Mean UMI per Cell":
                callcell_result['mean_UMI_per_cell'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Median UMI per Cell":
                callcell_result['median_UMI_per_cell'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Mean GeneFull_Ex50pAS per Cell":
                callcell_result['mean_gene_per_cell'] = lineinfo[1].strip()
            if lineinfo[0].strip() == "Median GeneFull_Ex50pAS per Cell":
                callcell_result['median_gene_per_cell'] = lineinfo[1].strip()
    f.close()
    return callcell_result


def get_total_qc(_starlog, _starsum, _callcellsum, _cellread_stat, _sample, _outdir):
    """
    Get total QC information from all QC files.
    """
    total_qc = {"sample_name": _sample}
    mapping_result = extract_qc_from_starmapping_log(_starlog)
    total_qc.update(mapping_result)
    cellread_result = extract_qc_from_cellread_stat(_cellread_stat, _outdir)
    total_qc.update(cellread_result)
    starsum_result = extract_qc_from_star_summary(_starsum)
    total_qc.update(starsum_result)
    callcell_result = extract_qc_from_callcell_summary(_callcellsum)
    total_qc.update(callcell_result)

    return total_qc


def export_total_qc_to_excel(_total_qc, _output):
    """
    Export total QC information to xls file.
    """
    ###manager the order of the columns
    _total_qc = pd.DataFrame(_total_qc, index=[0])
    _total_qc = _total_qc[["sample_name", "input_reads", "sequencing_saturation", \
                           "cbMatchreads", "cbMatch_rate", "Q30_bases_in_CB_UMI", "Q30_bases_in_RNA_read", \
                            "estimated_number_of_cells","mean_reads_per_cell","median_reads_per_cell", \
                            "mean_gene_per_cell", "median_gene_per_cell", "total_gene_detected",\
                            "UMIs_in_cells", "mean_UMI_per_cell","median_UMI_per_cell", \
                            "genome_mapping_reads", "genome_mapping_rate", "unique_genome_mapping_reads", "unique_genome_mapping_rate", \
                            "transcriptome_mapping_reads", "transcriptome_mapping_rate", \
                            "exon_mapping_reads", "exon_mapping_rate", "intron_mapping_reads", "intron_mapping_rate", \
                            "antisense_mapping_reads", "antisense_mapping_rate", "Intergenic_mapping_reads", "Intergenic_mapping_rate", \
                            "mito_mapping_reads", "mito_mapping_rate", \
                            "unique_reads_in_cells_mapped_to_gene", "fraction_of_unique_reads_in_cells", \
                            "multi_loci_mapped_reads", "multi_loci_mapped_rate", "too_many_loci_mapped_reads", "too_many_loci_mapped_rate",\
                            "too_many_mismatches_reads", "too_many_mismatches_rate", "too_short_reads", "too_short_rate",\
                            "other_unmapped_reads", "other_unmapped_rate", ]]

    _total_qc.to_csv(_output, index=False, sep="\t")

    return


def prepare_dataframes_html(diffgenesfile):

    diffgenes = pd.read_csv(diffgenesfile, sep="\t", header=0)
    colnames = ["gene", "cluster", "p_val", "avg_log2FC", "pct1", "pct2", "p_val_adj"]
    ##change the column name
    diffgenes.columns = colnames
    
    dataframes_html ='\t<div class="table-demo ">\n' + \
                    '\t\t<table id="table"></table>\n' + \
                    '\t</div>\n' + \
                    '<script>\n' + \
                    '//设置需要显示的列\n' + \
                    'var columns = [\n'
    
    for col in colnames:
        dataframes_html += f'\t\t{{field: "{col}", title: "{col}"}},\n'
    
    dataframes_html += '\t];\n'
    
    dataframes_html += '//bootstrap table初始化数据\n' + \
                    '$("#table").bootstrapTable({\n' + \
                    '\tcolumns: columns,\n' + \
                    '\tdata: getData(),\n' + \
                    '\tclasses: "table table-bordered table-striped table-sm", //设置表格样式\n' + \
                    '\theight:400,\n' + \
                    '\t//******前端分页设置****\n' + \
                    '\tpagination: true,\n' + \
                    '\tpageSize: 10,\n' + \
                    '\tpageList: [10, 25, 50, 100],\n' + \
                    '\tpageNumber:2,\n' + \
                    '\tsearch: true,\n' + \
                    '\tpaginationHAlign:"left",\n' + \
                    '\tpaginationDetailHAlign:"right"\n' + \
                    '});\n'
    
    dataframes_html += 'function getData() {\n' + \
                    '\tvar data = [];\n'
    
    for index, row in diffgenes.iterrows():
        dataframes_html += '\tdata.push({\n'
        for col in colnames:
            dataframes_html += f'\t\t{col}: "{row[col]}",\n'
        dataframes_html += '\t});\n'
    dataframes_html += '\treturn data;\n' + \
                    '}\n' + \
                    '</script>\n'

    return dataframes_html

def prepare_cb_sb_umi_plot(_sb_umi_file, _cb_umi_file):
    sb_umi, cb_umi = plot_sb_cb_umi_knee(_sb_umi_file, _cb_umi_file)
    plot_html = f"""
    <h4 class="fw-bold">Spatial Barcode UMI Distribution
        <button class="btn-sm btn-link btn-outline-light active" type="submit" data-bs-toggle="collapse" data-bs-target="#SB_UMI_info">
            <i class="bi bi-info-circle" style="font-size: 1rem;"></i>
        </button>
    </h4>
    <div class="collapse" id="SB_UMI_info">
        <div class="card card-body">
            <p>The display is limited to a random subset of Spots.</p>
            <p><span class="fw-bold" style="white-space: pre;">left:    </span>Total spatial UMI count for each spatial barcode plotted against its rank.</p>
            <p><span class="fw-bold" style="white-space: pre;">right:    </span>Mean spatial UMI count for each cell barcode plotted against its rank.
            Each cell barcode is colored according to its cluster type based on the DBSCAN clustering result.
            </p>
        </div>
    </div>
    <div class="col text-center p-1">
      {sb_umi.to_html(full_html=False)}
    </div>
    <div class="col text-center p-1">
      {cb_umi.to_html(full_html=False)}
    </div>
    """
    return plot_html


def prepare_plotly_html(_total_qc, _indir, _sample, dev=True):

    kneefile = os.path.join(_indir, f"{_sample}.barcode_rank.txt")
    if os.path.exists(kneefile):
        kneefig = plot_barcoderanks_rna(kneefile)
        _total_qc["knee_plot_bc"] = kneefig.to_html(full_html=False)
    else:
        print(f"{kneefile} not exits, please check the results.")

    saturationfile = os.path.join(_indir, f"{_sample}.saturation.xls")
    if os.path.exists(saturationfile):
        saturationfig = saturation_plot(saturationfile)
        _total_qc["saturation_plot_bc"] = saturationfig.to_html(full_html=False)
    else:
        print(f"{saturationfile} not exits, please check the results.")

    if int(_total_qc["estimated_number_of_cells"]) >= 100:
        clusterfile = os.path.join(_indir, f"{_sample}.UMAPpos.txt")
        if os.path.exists(clusterfile):

            distributionfig = distribution_violin(clusterfile, _sample)
            _total_qc["distribution_plot_bc"] = distributionfig.to_html(full_html=False)

            dim_umi,umap_umi = spatial_scatter(clusterfile, "UMI")
            _total_qc["dim_umi"] = dim_umi.to_html(full_html=False)
            _total_qc["umap_umi"] = umap_umi.to_html(full_html=False)

            dim_cluster,umap_cluster = spatial_scatter(clusterfile, "Cluster")
            _total_qc["dim_cluster"] = dim_cluster.to_html(full_html=False)
            _total_qc["umap_cluster"] = umap_cluster.to_html(full_html=False)
        
        diffgenes = os.path.join(_indir, f"{_sample}.topDEG.txt")
        if os.path.exists(diffgenes):
            _total_qc["DEGtable"] = prepare_dataframes_html(diffgenes)

    sb_umi_file = os.path.join(_indir, f"{_sample}.spatial_umi_knee_plot.txt")
    cb_umi_file = os.path.join(_indir, f"{_sample}.cb_umi_mean_knee_plot.txt")
    if dev and os.path.exists(cb_umi_file) and os.path.exists(sb_umi_file):
        _total_qc["sb_cb_umi"] = prepare_cb_sb_umi_plot(sb_umi_file, cb_umi_file)

    return _total_qc



def main():
    parser = argparse.ArgumentParser(description="Extract QC information from QC files.")
    parser.add_argument("-starlog", help="STAR mapping log file (XXXX_Log.final.out)", required=True)
    parser.add_argument("-starsum", help="STAR mapping summary file (Summary.csv) ", required=True)
    parser.add_argument("-cellread_stat", help="CellReads.stats file", required=True)
    parser.add_argument("-callcellsum", help="CallCell summary file", required=False)
    parser.add_argument("-output", help="Output excel file", required=True)
    parser.add_argument("-template", help="QC report template", required=True)
    parser.add_argument("-p", "--pipeline", help="Pipeline", required=True)
    parser.add_argument("-r", "--reference", help="Reference", required=True)
    parser.add_argument("-k", "--kit", help="Kit", required=True)
    parser.add_argument("-d", "--dev", help="development mode, defaulte True", default=True)
    args = parser.parse_args()

    indir = os.path.dirname(args.output)
    sample = os.path.basename(args.output).split(".")[0]
    htmloutput = os.path.join(indir, f"{sample}.QC_report.html")
    total_qc = get_total_qc(args.starlog, args.starsum, args.callcellsum, args.cellread_stat, sample, indir)
    export_total_qc_to_excel(total_qc, args.output)

    total_qc = prepare_plotly_html(total_qc, indir, sample, args.dev)
    total_qc["reference"] = args.reference
    total_qc["Pipeline"] = args.pipeline
    total_qc["Kit"] = args.kit

    templatedir = os.path.dirname(args.template)
    tempfile = os.path.basename(args.template)
    env = Environment(loader=FileSystemLoader(templatedir))
    template = env.get_template(tempfile)

    with open(htmloutput, 'w+', encoding='utf-8') as f:
        f.write(template.render(total_qc))
    f.close()

    return

if __name__ == "__main__":
    main()