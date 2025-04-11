import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
import math, collections
from plotly.subplots import make_subplots
from scipy.interpolate import make_interp_spline

def cross_plot(crossfile):
    crossdf = pd.read_csv(crossfile, header=0, sep="\t")
    fig = px.scatter(crossdf, x = "column_sum_in_matrix1", y = "column_sum_in_matrix2", color="species")
    fig.update_traces({'opacity': 1.0, 'marker':{'color':'lightskyblue','size':3.5}}, selector={'name': 'GRCh38'})
    fig.update_traces({'opacity': 1.0, 'marker':{'color':'mediumseagreen','size':3.5}}, selector={'name': 'mm10'})
    fig.update_traces({'opacity': 0.5, 'marker':{'color':'lightgray','size':3.5}}, selector={'name': 'Multiplet'})

    fig.update_layout(
        width=600, height=500,
        plot_bgcolor='white', ###set the background color
        margin=dict(l=20, r=20, t=30, b=20), ###set the margin of the plot
        title=dict(text="Cell UMI Counts", font=dict(size=15),x=0.5,y=0.99),
        xaxis_title="mm10", yaxis_title="GRCh38")
    ##update the axis style
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')
    return fig

# barcoderanks plot density
def segment_log_plot(y_data, x_start, x_end):
    log_max_x = np.log(len(y_data))
    log_max_y = np.log(max(y_data))
    segment_len = 0.0
    segment_idx = [x_start]
    for i in range(x_start, x_end):
        last_i = max(x_start, i-1)
        dx = (np.log(i) - np.log(last_i)) / log_max_x
        dy = (np.log(y_data[i]) - np.log(y_data[last_i])) / log_max_y
        segment_len += np.linalg.norm([dx, dy])
        if segment_len >= 0.02 and i > (segment_idx[-1] + 20):
            segment_idx.append(i+1)
            segment_len = 0.0
    if segment_idx[-1] != x_end:
        segment_idx.append(x_end)
    return segment_idx

# barcoderanks plot density
def plot_cmap(density):
    plot_colors =  [
        "#DDDDDD","#D6D9DC","#CFD6DB","#C8D3DA","#C1D0D9",
        "#BACDD9","#B3C9D8","#ACC6D7","#A5C3D6","#9EC0D6",
        "#97BDD5","#90BAD4","#89B6D3","#82B3D3","#7BB0D2",
        "#74ADD1","#6DAAD0","#66A6CF","#5FA3CF","#58A0CE",
        "#539DCC","#4F99CA","#4C95C8","#4992C6","#458EC3",
        "#428AC1","#3F87BF","#3B83BD","#3880BA","#347CB8",
        "#3178B6","#2E75B4","#2A71B1","#276DAF","#236AAD",
        "#2066AB","#1D62A8","#195FA6","#165BA4","#1358A2"]
    levels = len(plot_colors)
    ind = min(levels - 1, int(math.floor(levels * density)))
    return plot_colors[ind]

def plot_barcoderanks_rna(kneefile):

    dataframe_df = pd.read_csv(kneefile, header=0, sep = "\t")
    dataframe_df = dataframe_df.sort_values(by="UMI" , ascending=False)
    dataframe_df = dataframe_df.reset_index(drop=True)
    dataframe_df['New']=dataframe_df.index
    cell_bc = np.array(dataframe_df[dataframe_df['is_cell_barcode'] == 1].index)
    sorted_bc = np.array(dataframe_df.index)
    sorted_counts = np.array(dataframe_df['UMI'])
    total_bc = len(sorted_bc)
    ix1 = dataframe_df.drop_duplicates('is_cell_barcode',keep='first').index[1]-1
    ix2 = dataframe_df.drop_duplicates('is_cell_barcode',keep='last').index[0]
    plot_segments = []
    barcodeSegment = collections.namedtuple(
        'barcodeSegment', 
        ['start', 'end', 'density', 'legend'])

    plot_segments.append(barcodeSegment(start=0, end=ix1, density=1.0, legend=True))
    plot_segments.append(barcodeSegment(start=ix2+1, end=total_bc, density=0.0, legend=True))
    mixed_segments = segment_log_plot(sorted_counts, ix1, ix2)

    for i in range(len(mixed_segments) - 1):
        num_cells = sum([1 for i in range(mixed_segments[i], mixed_segments[i + 1]) if sorted_bc[i] in cell_bc])
        density = float(num_cells)/float(mixed_segments[i + 1]-mixed_segments[i])
        plot_segments.append(barcodeSegment(
                start=mixed_segments[i], end = mixed_segments[i + 1], density=density, legend=False))

    plot_data = []
    for plot_segment in plot_segments:
        start = max(0, plot_segment.start - 1)
        end = plot_segment.end
        selct_count = dataframe_df[start:end]
        dp_first = set(selct_count[selct_count[["UMI"]].duplicated(keep="first")].index)
        dp_last = set(selct_count[selct_count[["UMI"]].duplicated(keep="last")].index)
        dp_inter = dp_first & dp_last
        selct_count=selct_count.drop(list(dp_inter),axis=0)
        x = list(selct_count['New'])
        y = list(selct_count['UMI'])
        name = 'Cell' if plot_segment.density > 0 else 'Background'
        if plot_segment.density > 0:
            n_barcodes = plot_segment.end - plot_segment.start
            n_cells = int(round(plot_segment.density * n_barcodes))
            hover = "{:.0f}% Cell<br>({}/{})".format(100 * plot_segment.density, n_cells, n_barcodes)
        else:
            hover = "NOISE"

        data_dict = {
            "x": x,"y": y,"name": name, "hoverinfo": "text",
            "text": hover,"type": "scattergl","mode": "lines",
            "line": {"width": 3,"color": plot_cmap(plot_segment.density),},
            "showlegend": plot_segment.legend}
        plot_data.append(data_dict)

    plotly_data = [
        go.Scatter(
            x=dat['x'], y=dat['y'], name=dat['name'], mode=dat['mode'], 
            showlegend=dat['showlegend'],
            marker={'color': dat['line']['color']}, 
            line=dat['line'], text=dat['text']
            ) for dat in plot_data]
    
    layout = go.Layout(
        xaxis = dict(
            type="log", gridcolor="whitesmoke", title="Barcode in Rank-descending Order",
            color="black", showline=True, zeroline=True, linewidth=1, fixedrange= True,
            linecolor="black"),
        yaxis = dict(
            type="log", title="UMI counts", gridcolor="whitesmoke",
            linewidth=1, fixedrange= True, color="black", linecolor="black"),
        height= 500, width= 500,
        plot_bgcolor='rgba(0,0,0,0)',hovermode='closest',paper_bgcolor='white',
        legend = dict(
            x=0.75,y=0.99,traceorder="normal",
            font = dict(
                family="Arial",size=12,color="black"),
            bordercolor="Black",borderwidth=0),
        margin = dict(l=0,r=0,b=0,t=0,pad=0.1),
        title=dict(text="Barcode Rank Plot", font=dict(size=15),x=0.5,y=0.99),
        font = dict(size=10))
    
    fig = go.Figure(
        data=plotly_data, layout=layout)

    return fig


def saturation_plot(saturationfile):
    saturantion_df = pd.read_csv(saturationfile, header=0, sep = "\t")
    x, y = saturantion_df["Mean Reads per Cell"], saturantion_df["Sequencing Saturation"]/100
    if len(saturantion_df) > 2:
        xnew = np.linspace(x.min(),x.max(),50)
        ynew = make_interp_spline(x,y)(xnew)
    else:
        xnew = x
        ynew = y
        
    fig = px.line(x=xnew, y=ynew)
    fig.update_traces(line_color='cornflowerblue', line_width=3)
    fig.add_hline(y=0.9, line_width=1, line_dash="dash", line_color="black")
    fig.update_layout(
        width=500, height=500,
        plot_bgcolor='white', ###set the background color
        margin=dict(l=20, r=20, t=30, b=20),
        title=dict(text="Sequencing Saturation", font=dict(size=15),x=0.5,y=0.99),
        showlegend=False) ###set the margin of the plot
    # Update xaxis properties
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke',
        title_text = "Total mapped reads",)
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke',
        title_text = "Sequencing Saturation",
        range = [0,1])

    return fig


def prepare_cluster_data(infile):
    """
    Calculate each cluster's percentage
    sampling data to 20000 cells
    """
    clusterdf = pd.read_csv(infile, header=0, sep="\t")
    ###reset column names
    cols = ["CB", "sample", "UMI counts", "Gene counts", "Mito Percentage", "RNA_snn_res.0.5", "Cluster", "UMAP1", "UMAP2"]
    clusterdf.columns = cols
    clusternum = clusterdf["Cluster"].nunique()
    totalcells = clusterdf.shape[0]
    for i in range(clusternum):
        cluster = clusterdf[clusterdf["Cluster"]==i]
        cluster_percentage = str(round(cluster.shape[0]*100/totalcells, 2))+"%"
        clusterdf.loc[clusterdf["Cluster"]==i, "percentage"] = cluster_percentage
        clusterdf.loc[clusterdf["Cluster"]==i, "Clu"] = cluster_percentage
    ###subset the data to 20000 cells
    if totalcells > 20000:
        clusterdf = clusterdf.sample(n=20000, random_state=1)

    clusterdf["Cluster"] = clusterdf["Cluster"].astype("category")
    ##sort the clusterdf by the Cluster column
    clusterdf = clusterdf.sort_values(by="Cluster")
    return clusterdf

my36colors = ['#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
               '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
               '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
               '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
               '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175']

def cluster_plot(clusterfile):
    clusterdf = prepare_cluster_data(clusterfile)
    fig = make_subplots(rows=1, cols=2, 
                        subplot_titles=("Cells Colored by UMI Counts", "Cells Colored by Cluster"),
                        horizontal_spacing=0.2)
    
    fig.add_trace(px.scatter(clusterdf,
                                x = "UMAP1",
                                y = "UMAP2",
                                color="UMI counts").data[0], 
                                row=1, col=1)    
    
    fig = fig.add_traces(px.scatter(clusterdf,
                        x = "UMAP1",
                        y = "UMAP2",
                        color="Cluster",
                        hover_data=["percentage"],
                        color_discrete_map={i:my36colors[i] for i in range(clusterdf["Cluster"].nunique())},
                        ).data, 
                        rows=1, cols=2)
    
    fig.update_traces(marker_size=3)
    fig.update_layout(
        coloraxis={"colorbar": {"x": 0.4,"len": 1.0,"y": 0.5,'title': 'UMI counts'}})
    fig.update_layout(
        width=1200, height=500,
        plot_bgcolor='white', ###set the background color
        margin=dict(l=20, r=20, t=30, b=20))
    
    ##update the axis style
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')
    
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')
    
    return fig
    
def distribution_violin(clusterfile, samplename):

    clusterdf = pd.read_csv(clusterfile, header=0, sep="\t")
    cols = ["CB", "UMI counts", "Gene counts", "Mito Percentage", "xcoord", "ycoord", "UMAP1", "UMAP2", "Cluster"]
    clusterdf.columns = cols
    clusterdf["sample"] = samplename
    fig = make_subplots(rows=1, cols=3, 
                        subplot_titles=("Gene counts", "UMI counts", "Mito Percentage(%)"),
                        horizontal_spacing=0.1)
    fig.add_trace(go.Figure(data=go.Violin(y=clusterdf["Gene counts"], x=clusterdf["sample"], box_visible=True, line_color='black', meanline_visible=True, fillcolor='#E5D2DD', opacity=0.6, name='Gene counts')).data[0], row=1, col=1)
    fig.add_trace(go.Figure(data=go.Violin(y=clusterdf["UMI counts"], x=clusterdf["sample"], box_visible=True, line_color='black', meanline_visible=True, fillcolor='#53A85F', opacity=0.6, name='UMI counts')).data[0], row=1, col=2)
    fig.add_trace(go.Figure(data=go.Violin(y=clusterdf["Mito Percentage"], x=clusterdf["sample"], box_visible=True, line_color='black', meanline_visible=True, fillcolor='#F1BB72', opacity=0.6, name='Mito Percentage(%)')).data[0], row=1, col=3)
    
    fig.update_layout(
        width=1200, height=500,
        plot_bgcolor='white', ###set the background color
        margin=dict(l=20, r=20, t=30, b=20))
    ##update the axis style
    fig.update_xaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')
    
    fig.update_yaxes(
        mirror=True,
        ticks='outside',
        showline=True,
        linecolor='black',
        gridcolor='whitesmoke')    
    return fig

def spatial_scatter(_infile, _type = "UMI"):
    clusterdf = pd.read_csv(_infile, header=0, sep="\t")
    clusterdf["log_nCount_Spatial"] = np.log(clusterdf["nCount_Spatial"]+1)
    clusternum = clusterdf["Cluster"].nunique()
    totalcells = clusterdf.shape[0]
    for i in range(clusternum):
        cluster = clusterdf[clusterdf["Cluster"]==i]
        cluster_percentage = str(round(cluster.shape[0]*100/totalcells, 2))+"%"
        clusterdf.loc[clusterdf["Cluster"]==i, "percentage"] = cluster_percentage
        clusterdf.loc[clusterdf["Cluster"]==i, "Clu"] = cluster_percentage
    ###subset the data to 20000 cells
    if totalcells > 20000:
        clusterdf = clusterdf.sample(n=20000, random_state=1)

    clusterdf["Cluster"] = clusterdf["Cluster"].astype("category")
    ##sort the clusterdf by the Cluster column
    clusterdf = clusterdf.sort_values(by="Cluster")

    color_map = {}
    for i in range(len(clusterdf["Cluster"].unique())):
        color_map[clusterdf["Cluster"].unique()[i]] = my36colors[i]

    if _type == "Cluster":
        plot1 = px.scatter(clusterdf, 
                           x="xcoord", y="ycoord", 
                           color="Cluster", hover_data=["percentage"],
                           color_discrete_map = color_map,
                           title="Spots colored by Cluster")
        plot1.update_traces(marker_size=3)
        plot1.update_layout(
            autosize=False,
            width=500, height=500,
            plot_bgcolor='white', ###set the background color
            xaxis_tickvals=[0, 2000,4000,6000,8000],
            yaxis_tickvals=[0, 2000,4000,6000,8000],
            yaxis_title=None, xaxis_title=None,
            xaxis=dict(
                    range=[0, 9399],  # 根据实际数据范围调整
                    constrain="domain"),  # 限制轴范围
            yaxis=dict(
                scaleanchor="x",  # y轴比例锚定到x轴
                scaleratio=1,      # 比例1:1
                constrain="domain",  # 限制轴范围
                range=[0, 9399]),          
            margin=dict(l=20, r=20, t=30, b=20))
        ##update the axis style
        plot1.update_xaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
        plot1.update_yaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')

        plot2 = px.scatter(clusterdf,
                            x="UMAP1", y="UMAP2",
                            color="Cluster", hover_data=["percentage"],
                            color_discrete_map = color_map,
                            title="UMAP projection of spots colored by Cluster")
        plot2.update_traces(marker_size=3)
        plot2.update_layout(
            width=500, height=500,
            plot_bgcolor='white', ###set the background color
            margin=dict(l=20, r=20, t=30, b=20),
            xaxis=dict(
                    constrain="domain"),  # 限制轴范围
            yaxis=dict(
                scaleanchor="x",  # y轴比例锚定到x轴
                scaleratio=1,      # 比例1:1
                constrain="domain"),  # 限制轴范围                      
            )

        ##update the axis style
        plot2.update_xaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
        plot2.update_yaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
    elif _type == "UMI":
        plot1 = px.scatter(clusterdf, 
                           x="xcoord", y="ycoord", 
                           color="log_nCount_Spatial", 
                           color_continuous_scale=px.colors.sequential.RdBu_r,
                           title="Spots colored by UMI counts")
        plot1.update_traces(marker_size=3)
    
        plot1.update_layout(
            autosize=False,
            width=600, height=500,
            plot_bgcolor='white', ###set the background color
            xaxis_tickvals=[0, 2000,4000,6000,8000],
            yaxis_tickvals=[0, 2000,4000,6000,8000],
            yaxis_title=None, xaxis_title=None,
            xaxis=dict(
                    range=[0, 9399],  # 根据实际数据范围调整
                    constrain="domain"),  # 限制轴范围
            yaxis=dict(
                scaleanchor="x",  # y轴比例锚定到x轴
                scaleratio=1,      # 比例1:1
                constrain="domain",  # 限制轴范围
                range=[0, 9399]),          
            margin=dict(l=20, r=20, t=30, b=20))
        ##update the axis style
        plot1.update_xaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
        plot1.update_yaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')        
           
        plot2 = px.scatter(clusterdf, 
                           x="UMAP1", y="UMAP2", 
                           color="log_nCount_Spatial", 
                           color_continuous_scale=px.colors.sequential.RdBu_r,
                           title="UMAP projection of spots colored by UMI counts")
        plot2.update_traces(marker_size=3)
        plot2.update_layout(
            width=600, height=500,
            plot_bgcolor='white', ###set the background color
            margin=dict(l=20, r=20, t=30, b=20),
            xaxis=dict(
                    constrain="domain"),  # 限制轴范围
            yaxis=dict(
                scaleanchor="x",  # y轴比例锚定到x轴
                scaleratio=1,      # 比例1:1
                constrain="domain"),  # 限制轴范围                      
            )
        ##update the axis style
        plot2.update_xaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
        plot2.update_yaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
    else:
        print("Please input the correct type: UMI or Cluster")

    return plot1, plot2


def plot_sb_cb_umi_knee(_sb_umi_file, _cb_umi_file):

    sb_umi_df = pd.read_csv(_sb_umi_file, sep = "\t", header = 0)
    # 创建 Log-Log 散点图
    plot1 = px.line(sb_umi_df, x='rank', y='umi_count_sum', log_x=True, log_y=True, markers=True,
                    title='Log-Log Plot of UMI Count Sum by SB',
                    labels={'rank': 'Rank (log scale)', 'umi_count_sum': 'Sum of UMI Count (log scale)'})    
    plot1.update_layout(
        width=600, height=500,
        plot_bgcolor='white', ###set the background color
        yaxis_title=None, xaxis_title=None,
        margin=dict(l=20, r=20, t=30, b=20))
    ##update the axis style
    plot1.update_xaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
    plot1.update_yaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')

    cb_umi_df = pd.read_csv(_cb_umi_file, sep = "\t", header = 0)
    # 创建 Log-Log 散点图
    plot2 = px.line(cb_umi_df, x='rank', y='umi_count_mean', color='cluster', log_x=True, log_y=True, markers=True,
                    # color_discrete_map = {"blue": "single-cluster", "yellow": "multi-cluster", "red": "no-cluster"},
                    color_discrete_sequence = px.colors.qualitative.Pastel2,
                    title='Log-Log Plot of UMI Count Mean by CB',
                    labels={'rank': 'Rank (log scale)', 'umi_count_mean': 'Mean of UMI Count (log scale)'})    
    plot2.update_layout(
        width=600, height=500,
        plot_bgcolor='white', ###set the background color
        yaxis_title=None, xaxis_title=None,
        margin=dict(l=20, r=20, t=30, b=20))
    ##update the axis style
    plot2.update_xaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
    plot2.update_yaxes(mirror=True,ticks='outside',showline=True,linecolor='black', gridcolor='whitesmoke')
    
    return plot1, plot2
