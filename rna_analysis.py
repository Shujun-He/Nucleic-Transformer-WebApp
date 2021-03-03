import pandas as pd
import numpy as np
from bokeh.models import HoverTool
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.plotting import figure
from bokeh.palettes import OrRd3,GnBu3,Spectral5,Spectral4
from bokeh.transform import factor_cmap

import os
os.environ["ARNIEFILE"] = f"arnie.conf"
from arnie.bpps import bpps
from arnie.mfe import mfe

target_columns = ['reactivity', 'deg_Mg_pH10',
                  'deg_pH10', 'deg_Mg_50C', 'deg_50C']

colorlist = ["red", "orange", "yellow", "green", "blue", "black"]


def preprocess_inputs(df, preprocess_columns):
    input = np.transpose(
        np.array(
            df[preprocess_columns]
                .applymap(lambda seq: list(seq))
                .values
                .tolist()
        ),
        (0, 2, 1)
    )
    return input

def flag_correct_target(target_sequence,total_length):
    return list(np.ones(len(target_sequence))) + list(np.zeros(total_length - len(target_sequence)))

def correct_target(target_sequence,total_length):
    return list(target_sequence) + list(np.zeros(total_length - len(target_sequence)))


def get_unseq_data(df, target_columns):
    df["seqpos"] = df.sequence.apply(lambda x: np.arange(1, 1 + len(x)))  # get sequence based position info first.

    flag_columns = []
    for col in target_columns:  # fill unavailable target infos with zeros.
        df[col + '_flag'] = df.apply(lambda x: flag_correct_target(x[col], x["sequence_length"]), axis=1)
        df[col] = df.apply(lambda x: correct_target(x[col], x["sequence_length"]), axis=1)
        flag_columns.append(col + '_flag')

    preprocess_columns = ['seqpos'] + ['sequence', 'structure', 'predicted_loop_type'] + target_columns + flag_columns

    final_inputs = pd.DataFrame()
    for value in df["sequence_length"].unique():
        filter_ = df["sequence_length"] == value
        inputs = preprocess_inputs(df.loc[filter_].copy(), preprocess_columns)  # get processed data.
        inputs = pd.DataFrame(inputs.reshape(-1, len(preprocess_columns)),
                              columns=preprocess_columns)  # turn it into dataframe.
        final_inputs = pd.concat([final_inputs, inputs], axis=0)

    for col in ["seqpos"] + target_columns + flag_columns:
        final_inputs[col] = pd.to_numeric(final_inputs[col])  # make sure some columns to be numeric.

    return final_inputs

def plot_aw(aw,size=500):
    p = figure(title="", width=size, height=size)
    p.image(image=[aw], x=0, y=0, dw=2, dh=2, palette="Spectral11")


    html = file_html(p, CDN, "plot")
    return html

def position_based_plot(unseq_df,target_columns,size=650):
    p = figure(title="", x_axis_label='sequence position', y_axis_label='mean_by_seqposition',width=size, height=size)
    for i, col in enumerate(target_columns):
        g = unseq_df.groupby("seqpos")[col].mean()
        x = g.index
        y = g.values
        p.line(x, y, legend_label=col, line_width=2, line_color=colorlist[i % 6])

    p.add_tools(HoverTool(
        mode="mouse",
        point_policy="follow_mouse"))
    html = file_html(p, CDN, "plot")
    return html

def position_based_plot_single(values,target_columns,size=650):
    p = figure(title="", x_axis_label='sequence position', y_axis_label='mean_by_seqposition',width=int(size*2.5), height=size)
    for i, col in enumerate(target_columns):
        g = values[:,i]
        x = np.arange(len(g))
        y = g
        p.line(x, y, legend_label=col, line_width=2, line_color=colorlist[i % 6])

    p.add_tools(HoverTool(
        mode="mouse",
        point_policy="follow_mouse"))
    html = file_html(p, CDN, "plot")
    return html


def sequence_length_counts(df,seq_column="sequence"):
    seq_list = df[seq_column].apply(lambda x: len(x))
    g = seq_list.value_counts()
    index = list(g.index)
    index = [str(i) for i in index]
    counts = list(g.values)

    p = figure(x_range=index, plot_height=250, title="",
               toolbar_location=None, tools="")

    p.vbar(x=index, top=counts, width=0.23)


    html = file_html(p, CDN, "plot")
    return html



def feature_base_target_means(unseq_df,target_columns,features):

    print("location: feature_base_target_means")
    print(unseq_df.shape)
    # With following part, we make sure that we only take real target values into consideration.
    for col in [col+'_flag' for col in target_columns]:
        unseq_df = unseq_df[unseq_df[col]==1]
    print(unseq_df.shape)

    for f_ in features:
        data = {}
        data[f_] = list(unseq_df[f_].unique())
        for i, col in enumerate(target_columns):
            g = unseq_df.groupby(f_)[col].mean()
            data[col] = list(g.values)

        x = [(seq, target) for seq in data[f_] for target in target_columns]
        counts = []
        for i in range(len(data[list(data.keys())[0]])):
            for key in target_columns:
                counts.append(data[key][i])

        source = ColumnDataSource(data=dict(x=x, counts=counts))

        p = figure(x_range=FactorRange(*x), plot_height=250, title=f"by {f_}",
                   toolbar_location=None, tools="")

        p.vbar(x='x', top='counts', width=0.9, source=source, line_color="white",
               fill_color=factor_cmap('x', palette=Spectral4, factors=target_columns, start=1, end=2))
        p.y_range.start = 0
        p.x_range.range_padding = 0.1
        p.xaxis.major_label_orientation = 1
        p.xgrid.grid_line_color = None
        #p.add_tools(HoverTool(
        #    mode="mouse",
        #    point_policy="follow_mouse"))


    html = file_html(p, CDN, "plot2")
    return html

def arnie_feature(df, feature, package, sequence_column="sequence"):

    s = df[["id"]]
    if feature == "mfe":
        s[feature] = df[sequence_column].apply(lambda x: mfe(x,package=package))
    elif feature == "bpps":
        s[feature] = df[sequence_column].apply(lambda x: bpps(x,package=package))
    elif feature == "unpaired_proba":
        s[feature] = df[sequence_column].apply(lambda x: 1 - np.sum(bpps(x,package=package),axis=0))

    return s

def make_hist_plot(title, hist, edges,size=550):
    p = figure(title=title, tools='', background_fill_color="#fafafa",width=size, height=size)
    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
           fill_color="red", line_color="white", alpha=0.5)

    p.y_range.start = 0
    p.legend.location = "center_right"
    #p.xaxis.axis_label = 'x'
    #p.yaxis.axis_label = 'Pr(x)'
    p.grid.grid_line_color="white"
    html = file_html(p, CDN, "plot")
    return html
