from bokeh.io import output_file, show
from bokeh.plotting import figure
from bokeh.palettes import Plasma8
from bokeh.models import ColumnDataSource
import numpy as np

import pandas as pd
import numpy as np
from bokeh.models import HoverTool
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.plotting import figure
from bokeh.palettes import OrRd3,GnBu3,Spectral5,Spectral4
from bokeh.transform import factor_cmap

def plot_top_promoter_kmers(top_kmers,top_kmer_counts,top=10):
    #top=10
    top_indices=np.flip(np.argsort(top_kmer_counts))

    top_kmer_counts=top_kmer_counts[top_indices[:top]]
    top_kmers=top_kmers[top_indices[:top]]

    #output_file("colormapped_bars.html")

    # fruits = ['Apples', 'Pears', 'Nectarines', 'Plums', 'Grapes', 'Strawberries']
    # counts = [5, 3, 4, 2, 4, 6]

    source = ColumnDataSource(data=dict(top_kmers=top_kmers, top_kmer_counts=top_kmer_counts, color=np.flip(Plasma8)))

    p = figure(x_range=top_kmers, title="Top Promoter Kmers",
               toolbar_location=None, tools="", width=680, height=450)

    p.vbar(x='top_kmers', top='top_kmer_counts', width=0.75, color='color', source=source)

    p.xgrid.grid_line_color = None
    p.legend.orientation = "horizontal"
    p.legend.location = "top_center"

    html = file_html(p, CDN, "plot")
    #show(p)
    return html

def plot_top_kmers(top_kmers,top_kmer_counts,top=10):
    #top=10
    top_indices=np.flip(np.argsort(top_kmer_counts))

    top_kmer_counts=top_kmer_counts[top_indices[:top]]
    top_kmers=top_kmers[top_indices[:top]]

    #output_file("colormapped_bars.html")

    # fruits = ['Apples', 'Pears', 'Nectarines', 'Plums', 'Grapes', 'Strawberries']
    # counts = [5, 3, 4, 2, 4, 6]

    source = ColumnDataSource(data=dict(top_kmers=top_kmers, top_kmer_counts=top_kmer_counts, color=np.flip(Plasma8)))

    p = figure(x_range=top_kmers,
               toolbar_location=None, tools="", width=680, height=450)

    p.vbar(x='top_kmers', top='top_kmer_counts', width=0.75, color='color', source=source)

    p.xgrid.grid_line_color = None
    p.legend.orientation = "horizontal"
    p.legend.location = "top_center"

    html = file_html(p, CDN, "plot")
    #show(p)
    return html

def plot_promoter_percent(df):
    #seq_list = df[seq_column].apply(lambda x: len(x))
    #g = df.value_counts()

    categories=['Not Promoter','Promoter',]
    counts=np.zeros(2)
    for s in df.predictions:
        if s == 'not promoter':
            counts[0]+=1
        else:
            counts[1]+=1

    source = ColumnDataSource(data=dict(categories=categories, counts=counts, color=Spectral4))

    p = figure(x_range=categories, title="",
               toolbar_location=None, tools="", width=450, height=450)

    p.vbar(x='categories', top='counts', width=0.5, color='color', source=source)

    html = file_html(p, CDN, "plot")
    return html

def plot_virus_percent(df):
    #seq_list = df[seq_column].apply(lambda x: len(x))
    #g = df.value_counts()

    categories=['Not Viral','Viral',]
    counts=np.zeros(2)
    for s in df.predictions:
        if s == 'not viral':
            counts[0]+=1
        else:
            counts[1]+=1

    source = ColumnDataSource(data=dict(categories=categories, counts=counts, color=Spectral4))

    p = figure(x_range=categories, title="",
               toolbar_location=None, tools="", width=450, height=450)

    p.vbar(x='categories', top='counts', width=0.5, color='color', source=source)

    html = file_html(p, CDN, "plot")
    return html
