# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
app = dash.Dash()

from Bio import Phylo
import pandas as pd
from plotly.grid_objs import Column, Grid
from plotly.offline import download_plotlyjs, init_notebook_mode,  iplot, plot
init_notebook_mode(connected=True)

import plotly.figure_factory as ff
import numpy as np


def get_x_coordinates(tree):
    """Associates to  each clade an x-coord.
       returns dict {clade: x-coord}
    """
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree, and depth is the distance from root
    # to clade

    #  If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords


def get_y_coordinates(tree, dist=1.3):
    """
       returns  dict {clade: y-coord}
       The y-coordinates are  (float) multiple of integers (i*dist below)
       dist depends on the number of tree leafs
    """
    maxheight = tree.count_terminals()  # Counts the number of tree leafs.
    # Rows are defined by the tips/leafs
    ycoords = dict((leaf, maxheight - i * dist) for i, leaf in enumerate(reversed(tree.get_terminals())))

    def calc_row(clade):
        for subclade in clade:
            if subclade not in ycoords:
                calc_row(subclade)
        ycoords[clade] = (ycoords[clade.clades[0]] +
                          ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
    return ycoords


def get_clade_lines(orientation='horizontal', y_curr=0, x_start=0, x_curr=0, y_bot=0, y_top=0,
                    line_color='rgb(25,25,25)', line_width=0.5):
    """define a shape of type 'line', for branch
    """
    branch_line = dict(type='line',
                       layer='below',
                       line=dict(color=line_color,
                                 width=line_width)
                       )
    if orientation == 'horizontal':
        branch_line.update(x0=x_start,
                           y0=y_curr,
                           x1=x_curr,
                           y1=y_curr)
    elif orientation == 'vertical':
        branch_line.update(x0=x_curr,
                           y0=y_bot,
                           x1=x_curr,
                           y1=y_top)
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")

    return branch_line


def draw_clade(clade, x_start, line_shapes, line_color='rgb(15,15,15)', line_width=1):
    """Recursively draw the tree branches, down from the given clade"""

    x_curr = x_coords[clade]
    y_curr = y_coords[clade]

    # Draw a horizontal line from start to here
    branch_line = get_clade_lines(orientation='horizontal', y_curr=y_curr, x_start=x_start, x_curr=x_curr,
                                  line_color=line_color, line_width=line_width)

    line_shapes.append(branch_line)

    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]

        line_shapes.append(get_clade_lines(orientation='vertical', x_curr=x_curr, y_bot=y_bot, y_top=y_top,
                                           line_color=line_color, line_width=line_width))

        # Draw descendants
        for child in clade:
            draw_clade(child, x_curr, line_shapes)



tree = Phylo.read('nextstrain_zika_tree.new', "newick")
x_coords = get_x_coordinates(tree)
y_coords = get_y_coordinates(tree)
line_shapes=[]
draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1)
my_tree_clades=x_coords.keys()
X=[]
Y=[]
text=[]

for cl in my_tree_clades:
    X.append(x_coords[cl])
    Y.append(y_coords[cl])
    text.append(cl.name)

df=pd.read_csv('nextstrain_zika_metadata.csv')
df.columns
nb_genome = len(df)

species=['avian','dengue','ebola','flu','lassa','measles','mumps','zika']

graph_title = "Phylogeny of "+species[0]+" Virus<br>"+ str(nb_genome)+" genomes colored according to region and country"
intermediate_node_color='rgb(100,100,100)'

NA_color={'Cuba': 'rgb(252, 196, 174)',#from cm.Reds color 0.2, ... 0.8
 'Dominican Republic': 'rgb(201, 32, 32)',
 'El Salvador': 'rgb(253, 202, 181)',
 'Guadeloupe': 'rgb(253, 202, 181)',
 'Guatemala': 'rgb(252, 190, 167)',
 'Haiti': 'rgb(252, 145, 114)',
 'Honduras': 'rgb(239, 66, 49)',
 'Jamaica': 'rgb(252, 185, 161)',
 'Martinique': 'rgb(252, 190, 167)',
 'Mexico': 'rgb(247, 109, 82)',
 'Nicaragua': 'rgb(249, 121, 92)',
 'Panama': 'rgb(252, 185, 161)',
 'Puerto Rico': 'rgb(252, 174, 148)',
 'Saint Barthelemy': 'rgb(253, 202, 181)',
 'USA': 'rgb(188, 20, 26)',
 'USVI': 'rgb(206, 36, 34)'}


SAmer_color={'Brazil': 'rgb(21, 127, 59)',# from cm.Greens colors 0.2, 0.4, 0.6, 0.8
 'Colombia': 'rgb(153, 213, 149)',
 'Ecuador': 'rgb(208, 237, 202)',
 'French Guiana': 'rgb(211, 238, 205)',
 'Peru': 'rgb(208, 237, 202)',
 'Suriname': 'rgb(206, 236, 200)',
 'Venezuela': 'rgb(202, 234, 196)'}


SAsia_color={'Singapore': '#0000EE', 'Vietnam': '#1E90FF'}
pl_SAsia=[[0.0, '#1E90FF'], [0.5, '#1E90FF'], [0.5, '#0000EE'], [1.0,'#0000EE' ]]


Oceania_color={'American Samoa': 'rgb(209,95,238)',
 'Fiji': 'rgb(238,130, 238)',
 'French Polynesia': 'rgb(148,0,211)',
 'Tonga': 'rgb(238,130, 238)'}


China_color={'China': 'rgb(255,185,15'}

JapanKorea_color={'Japan': '#fcdd04'}

country = []
region = []
color = [intermediate_node_color] * len(X)

for k, strain in enumerate(df['Strain']):

    i = text.index(strain)

    text[i] = text[i] + '<br>Country: ' + '{:s}'.format(df.loc[k, 'Country']) + '<br>Region: ' + '{:s}'.format(
        df.loc[k, 'Region']) + \
              '<br>Collection date: ' + '{:s}'.format(df.loc[k, 'Date']) + \
              '<br>Journal: ' + '{:s}'.format(df.loc[k, 'Journal']) + '<br>Authors: ' + '{:s}'.format(
        df.loc[k, 'Authors'])
    country.append(df.loc[k, 'Country'])
    region.append(df.loc[k, 'Region'])
    if df.loc[k, 'Region'] == 'North America':
        color[i] = NA_color[df.loc[k, 'Country']]
    elif df.loc[k, 'Region'] == 'South America':
        color[i] = SAmer_color[df.loc[k, 'Country']]
    elif df.loc[k, 'Region'] == 'Southeast Asia':
        color[i] = SAsia_color[df.loc[k, 'Country']]
    elif df.loc[k, 'Region'] == 'Oceania':
        color[i] = Oceania_color[df.loc[k, 'Country']]
    elif df.loc[k, 'Region'] == 'China':
        color[i] = '#fecc00'
    elif df.loc[k, 'Region'] == 'Japan Korea':
        color[i] = '#dc7928'
    else:
        pass


axis=dict(showline=False,
          zeroline=False,
          showgrid=False,
          showticklabels=False,
          title='' #y title
          )

nodes=dict(type='scatter',
           x=X,
           y=Y,
           mode='markers',
           marker=dict(color=color,
                       size=5),
           text=text, #vignet information of each node
           hoverinfo='')

layout=dict(title=graph_title,
            font=dict(family='Balto',size=14),
            width=1000,
            height=3000,
            autosize=False,
            showlegend=False,
            xaxis=dict(showline=True,
                       zeroline=False,
                       showgrid=False,
                       ticklen=4,
                       showticklabels=True,
                       title='branch length'),
            yaxis=axis,
            hovermode='closest',
            shapes=line_shapes,
            plot_bgcolor='rgb(250,250,250)',
            margin=dict(l=10)
           )

fig=dict(data=[nodes], layout=layout)
iplot(fig)


app.layout = html.Div(children=[
    html.H1(children='Phyloegnetic'),
    html.Div(children=[
        dcc.Dropdown(
            options=[
                {'label': 'avian', 'value': 'avian'},
                {'label': 'dengue', 'value': 'dengue'},
                {'label': 'ebola', 'value': 'ebola'},
                {'label': 'flu', 'value': 'flu'},
                {'label': 'lassa', 'value': 'lassa'},
                {'label': 'measles', 'value': 'measles'},
                {'label': 'mumps', 'value': 'mumps'},
                {'label': 'zika', 'value': 'zika'}
            ],
            value='zika',
        )]),

    html.Div(children=[
        dcc.Graph(
            id='example-graph',
            figure=fig
        )])
])


if __name__ == '__main__':
    app.run_server(debug=True, port=5556)

