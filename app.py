# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import csv
from Bio import Phylo
import pandas as pd
import numpy as np
from geopy.geocoders import Nominatim
import base64


virus_name = "zika"
species = ['avian', 'dengue', 'ebola', 'flu', 'lassa', 'measles', 'mumps', 'zika']
static_image_route = ''

tree_fig = {}

image_filename = 'img/forum_logo.png'
encoded_image = base64.b64encode(open(image_filename, 'rb').read())


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


def draw_clade(clade, x_start, line_shapes, line_color='rgb(15,15,15)', line_width=1, x_coords=0, y_coords=0):
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
            draw_clade(child, x_curr, line_shapes, x_coords=x_coords, y_coords=y_coords)


def read_treefile(filename):
    tree = Phylo.read(filename, "newick")
    return tree


def read_metadata(filename):
    df = pd.read_csv(filename)
    return df


def create_title(virus, nb_genome):
    graph_title = "Phylogeny of " + virus + " Virus<br>" + str(
        nb_genome) + " genomes colored according to region and country"
    return graph_title


def compute_expensive_data(chemin):
    dir = dir + chemin
    return dir


def create_map():
    df_airports = pd.read_csv(
        'https://raw.githubusercontent.com/plotly/datasets/master/2011_february_us_airport_traffic.csv')
    df_airports.head()

    df_flight_paths = pd.read_csv(
        'https://raw.githubusercontent.com/plotly/datasets/master/2011_february_aa_flight_paths.csv')
    df_flight_paths.head()

    airports = [dict(
        type='scattergeo',
        locationmode='USA-states',
        lon=df_airports['long'],
        lat=df_airports['lat'],
        hoverinfo='text',
        text=df_airports['airport'],
        mode='markers',
        marker=dict(
            size=2,
            color='rgb(255, 0, 0)',
            line=dict(
                width=3,
                color='rgba(68, 68, 68, 0)'
            )
        ))]

    flight_paths = []
    for i in range(len(df_flight_paths)):
        flight_paths.append(
            dict(
                type='scattergeo',
                locationmode='USA-states',
                lon=[df_flight_paths['start_lon'][i], df_flight_paths['end_lon'][i]],
                lat=[df_flight_paths['start_lat'][i], df_flight_paths['end_lat'][i]],
                mode='lines',
                line=dict(
                    width=1,
                    color='red',
                ),
                opacity=float(df_flight_paths['cnt'][i]) / float(df_flight_paths['cnt'].max()),
            )
        )

    layout = dict(
        title='Feb. 2011 American Airline flight paths<br>(Hover for airport names)',
        showlegend=False,
        geo=dict(
            scope='north america',
            projection=dict(type='azimuthal equal area'),
            showland=True,
            landcolor='rgb(243, 243, 243)',
            countrycolor='rgb(204, 204, 204)',
        ),
    )

    fig_map = dict(data=flight_paths + airports, layout=layout)
    return fig_map


def get_lon_lat(city):
    '''
    Example:
    location = geolocator.geocode("Chicago Illinois")
    return:
    Chicago, Cook County, Illinois, United States of America
    location.address    location.altitude   location.latitude   location.longitude  location.point      location.raw
    '''
    geolocator = Nominatim()
    location = geolocator.geocode(city)
    return location.longitude, location.latitude


def get_lon(city):
    '''
    Example:
    location = geolocator.geocode("Chicago Illinois")
    return:
    Chicago, Cook County, Illinois, United States of America
    location.address    location.altitude   location.latitude   location.longitude  location.point      location.raw
    '''
    geolocator = Nominatim()
    location = geolocator.geocode(city)
    return location.longitude


def get_lat(city):
    '''
    Example:
    location = geolocator.geocode("Chicago Illinois")
    return:
    Chicago, Cook County, Illinois, United States of America
    location.address    location.altitude   location.latitude   location.longitude  location.point      location.raw
    '''
    geolocator = Nominatim()
    location = geolocator.geocode(city)
    return location.latitude


def create_map_bubble(metadata_file_stat):
    df = pd.read_csv(metadata_file_stat)
    df.head()

    df['text'] = df['name'] + '<br>Virus proportion ' + (df['pop']).astype(str)
    limits = [(0, 2), (3, 10), (11, 20), (21, 50), (50, 3000)]
    colors = ["rgb(0,116,217)", "rgb(255,65,54)", "rgb(133,20,75)", "rgb(255,133,27)", "lightgrey"]
    cities = []
    scale = 1

    for i in range(len(limits)):
        lim = limits[i]
        df_sub = df[lim[0]:lim[1]]
        city = dict(
            type='scattergeo',
            lon=df_sub['lon'],
            lat=df_sub['lat'],
            text=df_sub['text'],
            marker=dict(
                size=df_sub['pop'] / scale,
                color=colors[i],
                line=dict(width=0.5, color='rgb(40,40,40)'),
                sizemode='area'
            ),
            name='{0} - {1}'.format(lim[0], lim[1]))
        cities.append(city)

    layout = dict(
        title='2018 Dispersion of virus<br>(Click legend to toggle traces)',
        showlegend=True,
        geo=dict(
            showland=True,
            scope='world',
            landcolor='rgb(217, 217, 217)',
            subunitwidth=1,
            countrywidth=1,
            subunitcolor="rgb(255, 255, 255)",
            countrycolor="rgb(255, 255, 255)"
        ),
    )

    fig_map_bubble = dict(data=cities, layout=layout)
    return fig_map_bubble


def create_fig(virus_name, tree_file, metadata_file):
    tree = read_treefile(tree_file)
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)
    line_shapes = []
    draw_clade(tree.root, 0, line_shapes, line_color='rgb(25,25,25)', line_width=1, x_coords=x_coords, y_coords=y_coords)
    my_tree_clades = x_coords.keys()
    X = []
    Y = []
    text = []

    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        text.append(cl.name)

    df = read_metadata(metadata_file)
    data_metadata_stat_csv = df.groupby('Country')['Strain'].count()

    #for index_val, series_val in data_metadata_stat_csv.iteritems():
        #print(index_val, ",", series_val, ",", get_lat(index_val), ",", get_lon(index_val))

    #print(data_metadata_stat_csv)
    #print(type(data_metadata_stat_csv))
    df.columns
    nb_genome = len(df)

    graph_title = create_title(virus_name, nb_genome)
    intermediate_node_color ='rgb(100,100,100)'

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
              'Canada': 'rgb(188, 20, 26)',
              'USVI': 'rgb(206, 36, 34)',
              'Dominica': 'rgb(209, 39, 37)',
              'Trinidad And Tobago': 'rgb(201, 19, 51)',
              'Belize': 'rgb(241, 49, 61)',
              'Grenada': 'rgb(222, 32, 38)',
              'Costa Rica': 'rgb(211, 12, 48)',
              'Bermuda': 'rgb(231, 53, 24)'
              }


    SAmer_color={'Brazil': 'rgb(21, 127, 59)',# from cm.Greens colors 0.2, 0.4, 0.6, 0.8
                 'Colombia': 'rgb(153, 213, 149)',
                 'Ecuador': 'rgb(208, 237, 202)',
                 'French Guiana': 'rgb(211, 238, 205)',
                 'Peru': 'rgb(208, 237, 202)',
                 'Suriname': 'rgb(206, 236, 200)',
                 'Venezuela': 'rgb(202, 234, 196)',
                 'Puerto Rico': 'rgb(201, 235, 199)',
                 'Argentina': 'rgb(203, 225, 185)',
                 'Bolivia': 'rgb(217, 197, 165)',
                 'Paraguay': 'rgb(154, 217, 195)',
                 'Chile': 'rgb(231, 85, 168)',
                 'Aruba': 'rgb(191, 35, 128)'
                 }


    SAsia_color={'Singapore': '#0000EE',
                 'Vietnam': '#1E90FF',
                 'Malaysia': '#1E90AF',
                 'Philippines': '#1E90AE',
                 'Thailand': '#1E90AB',
                 'Myanmar': '#1E90AC',
                 'Cambodia': '#1E90AA',
                 'Indonesia': '#1E90AA',
                 'Brunei': '#1E90BA',
                 'Laos': '#1E90BF'
                 }

    pl_SAsia=[[0.0, '#1E90FF'], [0.5, '#1E90FF'], [0.5, '#0000EE'], [1.0, '#0000EE']]


    Oceania_color={'American Samoa': 'rgb(209,95,238)',
                   'Fiji': 'rgb(238,130, 238)',
                   'French Polynesia': 'rgb(148,0,211)',
                   'Tonga': 'rgb(238,130, 238)',
                   'Australia': 'rgb(233,125, 235)',
                   'Micronesia': 'rgb(231,123, 235)',
                   'New Caledonia': 'rgb(229,119, 233)',
                   'Marshall Islands': 'rgb(227,117, 231)',
                   'Guam': 'rgb(267,137, 251)',
                   'Papua New Guinea': 'rgb(277,187, 291)',
                   'Solomon Islands': 'rgb(22,167, 251)',
                   'Cook Islands': 'rgb(20,187, 211)',
                   'Samoa': 'rgb(50,127, 221)',
                   'Nauru': 'rgb(34, 92, 98)',
                   'Palau': 'rgb(214, 132, 238)',
                   'Vanuatu': 'rgb(252, 42, 128)',
                   'Niue': 'rgb(272, 52, 158)',
                   'New Zealand': 'rgb(242, 71, 133)'
                   }

    China_color={'China': 'rgb(255,185,15'}

    JapanKorea_color={'Japan': '#fcdd04'}

    SubsaharanAfrica_color={'Guinea': 'rgb(209,95,238)',
                            'Liberia': 'rgb(238,130, 238)',
                            'Sierra Leone': 'rgb(148,0,211)',
                            'Cote D Ivoire': 'rgb(145,0,209)',
                            'Angola': 'rgb(143,0,207)',
                            'Seychelles': 'rgb(145,10,217)',
                            'Comoros': 'rgb(141,5,203)',
                            'Madagascar': 'rgb(233,60, 281)',
                            'Eritrea': 'rgb(202, 122, 118)',
                            'Somalia': 'rgb(115,51,222)',
                            'Djibouti': 'rgb(203,57, 211)',
                            'Burkina Faso': 'rgb(141,21,239)',
                            'Ghana': 'rgb(102,57,232)',
                            'Tanzania': 'rgb(217,37, 291)',
                            'Mozambique': 'rgb(213,17, 231)',
                            'Senegal': 'rgb(231,133, 219)',
                            'Togo': 'rgb(121,21,198)',
                            }

    Africa_color={'Sudan': 'rgb(209,95,238)',
                  'Gambia': 'rgb(238,130, 238)',
                  'Nigeria': 'rgb(235,135, 233)',
                  'Mali': 'rgb(235,131, 229)',
                  'Senegal': 'rgb(231,133, 219)',
                  'Cote D Ivoire': 'rgb(145,0,209)',
                  'Burkina Faso': 'rgb(141,21,239)',
                  'Seychelles': 'rgb(145,10,217)',
                  'Somalia': 'rgb(115,51,222)',
                  'Ghana': 'rgb(102,57,232)',
                  'Tanzania': 'rgb(217,37, 291)',
                  'Mozambique': 'rgb(213,17, 231)',
                  'Djibouti': 'rgb(203,57, 211)',
                  'Madagascar': 'rgb(233,60, 281)',
                  'Comoros': 'rgb(141,5,203)',
                  'Togo': 'rgb(121,21,198)',
                  'Angola': 'rgb(212, 92, 138)',
                  'Eritrea': 'rgb(202, 122, 118)',
                  'Guinea': 'rgb(209,95,238)',
                  'Sierra Leone': 'rgb(148,0,211)',
                  'Liberia': 'rgb(238,130, 238)',
                  'Tunisia': 'rgb(228,99, 298)',
                  'Cameroon': 'rgb(207,78, 199)',
                  'South Africa': 'rgb(222,22, 222)',
                  'Congo': 'rgb(231,41, 172)',
                  'Algeria': 'rgb(237,35, 168)',
                  'Morocco': 'rgb(223,27, 165)',
                  'Zambia': 'rgb(218,62, 265)',
                  'Kenya': 'rgb(118,2, 215)',
                  'Uganda': 'rgb(128,35, 265)',
                  'Egypt': 'rgb(143,52, 265)',
                  'Ethiopia': 'rgb(206,36, 265)',
                  'Niger': 'rgb(121,52, 187)',
                  'Mayotte': 'rgb(101,32,165)',
                  'Rwanda': 'rgb(325,25,144)',
                  'Gabon': 'rgb(319,7,197)'
                  }


    Europe_color={'France': 'rgb(209,95,238)',
                  'Germany': 'rgb(238,130, 238)',
                  'Italy': 'rgb(238,130, 238)',
                  'United Kingdom': 'rgb(238,130, 238)',
                  'Netherlands': 'rgb(148,0,211)',
                  'Spain': 'rgb(141,7,221)',
                  'Portugal': 'rgb(139,11,219)',
                  'Ireland': 'rgb(128,15,279)',
                  'Slovakia': 'rgb(121,25,209)',
                  'Romania': 'rgb(171,45,197)',
                  'Sweden': 'rgb(135,96,208)',
                  'Norway': 'rgb(138,56,213)',
                  'Slovenia': 'rgb(138,45,265)',
                  'Denmark': 'rgb(258,25,265)',
                  'Iceland': 'rgb(138,7,185)',
                  'Ukraine': 'rgb(298,65,265)',
                  'Czech Republic': 'rgb(226,96,128)',
                  'Albania': 'rgb(111,20,201)',
                  'Greece': 'rgb(108,63,265)',
                  'Latvia': 'rgb(121,35,299)'
                  }

    country = []
    region = []
    color = [intermediate_node_color] * len(X)
    #print(set(list(df['Region'])))
    #print(set(list(df['Country'])))

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
        elif df.loc[k, 'Region'] == 'Subsaharan Africa':
            color[i] = SubsaharanAfrica_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'Africa':
            color[i] = Africa_color[df.loc[k, 'Country']]
        elif df.loc[k, 'Region'] == 'Europe':
            color[i] = Europe_color[df.loc[k, 'Country']]
        else:
            pass

    axis = dict(showline=False,
              zeroline=False,
              showgrid=False,
              showticklabels=False,
              title='' #y title
              )

    label_legend = set(list(df['Country']))
    nodes = []

    for elt in label_legend:
        node = dict(type='scatter',
                   x=X,
                   y=Y,
                   mode='markers',
                   marker=dict(color=color,
                               size=5),
                   text=text, #vignet information of each node
                   hoverinfo='',
                   name=elt
                   )
        nodes.append(node)

    layout = dict(title=graph_title,
                font=dict(family='Balto', size=14),
                width=1000,
                height=3000,
                autosize=True,
                showlegend=True,
                xaxis=dict(showline=True,
                           zeroline=False,
                           showgrid=True, #To visualize the vertical lines
                           ticklen=4,
                           showticklabels=True,
                           title='branch length'),
                yaxis=axis,
                hovermode='closest',
                shapes=line_shapes,
                plot_bgcolor='rgb(250,250,250)',
                margin=dict(l=10)
               )

    fig = dict(data=nodes, layout=layout)
    return fig

#TO DO validation file and directory exist
def create_paths_file(virus_name, level1="", level2="", level3=""):
    dir = "data/" + virus_name + "/"
    if level1 == "" and level2 == "" and level3 == "":
        tree_file = dir + "nextstrain_" + virus_name + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_metadata.csv"
        stat_file = dir + "stat_nextstrain_" + virus_name + "_metadata.csv"
        return tree_file, metadata_file, stat_file
    elif level2 == "" and level3 == "":
        dir = dir + "/"+level1+"/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_metadata.csv"
        stat_file = dir + "stat_nextstrain_" + virus_name + "_" + level1 + "_metadata.csv"
        return tree_file, metadata_file, stat_file
    elif level3 == "":
        dir = dir + "/" + level1 + "/"+level2+"/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_metadata.csv"
        stat_file = dir + "stat_nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_metadata.csv"
        return tree_file, metadata_file, stat_file
    else:
        dir = dir + "/" + level1 + "/"+level2+"/"+level3+"/"
        tree_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_" + level3 + "_tree.new"
        metadata_file = dir + "nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_" + level3 + "_metadata.csv"
        stat_file = dir + "stat_nextstrain_" + virus_name + "_" + level1 + "_" + level2 + "_" + level3 + "_metadata.csv"
        return tree_file, metadata_file, stat_file



app = dash.Dash()

virus_name = "zika"
species = ['avian', 'dengue', 'ebola', 'flu', 'lassa', 'measles', 'mumps', 'zika']
tree_file, metadata_file, metadata_file_stat = create_paths_file(virus_name, level1="", level2="", level3="")

print(tree_file)
fig = create_fig(virus_name, tree_file, metadata_file)
tree_fig[tree_file] = fig

#fig_map = Virus.create_map()
fig_map_bubble = create_map_bubble(metadata_file_stat)

def serve_layout():
    return html.Div([
        html.Div(
            className="row",
            children=[
                html.Div(
                    className="one columns"
                ),
                html.Div(
                    className="four columns",
                    children=[
                        html.Div(
                            children=html.Div([
                                html.H1(children='Criterion'),
                                html.H1(children=''),
                                html.H6(children='Dataset'),
                                dcc.Dropdown(
                                    id='my-dropdown1',
                                    options=[{'label': species[i], 'value': species[i]} for i in range(len(species))],
                                    value='zika',
                                ),
                                html.Div(id='output-container'),

                                html.Div(id='controls-container_mumps', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown2',
                                        options=[{'label': i, 'value': i} for i in ['global', 'na']],
                                        value='global',
                                    ),
                                ]),

                                html.Div(id='controls-container_dengue', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown3',
                                        options=[{'label': i, 'value': i} for i in ['all', 'denv1', 'denv2', 'denv3', 'denv4']],
                                        value='all',
                                    ),
                                ]),

                                html.Div(id='controls-container_lassa', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown4',
                                        options=[{'label': i, 'value': i} for i in ['s', 'l']],
                                        value='s',
                                    ),
                                ]),

                                html.Div(id='controls-container_avian', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown5',
                                        options=[{'label': i, 'value': i} for i in ['h7n9']],
                                        value='h7n9',
                                    ),
                                    dcc.Dropdown(
                                        id='my-dropdown6',
                                        options=[{'label': i, 'value': i} for i in ['ha', 'mp', 'na', 'ns', 'np', 'pa', 'pb2', 'pb1']],
                                        value='ha',
                                    ),
                                ]),

                                html.Div(id='controls-container_flu', children=[
                                    dcc.Dropdown(
                                        id='my-dropdown7',
                                        options=[{'label': i, 'value': i} for i in ['h3n2', 'h1n1pdm', 'vic', 'yam']],
                                        value='h3n2',
                                    ),
                                    dcc.Dropdown(
                                        id='my-dropdown8',
                                        options=[{'label': i, 'value': i} for i in
                                                 ['ha', 'na']],
                                        value='ha',
                                    ),
                                    dcc.Dropdown(
                                        id='my-dropdown9',
                                        options=[{'label': i, 'value': i} for i in
                                                 ['2y', '3y', '6y', '12y']],
                                        value='3y',
                                    ),
                                ]),

                                html.H1(children=''),
                                html.H1(children=''),
                                html.H6(children='Date Range'),
                                dcc.RangeSlider(
                                    count=1,
                                    min=0,
                                    max=10,
                                    step=0.5,
                                    marks={
                                        0: '0 °F',
                                        3: '3 °F',
                                        5: '5 °F',
                                        7.65: '7.65 °F',
                                        10: '10 °F'
                                    },
                                    value=[3, 7.65]
                                ),
                                html.H1(children=''),
                                html.H1(children=''),
                                html.H1(children=''),
                                html.H6(children='Color by'),
                                dcc.Dropdown(
                                    id='my-dropdown10',
                                    options=[{'label': i, 'value': i} for i in ['genetype', 'country', 'region', 'authors', 'data']],
                                    value='country',
                                ),

                                dcc.Graph(
                                    id='right-mid-graph',
                                    figure=fig_map_bubble
                                )
                            ])
                        )
                    ]
                ),
                html.Div(
                    className="seven columns",
                    children=html.Div([
                        dcc.Graph(
                            id='right-top-graph',
                            figure=fig
                        ),
                        html.Img(id='image',
                                 src='data:image/png;base64,{}'.format(encoded_image)),
                        dcc.Graph(
                            id='right-bottom-graph',
                            figure={
                                'data': [{
                                    'x': [1, 2, 3],
                                    'y': [3, 1, 2],
                                    'type': 'bar'
                                }],
                                'layout': {
                                    'height': 400,
                                    'margin': {'l': 10, 'b': 20, 't': 0, 'r': 0}
                                }
                            }
                        )

                    ])
                ),

            ]
        )
    ])


app.layout = serve_layout()


app.css.append_css({
    'external_url': 'https://codepen.io/chriddyp/pen/bWLwgP.css'
})


@app.callback(
    dash.dependencies.Output('output-container', 'children'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def _update_output(virus_name):
    return 'You have selected "{}" virus'.format(virus_name)


@app.callback(
    dash.dependencies.Output('controls-container_mumps', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def _update_output(virus_name):
    if virus_name == "mumps":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('controls-container_dengue', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def _update_output(virus_name):
    if virus_name == "dengue":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('controls-container_lassa', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def _update_output(virus_name):
    if virus_name == "lassa":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('controls-container_avian', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def _update_output(virus_name):
    if virus_name == "avian":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('controls-container_flu', 'style'),
    [dash.dependencies.Input('my-dropdown1', 'value')])
def _update_output(virus_name):
    if virus_name == "flu":
        return {'display': 'block'}
    else:
        return {'display': 'none'}


@app.callback(
    dash.dependencies.Output('right-top-graph', 'figure'),
    [dash.dependencies.Input('my-dropdown1', 'value'),
     dash.dependencies.Input('my-dropdown2', 'value'),
     dash.dependencies.Input('my-dropdown3', 'value'),
     dash.dependencies.Input('my-dropdown4', 'value'),
     dash.dependencies.Input('my-dropdown5', 'value'), dash.dependencies.Input('my-dropdown6', 'value'),
     dash.dependencies.Input('my-dropdown7', 'value'), dash.dependencies.Input('my-dropdown8', 'value'), dash.dependencies.Input('my-dropdown9', 'value')])
def _update_fig(virus_name, mumps, dengue, lassa, avian_opt1, avian_opt2, flu_opt1, flu_opt2, flu_opt3):
    '''
    I use the underscore in front of this function to suggest that it should not be called. In fact, only changes to the dropdown values should trigger its execution.
    :param virus_name:
    :param mumps:
    :param dengue:
    :param lassa:
    :param avian_opt1:
    :param avian_opt2:
    :param flu_opt1:
    :param flu_opt2:
    :param flu_opt3:
    :return: phylogeny tree
    '''
    if virus_name == "ebola" or virus_name == "zika" or virus_name == "measles":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1="", level2="", level3="")
        return create_fig(virus_name, tree_file_filtred, metadata_file_filtred)
    elif virus_name == "mumps":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=mumps, level2="", level3="")
        return create_fig(virus_name, tree_file_filtred, metadata_file_filtred)
    elif virus_name == "dengue":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat = create_paths_file(virus_name, level1=dengue, level2="", level3="")
        return create_fig(virus_name, tree_file_filtred, metadata_file_filtred)
    elif virus_name == "lassa":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=lassa, level2="", level3="")
        return create_fig(virus_name, tree_file_filtred, metadata_file_filtred)
    elif virus_name == "avian":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=avian_opt1, level2=avian_opt2, level3="")
        return create_fig(virus_name, tree_file, metadata_file)
    elif virus_name == "flu":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=flu_opt1, level2=flu_opt2, level3=flu_opt3)
        return create_fig(virus_name, tree_file_filtred, metadata_file_filtred)


@app.callback(
    dash.dependencies.Output('right-mid-graph', 'figure'),
    [dash.dependencies.Input('my-dropdown1', 'value'),
     dash.dependencies.Input('my-dropdown2', 'value'),
     dash.dependencies.Input('my-dropdown3', 'value'),
     dash.dependencies.Input('my-dropdown4', 'value'),
     dash.dependencies.Input('my-dropdown5', 'value'), dash.dependencies.Input('my-dropdown6', 'value'),
     dash.dependencies.Input('my-dropdown7', 'value'), dash.dependencies.Input('my-dropdown8', 'value'), dash.dependencies.Input('my-dropdown9', 'value')])
def _update_map(virus_name, mumps, dengue, lassa, avian_opt1, avian_opt2, flu_opt1, flu_opt2, flu_opt3):
    if virus_name == "ebola" or virus_name == "zika" or virus_name == "measles":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1="", level2="", level3="")
        return create_map_bubble(metadata_file_stat_filtred)
    elif virus_name == "mumps":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=mumps, level2="", level3="")
        return create_map_bubble(metadata_file_stat_filtred)
    elif virus_name == "dengue":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=dengue, level2="", level3="")
        return create_map_bubble(metadata_file_stat_filtred)
    elif virus_name == "lassa":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=lassa, level2="", level3="")
        return create_map_bubble(metadata_file_stat_filtred)
    elif virus_name == "avian":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=avian_opt1, level2=avian_opt2, level3="")
        return create_map_bubble(metadata_file_stat_filtred)
    elif virus_name == "flu":
        tree_file_filtred, metadata_file_filtred, metadata_file_stat_filtred = create_paths_file(virus_name, level1=flu_opt1, level2=flu_opt2, level3=flu_opt3)
        return create_map_bubble(metadata_file_stat_filtred)


if __name__ == '__main__':
    app.run_server(debug=True, port=5555)
