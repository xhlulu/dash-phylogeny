import os

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import numpy as np
import pandas as pd
from PIL import ImageColor
import plotly.graph_objs as go
import plotly.figure_factory as ff

import video_engine as rpd
from utils.coco_colors import STANDARD_COLORS

DEBUG = True

app = dash.Dash(__name__)
server = app.server

# Custom Script for Heroku
if 'DYNO' in os.environ:
    app.scripts.config.serve_locally = False
    app.scripts.append_script({
        'external_url': 'https://cdn.rawgit.com/chriddyp/ca0d8f02a1659981a0ea7f013a378bbd/raw/e79f3f789517deec58f41251f7dbb6bee72c44ab/plotly_ga.js'
    })

app.scripts.config.serve_locally = True
app.config['suppress_callback_exceptions'] = True


def load_data(path):
    """Load data about a specific footage (given by the path). It returns a dictionary of useful variables such as
    the dataframe containing all the detection and bounds localization, the number of classes inside that footage,
    the matrix of all the classes in string, the given class with padding, and the root of the number of classes,
    rounded."""

    # Load the dataframe containing all the processed object detections inside the video
    video_info_df = pd.read_csv(path)

    # The list of classes, and the number of classes
    classes_list = video_info_df["class_str"].value_counts().index.tolist()
    n_classes = len(classes_list)

    # Gets the smallest value needed to add to the end of the classes list to get a square matrix
    root_round = np.ceil(np.sqrt(len(classes_list)))
    total_size = root_round ** 2
    padding_value = int(total_size - n_classes)
    classes_padded = np.pad(classes_list, (0, padding_value), mode='constant')

    # The padded matrix containing all the classes inside a matrix
    classes_matrix = np.reshape(classes_padded, (int(root_round), int(root_round)))

    # Flip it for better looks
    classes_matrix = np.flip(classes_matrix, axis=0)

    data_dict = {
        "video_info_df": video_info_df,
        "n_classes": n_classes,
        "classes_matrix": classes_matrix,
        "classes_padded": classes_padded,
        "root_round": root_round
    }

    if DEBUG:
        print(f'{path} loaded.')

    return data_dict


######################################### MAIN APP #########################################
app.layout = html.Div([
    # Banner display
    html.Div([
        html.H2(
            'Object Detection Application',
            id='title'
        ),
        html.Img(
            src="https://s3-us-west-1.amazonaws.com/plotly-tutorials/logo/new-branding/dash-logo-by-plotly-stripe-inverted.png"
        )
    ],
        className="banner"
    ),

    # Body
    html.Div([
        html.Div([
            html.Div([
                html.Div([
                    rpd.my_Player()
                ],
                    id='div-video-player',
                    style={
                        'color': 'rgb(255, 255, 255)',
                        'margin-bottom': '-30px'
                    }
                ),

                html.Div([
                    "Minimum Confidence Threshold:",
                    dcc.Slider(
                        min=20,
                        max=80,
                        marks={i: f'{i}%' for i in range(20, 81, 10)},
                        value=50,
                        id='slider-minimum-confidence-threshold'
                    )
                ],
                    style={'margin': '15px 30px 30px 30px'}  # top right bottom left
                ),

                html.Div([
                    "Footage Selection:",
                    dcc.Dropdown(
                        options=[
                            {'label': 'James Bond', 'value': 'james_bond'},
                            {'label': 'Lion fighting Zebras', 'value': 'zebra'},
                            {'label': 'Man driving expensive car', 'value': 'car_footage'},
                            {'label': 'Drone recording of car festival', 'value': 'car_show_drone'}
                        ],
                        value='james_bond',
                        id="dropdown-footage-selection",
                        clearable=False
                    )
                ],
                    style={'margin': '30px 20px 15px 20px'}  # top right bottom left
                ),

                html.Div([
                    "Video Display Mode:",
                    dcc.Dropdown(
                        options=[
                            {'label': 'Regular Display', 'value': 'regular'},
                            {'label': 'Display with Bounding Boxes', 'value': 'bounding_box'},
                        ],
                        value='regular',
                        id="dropdown-video-display-mode",
                        searchable=False,
                        clearable=False
                    )
                ],
                    style={'margin': '15px 20px 15px 20px'}  # top right bottom left
                ),

                html.Div([
                    "Graph View Mode:",
                    dcc.Dropdown(
                        options=[
                            {'label': 'Visual Mode', 'value': 'visual'},
                            {'label': 'Detection Mode', 'value': 'detection'}
                        ],
                        value='visual',
                        id="dropdown-graph-view-mode",
                        searchable=False,
                        clearable=False
                    )
                ],
                    style={'margin': '15px 20px 15px 20px'}  # top right bottom left
                ),
            ],
                className="six columns",
                style={'margin-bottom': '20px'}
            ),

            html.Div(id="div-visual-mode", className="six columns"),

            html.Div(id="div-detection-mode", className="six columns")
        ],
            className="row"
        )
    ],
        className="container"
    )
])


######################################### DATA LOADING #########################################
@app.server.before_first_request
def load_all_footage():
    global data_dict, url_dict

    # Load the dictionary containing all the variables needed for analysis
    data_dict = {
        'james_bond': load_data("data/video_info.csv"),
        'zebra': load_data("data/Zebra_object_data.csv"),
        'car_show_drone': load_data("data/CarShowDrone_object_data.csv"),
        'car_footage': load_data("data/CarFootage_object_data.csv")
    }

    url_dict = {
        'regular': {
            'james_bond': 'https://www.youtube.com/watch?v=g9S5GndUhko',
            'zebra': 'https://www.youtube.com/watch?v=TVvtD3AVt10',
            'car_show_drone': 'https://www.youtube.com/watch?v=gPtn6hD7o8g',
            'car_footage': 'https://www.youtube.com/watch?v=qX3bDxHuq6I'
        },

        'bounding_box': {
            'james_bond': 'https://www.youtube.com/watch?v=g9S5GndUhko',
            'zebra': 'https://www.youtube.com/watch?v=G2pbZgyWQ5E',
            'car_show_drone': 'https://www.youtube.com/watch?v=9F5FdcVmLOY',
            'car_footage': 'https://www.youtube.com/watch?v=EhnNosq1Lrc'
        }
    }


######################################### FOOTAGE SELECTION #########################################
@app.callback(Output("div-video-player", "children"),
              [Input('dropdown-footage-selection', 'value'),
               Input('dropdown-video-display-mode', 'value')])
def select_footage(footage, display_mode):
    url = url_dict[display_mode][footage]  # Find desired footage

    return [
        rpd.my_Player(
            id='video-display',
            url=url,
            width='100%',
            height='50vh',
            controls=True,
            playing=False,
            seekTo=0,
            volume=1
        )
    ]


######################################### GRAPH VIEW SELECTION #########################################
@app.callback(Output("div-visual-mode", "children"),
              [Input("dropdown-graph-view-mode", "value")])
def update_visual_mode(value):
    if value == "visual":
        return [
            dcc.Interval(
                id="interval-visual-mode",
                interval=800,
                n_intervals=0
            ),

            dcc.Graph(
                style={'height': '55vh'},
                id="heatmap-confidence"
            ),

            dcc.Graph(
                style={'height': '40vh'},
                id="pie-object-count"
            )
        ]

    else:
        return []


@app.callback(Output("div-detection-mode", "children"),
              [Input("dropdown-graph-view-mode", "value")])
def update_detection_mode(value):
    if value == "detection":
        return [
            dcc.Interval(
                id="interval-detection-mode",
                interval=800,
                n_intervals=0
            ),

            dcc.Graph(
                style={'height': '50vh'},
                id="bar-score-graph"
            )
        ]
    else:
        return []


######################################### UPDATING FIGURES #########################################
@app.callback(Output("bar-score-graph", "figure"),
              [Input("interval-detection-mode", "n_intervals")],
              [State("video-display", "currTime"),
               State('dropdown-footage-selection', 'value'),
               State('slider-minimum-confidence-threshold', 'value')])
def update_score_bar(n, current_time, footage, threshold):
    layout = go.Layout(
        title='Detection Score of Most Probable Objects',
        showlegend=False,
        margin=go.Margin(l=70, r=40, t=50, b=30),
        yaxis={'title': 'Score'}
    )

    if current_time is not None:
        current_frame = round(current_time * 23.98)

        if n > 0 and current_frame > 0:
            video_info_df = data_dict[footage]["video_info_df"]

            # Select the subset of the dataset that correspond to the current frame
            frame_df = video_info_df[video_info_df["frame"] == current_frame]

            # Select only the frames above the threshold
            threshold_dec = threshold/100  # Threshold in decimal
            frame_df = frame_df[frame_df["score"] > threshold_dec]

            # Select up to 8 frames with the highest scores
            frame_df = frame_df[:min(8, frame_df.shape[0])]

            # Add count to object names (e.g. person --> person 1, person --> person 2)
            objects = frame_df["class_str"].tolist()
            object_count_dict = {x: 0 for x in set(objects)}  # Keeps count of the objects
            objects_wc = []  # Object renamed with counts
            for object in objects:
                object_count_dict[object] += 1  # Increment count
                objects_wc.append(f"{object} {object_count_dict[object]}")

            # Add text information
            y_text = [f"{round(value*100)}% confidence" for value in frame_df["score"].tolist()]

            # Convert color into rgb
            color_map = lambda class_id: str(ImageColor.getrgb(STANDARD_COLORS[class_id % len(STANDARD_COLORS)]))
            colors = ["rgb" + color_map(class_id) for class_id in frame_df["class"].tolist()]

            bar = go.Bar(
                x=objects_wc,
                y=frame_df["score"].tolist(),
                text=y_text,
                name="Detection Scores",
                hoverinfo="x+text",
                marker=go.Marker(
                    color=colors,
                    line=dict(
                        color='rgb(79, 85, 91)',
                        width=1
                    )
                )
            )

            return go.Figure(data=[bar], layout=layout)

    return go.Figure(data=[go.Bar()], layout=layout)  # Returns empty bar


@app.callback(Output("pie-object-count", "figure"),
              [Input("interval-visual-mode", "n_intervals")],
              [State("video-display", "currTime"),
               State('dropdown-footage-selection', 'value'),
               State('slider-minimum-confidence-threshold', 'value')])
def update_object_count_pie(n, current_time, footage, threshold):
    layout = go.Layout(
        title='Object Count',
        showlegend=True,
        margin=go.Margin(l=50, r=30, t=30, b=40)
    )

    if current_time is not None:
        current_frame = round(current_time * 23.98)

        if n > 0 and current_frame > 0:
            video_info_df = data_dict[footage]["video_info_df"]

            # Select the subset of the dataset that correspond to the current frame
            frame_df = video_info_df[video_info_df["frame"] == current_frame]

            # Select only the frames above the threshold
            threshold_dec = threshold/100  # Threshold in decimal
            frame_df = frame_df[frame_df["score"] > threshold_dec]

            # Get the count of each object class
            class_counts = frame_df["class_str"].value_counts()

            classes = class_counts.index.tolist()  # List of each class
            counts = class_counts.tolist()  # List of each count

            text = [f"{count} detected" for count in counts]

            pie = go.Pie(
                labels=classes,
                values=counts,
                text=text,
                hoverinfo="text+percent",
                textinfo="label+percent"
            )

            return go.Figure(data=[pie], layout=layout)

    return go.Figure(data=[go.Pie()], layout=layout)  # Returns empty pie chart


@app.callback(Output("heatmap-confidence", "figure"),
              [Input("interval-visual-mode", "n_intervals")],
              [State("video-display", "currTime"),
               State('dropdown-footage-selection', 'value'),
               State('slider-minimum-confidence-threshold', 'value')])
def update_heatmap_confidence(n, current_time, footage, threshold):
    layout = go.Layout(
        title="Confidence Level of Object Presence",
        margin=go.Margin(l=20, r=20, t=57, b=30)
    )

    if current_time is not None:
        current_frame = round(current_time * 23.98)

        if n > 0 and current_frame > 0:
            # Load variables from the data dictionary
            video_info_df = data_dict[footage]["video_info_df"]
            classes_padded = data_dict[footage]["classes_padded"]
            root_round = data_dict[footage]["root_round"]
            classes_matrix = data_dict[footage]["classes_matrix"]

            # Select the subset of the dataset that correspond to the current frame
            frame_df = video_info_df[video_info_df["frame"] == current_frame]

            # Select only the frames above the threshold
            threshold_dec = threshold / 100
            frame_df = frame_df[frame_df["score"] > threshold_dec]

            # Remove duplicate, keep the top result
            frame_no_dup = frame_df[["class_str", "score"]].drop_duplicates("class_str")
            frame_no_dup.set_index("class_str", inplace=True)

            # The list of scores
            score_list = []
            for el in classes_padded:
                if el in frame_no_dup.index.values:
                    score_list.append(frame_no_dup.loc[el][0])
                else:
                    score_list.append(0)

            # Generate the score matrix, and flip it for visual
            score_matrix = np.reshape(score_list, (-1, int(root_round)))
            score_matrix = np.flip(score_matrix, axis=0)

            # We set the color scale to white if there's nothing in the frame_no_dup
            if frame_no_dup.shape != (0, 1):
                colorscale = [[0, '#ffffff'], [1, '#f71111']]
                font_colors = ['#3c3636', '#efecee']
            else:
                colorscale = [[0, '#ffffff'], [1, '#ffffff']]
                font_colors = ['#3c3636']

            hover_text = [f"{score * 100:.2f}% confidence" for score in score_list]
            hover_text = np.reshape(hover_text, (-1, int(root_round)))
            hover_text = np.flip(hover_text, axis=0)

            pt = ff.create_annotated_heatmap(
                score_matrix,
                annotation_text=classes_matrix,
                colorscale=colorscale,
                font_colors=font_colors,
                hoverinfo='text',
                text=hover_text
            )

            pt.layout.title = layout.title
            pt.layout.margin = layout.margin

            return pt

    # Returns empty figure
    return go.Figure(data=[go.Pie()], layout=layout)


######################################### CSS #########################################
external_css = [
    "https://cdnjs.cloudflare.com/ajax/libs/normalize/7.0.0/normalize.min.css",  # Normalize the CSS
    "https://fonts.googleapis.com/css?family=Open+Sans|Roboto"  # Fonts
    "https://maxcdn.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css",
    "https://cdn.rawgit.com/xhlulu/dash-object-detection/master/stylesheet.css"
]

for css in external_css:
    app.css.append_css({"external_url": css})


# Running the server
if __name__ == '__main__':
    app.run_server(debug=DEBUG)