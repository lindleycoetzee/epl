import pandas as pd
import numpy as np
from dash import Dash, html, dcc, Output, Input, callback
import dash_bootstrap_components as dbc
import dash_ag_grid as dag
import plotly.graph_objects as go

url = "https://en.wikipedia.org/wiki/2024%E2%80%9325_Premier_League"

## get html epl table from wikipedia
dfs = pd.read_html(url)
log = dfs[4]
log = log[log.columns[:-1]]

log_grid = dag.AgGrid(
    rowData = log.to_dict("records"),
    columnSize="responsiveSizeToFit",
    dashGridOptions = {"domLayout": "autoHeight",
                       "rowHeight": 30},
    style = {"height": None},
    columnDefs=[{"field": i} for i in log.columns],
)

## get all epl stats from local file
full_stats = pd.read_csv("eplstats.csv")
full_stats_grid = dag.AgGrid(
    rowData = full_stats.to_dict("records"),
    columnSize="responsiveSizeToFit",
        dashGridOptions = {"domLayout": "autoHeight",
                       "rowHeight": 30},
    columnDefs = [{"field" : i} for i in full_stats.columns],    
    )

## get full match data for all games
match_data = pd.read_csv("PremierLeague.csv")

## create dropdown menus for head to head
team1 = dcc.Dropdown(sorted(match_data["HomeTeam"].unique()), id = "dd_team1", value = "Man United")
team2 = dcc.Dropdown(sorted(match_data["HomeTeam"].unique()), id = "dd_team2", value = "Arsenal")


## create grids for log and full stats
log_table = dbc.Container([
        dbc.Row([log_grid]),
        ])
full_stats_table = dbc.Row([full_stats_grid])


pd.options.mode.chained_assignment = None  # default='warn'

## create layout for head to head tab
h2h_tab = dbc.Container([
    dbc.Row([
        html.Center([html.H1("Selects teams")]),
        dbc.Col([team1,]),
        dbc.Col([team2]),
        ]),
    dbc.Row([
        html.Div(id = "head2head_grid"),
        
        ]),
    html.Br(),
    dbc.Row([
        
        dbc.Col([html.Center(html.H2(html.Div(id = "head2head_team1_name"))),
                 html.Div(id = "head2head_team1_home_data"),
                 html.Div(id = "head2head_team1_away_data"),
                 html.Div(id = "head2head_team1_stats"),
                 ]),
        dbc.Col([html.Center(html.H2(html.Div(id = "head2head_team2_name"))),
                 html.Div(id = "head2head_team2_home_data"),
                 html.Div(id = "head2head_team2_stats"),
                 ]), 
        ]),
    ])


## create app and layout
app = Dash(__name__, suppress_callback_exceptions=True, external_stylesheets=[dbc.themes.JOURNAL],
                meta_tags=[{"name":"viewport",
                            "content": "width=device-width,initial-scale=1.0"}])
app.layout = html.Div([

    dbc.Row([
        html.H1("EPL score predictor")]),

    dbc.Tabs([
            dbc.Tab(log_table, label="Log"),
            dbc.Tab(full_stats_table, label="All Stats"),
            dbc.Tab(h2h_tab, label="Head 2 Head"),
        ]),

    dbc.Row([]),

    ])

@callback(Output("head2head_grid", "children"),
          Output("head2head_team1_name", "children"),
          Output("head2head_team2_name", "children"),
          Output("head2head_team1_stats", "children"),
          Output("head2head_team1_home_data", "children"),
          Output("head2head_team1_away_data", "children"),
          Output("head2head_team2_home_data", "children"),
          Output("head2head_team2_stats", "children"), 
          Input("dd_team1", "value"),
          Input("dd_team2", "value"),)

## create head to head dropdown function
def head2head(t1, t2):
    dff = match_data.copy()
## create head to head table    
    dff = dff[["Season", "Date", "HomeTeam", "AwayTeam", "FullTimeHomeTeamGoals", "FullTimeAwayTeamGoals", "FullTimeResult"]].sort_values("Date", ascending = False)
    dff = dff[((dff["HomeTeam"] == t1) | (dff["AwayTeam"] == t1)) & ((dff["HomeTeam"] == t2) | (dff["AwayTeam"] == t2)) ]
    h2h_grid = dag.AgGrid(
        rowData = dff.to_dict("records"),
        columnSize = "responsiveSizeToFit",
        style={"height": 250, "width": "100%"},
        columnDefs = [{"field" : i} for i in dff.columns])
## calculate how many wins, losses and draws
    dff2 = match_data.copy()
    dff2 = dff2[["Season", "HomeTeam", "AwayTeam", "FullTimeResult"]]
    dff2 = dff2[((dff2["HomeTeam"] == t1) | (dff2["AwayTeam"] == t1)) & ((dff2["HomeTeam"] == t2) | (dff2["AwayTeam"] == t2))]
    pivot_dff = pd.pivot_table(dff2, index = ["HomeTeam", "AwayTeam"], columns = ["FullTimeResult"], values = ["FullTimeResult"], aggfunc = ["count"])
    pivot_dff.columns = pivot_dff.columns.droplevel(0)
    pivot_dff.columns = pivot_dff.columns.droplevel(0)
    pivot_dff = pivot_dff.reset_index()
    pivot_dff = pivot_dff.fillna(0)

## add missing columns for teams that have no wins, losses or draws    
    def add_missing_columns(df, all_columns):
        # Determine the missing columns
        missing_columns = [col for col in all_columns if col not in df.columns]
        
        # Add missing columns with default value of 0
        for col in missing_columns:
##            df[col] = 0
            df = df.copy()
            df[col] = 0

        return df

    all_columns = ['HomeTeam', 'A', 'H', 'D']

## team 1 home wins, losses and draws
    ht = pivot_dff[pivot_dff["HomeTeam"] == t1]   
    ht = add_missing_columns(ht, all_columns)
    del ht["AwayTeam"]

    ht_hg = int(ht.sum(axis = 1))
    ht_hw = int(ht["H"])
    ht_hwp = round(ht_hw / ht_hg, 2)
    ht_hd = int(ht["D"])
    ht_hdp = round(ht_hd / ht_hg, 2)
    ht_hl = int(ht["A"])
    ht_hlp = round(ht_hl / ht_hg, 2)

    text_ht1 = dcc.Markdown(f'''
    #### **Home games : {ht_hg}**
    #### **Home wins : {ht_hw} with a win rate of {round(ht_hwp*100,2)}%**
    #### **Home draws : {ht_hd} with a draw rate of {round(ht_hdp*100,2)}%**
    #### **Home losses : {ht_hl} with a loss rate of {round(ht_hlp*100,2)}%**
    ''')

    win_loss_grid_ht = dag.AgGrid(
        rowData = ht.to_dict("records"),
        columnSize = "responsiveSizeToFit",
        style={"height": 100, "width": "100%"},
        columnDefs = [{"field" : i} for i in ht.columns])

## team 2 home wins, losses and draws
    at = pivot_dff[pivot_dff["HomeTeam"] == t2]
    at = add_missing_columns(at, all_columns)
    del at["AwayTeam"]

    win_loss_grid_at = dag.AgGrid(
        rowData = at.to_dict("records"),
        columnSize = "responsiveSizeToFit",
        style={"height": 100, "width": "100%"},
        columnDefs = [{"field" : i} for i in at.columns])

    at_ag = int(at.sum(axis = 1))
    at_aw = int(at["H"])
    at_awp = round(at_aw / at_ag, 2)
    at_ad = int(at["D"])
    at_adp = round(at_ad / at_ag, 2)
    at_al = int(at["A"])
    at_alp = round(at_al / at_ag, 2)

    text_at = dcc.Markdown(f'''
    #### **Home games : {at_ag}**
    #### **Home wins : {at_aw} with a win rate of {round(at_awp*100,2)}%**
    #### **Home draws : {at_ad} with a draw rate of {round(at_adp*100,2)}%**
    #### **Home losses : {at_al} with a loss rate of {round(at_alp*100,2)}%**
    ''')


## team 1 away wins, losses and draws
    ht_ag = int(at.select_dtypes(include='number').sum(axis=1).sum())
    ht_aw = int(at["A"])
    ht_awp = round(ht_aw / ht_ag, 2)
    ht_ad = int(at["D"])
    ht_adp = round(ht_ad / ht_ag, 2)
    ht_al = int(at["H"])
    ht_alp = round(ht_al / ht_ag, 2)

    text_ht2 = dcc.Markdown(f'''
    #### **Away games : {ht_ag}**
    #### **Away wins : {ht_aw} with a win rate of {round(ht_awp*100,2)}%**
    #### **Away draws : {ht_ad} with a draw rate of {round(ht_adp*100,2)}%**
    #### **Away losses : {ht_al} with a loss rate of {round(ht_alp*100,2)}%**
    ''')
    
        
## team 1 all goal stats
    pt1 = pd.pivot_table(dff, index = ["HomeTeam"], values = ["FullTimeHomeTeamGoals"], aggfunc = ["sum", "count"])
    pt1.columns = pt1.columns.droplevel(0)
    pt1 = pt1.reset_index()
    pt1.columns = ["Team", "HomeGoals", "Matches"]
    
    pt2 = pd.pivot_table(dff, index = ["AwayTeam"], values = ["FullTimeAwayTeamGoals"], aggfunc = ["sum", "count"])
    pt2.columns = pt2.columns.droplevel(0)
    pt2 = pt2.reset_index()
    pt2.columns = ["Team", "AwayGoals", "Matches"]

    dft = pd.merge(pt1, pt2, on = "Team")
    dft.columns = ["Team", "HGoals", "HGames", "AGoals", "AGames"]
    
    dft1 = dft[dft["Team"] == t1 ]
    dft1["TGoals"] = dft1["HGoals"] + dft1["AGoals"]
    dft1["TGames"] = dft1["HGames"] + dft1["AGames"]
    dft1["AveGG"] = round(dft1["TGoals"] / dft1["TGames"], 2)
    dft1["AveHGG"] = round(dft1["HGoals"] / dft1["HGames"], 2)
    dft1["AveAGG"] = round(dft1["AGoals"] / dft1["AGames"], 2)
    dft1_f = dft1[dft1.columns[1:-3]]

    dft1_dag = dag.AgGrid(
        rowData = dft1.to_dict("records"),
        columnSize = "responsiveSizeToFit",
        style={"height": 10, "width": "100%"},
        columnDefs = [{"field" : i} for i in dft1.columns])

    gg_t1 = dft1.iloc[0,7]
    hgg_t1 = dft1.iloc[0,8]
    agg_t1 = dft1.iloc[0,9]

    text_t1 = dcc.Markdown(f'''
    #### **Average goals per game : {gg_t1}**
    #### **Average home goals per game : {hgg_t1}**
    #### **Average away goals per game : {agg_t1}**
    ''')

    dt1 = dag.AgGrid(
        rowData = dft1_f.to_dict("records"),
        columnSize = "responsiveSizeToFit",
        style={"height": 100, "width": "100%"},
        columnDefs = [{"field" : i} for i in dft1_f.columns])
    
## team 2 all goal stats
    dft2 = dft[dft["Team"] == t2 ]
    dft2["TGoals"] = dft2["HGoals"] + dft2["AGoals"]
    dft2["TGames"] = dft2["HGames"] + dft2["AGames"]
    dft2["AveGG"] = round(dft2["TGoals"] / dft2["TGames"], 2)
    dft2["AveHGG"] = round(dft2["HGoals"] / dft2["HGames"], 2)
    dft2["AveAGG"] = round(dft2["AGoals"] / dft2["AGames"], 2)
    dft2_f = dft2[dft2.columns[1:-3]]

    gg_t2 = dft2.iloc[0,7]
    hgg_t2= dft2.iloc[0,8]
    agg_t2= dft2.iloc[0,9]

    text_t2 = dcc.Markdown(f'''
    #### **Average goals per game : {gg_t2}**
    #### **Average home goals per game : {hgg_t2}**
    #### **Average away goals per game : {agg_t2}**
    ''')
  

    dt2 = dag.AgGrid(
        rowData = dft2_f.to_dict("records"),
        columnSize = "responsiveSizeToFit",
        style={"height": 100, "width": "100%"},
        columnDefs = [{"field" : i} for i in dft2_f.columns])

    return h2h_grid, t1, t2, text_t1, text_ht1, text_ht2, text_at, text_t2

if __name__ == "__main__":
    app.run(debug = True)
