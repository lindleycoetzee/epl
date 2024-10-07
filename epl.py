import pandas as pd
import numpy as np
from dash import Dash, html, dcc, Output, Input, callback
import dash_bootstrap_components as dbc
import dash_ag_grid as dag
import plotly.graph_objects as go
from functools import reduce

## give a probabilty score
## error handling
## styling

url = "https://en.wikipedia.org/wiki/2024%E2%80%9325_Premier_League"

## get html epl table from wikipedia
dfs = pd.read_html(url)
log = dfs[4]
log = log[log.columns[:-1]]

columnDefs = [
    {"field": "Pos", "width": 40},
    {"field": "Team", "width": 100},
    {"field": "Pld", "width": 40},
    {"field": "W", "width": 40},
    {"field": "D", "width": 40},
    {"field": "L", "width": 40},
    {"field": "GF", "width": 40},
    {"field": "GA", "width": 40},
    {"field": "GD", "width": 40},
    {"field": "Pts", "width": 40},
]

log_grid = dag.AgGrid(
    rowData=log.to_dict("records"),
    columnSize="responsiveSizeToFit",
    dashGridOptions={"domLayout": "autoHeight", "rowHeight": 30},
    style={"height": None},
    ##    columnDefs=[{"field": i} for i in log.columns],
    columnDefs=columnDefs,
)


## get all epl stats from local file
full_stats = pd.read_csv("eplstats.csv")
full_stats_grid = dag.AgGrid(
    rowData=full_stats.to_dict("records"),
    columnSize="responsiveSizeToFit",
    style={"height": 800, "width": "100%"},
    columnDefs=[{"field": i} for i in full_stats.columns],
)

## get all epl goals from local file
goals_data = pd.read_csv("eplgoals.csv")
goals_data_grid = dag.AgGrid(
    rowData=goals_data.to_dict("records"),
    columnSize="responsiveSizeToFit",
    style={"height": 800, "width": "100%"},
    columnDefs=[{"field": i} for i in goals_data.columns],
)

## get full match data for all games and goals
match_data = pd.read_csv("PremierLeague(2).csv")

## create dropdown menus for head to head
team1 = dcc.Dropdown(
    sorted(match_data["HomeTeam"].unique()), id="dd_team1", value="Man United"
)
team2 = dcc.Dropdown(
    sorted(match_data["HomeTeam"].unique()), id="dd_team2", value="Arsenal"
)

## create dropdown menus for season(s)
seasons = [
    "1993-1994",
    "1994-1995",
    "1995-1996",
    "1996-1997",
    "1997-1998",
    "1998-1999",
    "1999-2000",
    "2000-2001",
    "2001-2002",
    "2002-2003",
    "2003-2004",
    "2004-2005",
    "2005-2006",
    "2006-2007",
    "2007-2008",
    "2008-2009",
    "2009-2010",
    "2010-2011",
    "2011-2012",
    "2012-2013",
    "2013-2014",
    "2014-2015",
    "2015-2016",
    "2016-2017",
    "2017-2018",
    "2018-2019",
    "2019-2020",
    "2020-2021",
    "2021-2022",
    "2022-2023",
    "2023-2024",
    "2024-2025",
]
season = dcc.Dropdown(
    match_data["Season"].unique(), id="season", value=seasons, multi=True
)


## create grids for log, full stats and goals
log_table = dbc.Container(
    [
        dbc.Row([log_grid]),
    ]
)
full_stats_table = dbc.Row([full_stats_grid])

goals_table = dbc.Container(
    [
        dbc.Row([goals_data_grid]),
    ]
)


pd.options.mode.chained_assignment = None  # default='warn'

## create layout for head to head tab
h2h_tab = html.Div(
    [
        dbc.Container(
            [
                dbc.Row(
                    [
                        html.Center([html.H4("Select teams")]),
                        dbc.Col(
                            [
                                team1,
                            ]
                        ),
                        dbc.Col([team2]),
                    ]
                ),
                dbc.Row(
                    [
                        html.Center([html.H4("Select season(s)")]),
                        season,
                    ]
                ),
                dbc.Row(
                    [
                        html.Div(id="head2head_grid"),
                    ]
                ),
                html.Br(),
            ]
        ),
        dbc.Container(
            [
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                html.Center(html.H5("Versus opponent")),
                                html.Div(id="head2head_team1_home_data"),
                                html.Div(id="head2head_team1_away_data"),
                                html.Div(id="head2head_team1_total_data"),
                                html.Div(id="head2head_team1_stats"),
                            ]
                        ),
                        dbc.Col(
                            [
                                html.Center(html.H5("Overall")),
                                html.Div(id="ovr_hg_team1_total_data"),
                                html.Div(id="ovr_ag_team1_total_data"),
                                html.Div(id="ovr_team1_total_data"),
                                html.Div(id="ovr_team1_goal_data"),
                            ]
                        ),
                        dbc.Col(
                            [
                                html.Center(html.H5("Versus opponent")),
                                html.Div(id="head2head_team2_home_data"),
                                html.Div(id="head2head_team2_away_data"),
                                html.Div(id="head2head_team2_total_data"),
                                html.Div(id="head2head_team2_stats"),
                            ]
                        ),
                        dbc.Col(
                            [
                                html.Center(html.H5("Overall")),
                                html.Div(id="ovr_hg_team2_total_data"),
                                html.Div(id="ovr_ag_team2_total_data"),
                                html.Div(id="ovr_team2_total_data"),
                                html.Div(id="ovr_team2_goal_data"),
                            ]
                        ),
                    ]
                ),
            ]
        ),
    ]
)


## create app and layout
app = Dash(
    __name__,
    suppress_callback_exceptions=True,
    external_stylesheets=[dbc.themes.UNITED],
    meta_tags=[{"name": "viewport", "content": "width=device-width,initial-scale=1.0"}],
)
app.layout = html.Div(
    [
        dbc.Row([html.H1("EPL score predictor")]),
        dbc.Tabs(
            [
                dbc.Tab(log_table, label="Log"),
                dbc.Tab(full_stats_table, label="All Win Loss Draw Stats"),
                dbc.Tab(goals_table, label="Goals Stats"),
                dbc.Tab(h2h_tab, label="Head 2 Head"),
            ]
        ),
        dbc.Row([]),
    ]
)


@callback(
    Output("head2head_grid", "children"),
    Output("head2head_team1_stats", "children"),
    Output("head2head_team1_home_data", "children"),
    Output("head2head_team1_away_data", "children"),
    Output("head2head_team1_total_data", "children"),
    Output("ovr_hg_team1_total_data", "children"),
    Output("ovr_ag_team1_total_data", "children"),
    Output("ovr_team1_total_data", "children"),
    Output("ovr_team1_goal_data", "children"),
    Output("ovr_hg_team2_total_data", "children"),
    Output("ovr_ag_team2_total_data", "children"),
    Output("ovr_team2_total_data", "children"),
    Output("head2head_team2_home_data", "children"),
    Output("head2head_team2_away_data", "children"),
    Output("head2head_team2_total_data", "children"),
    Output("head2head_team2_stats", "children"),
    Output("ovr_team2_goal_data", "children"),
    Input("dd_team1", "value"),
    Input("dd_team2", "value"),
    Input("season", "value"),
)

## create head to head dropdown function
def head2head(t1, t2, season):
    dff = match_data.copy()
    ## create head to head table
    dff = dff[
        [
            "Season",
            "Date",
            "HomeTeam",
            "AwayTeam",
            "FullTimeHomeTeamGoals",
            "FullTimeAwayTeamGoals",
            "FullTimeResult",
        ]
    ].sort_values("Date", ascending=False)
    dff = dff[
        ((dff["HomeTeam"] == t1) | (dff["AwayTeam"] == t1))
        & ((dff["HomeTeam"] == t2) | (dff["AwayTeam"] == t2))
    ]
    dff = dff[dff["Season"].isin(season)]
    dff = dff.sort_values(by="Season", ascending=False)

    h2h_grid = dag.AgGrid(
        rowData=dff.to_dict("records"),
        columnSize="responsiveSizeToFit",
        style={"height": 250, "width": "100%"},
        columnDefs=[{"field": i} for i in dff.columns],
    )
    ## calculate how many wins, losses and draws
    dff2 = match_data.copy()
    dff2 = dff2[["Season", "HomeTeam", "AwayTeam", "FullTimeResult"]]
    dff2 = dff2[
        ((dff2["HomeTeam"] == t1) | (dff2["AwayTeam"] == t1))
        & ((dff2["HomeTeam"] == t2) | (dff2["AwayTeam"] == t2))
    ]
    dff2 = dff2[dff2["Season"].isin(season)]
    pivot_dff = pd.pivot_table(
        dff2,
        index=["HomeTeam", "AwayTeam"],
        columns=["FullTimeResult"],
        values=["FullTimeResult"],
        aggfunc=["count"],
    )
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
            df = df.copy()
            df[col] = 0

        return df

    all_columns = ["HomeTeam", "A", "H", "D"]

    ## team 1 home wins, losses and draws
    ht = pivot_dff[pivot_dff["HomeTeam"] == t1]
    ht = add_missing_columns(ht, all_columns)
    del ht["AwayTeam"]

    ht_hg = int(ht.sum(axis=1))
    ht_hw = int(ht["H"])
    ht_hwp = round(ht_hw / ht_hg * 100, 2)
    ht_hd = int(ht["D"])
    ht_hdp = round(ht_hd / ht_hg * 100, 2)
    ht_hl = int(ht["A"])
    ht_hlp = round((ht_hl / ht_hg * 100), 2)

    text_ht1 = dcc.Markdown(
        f"""
    ###### **Home games : {ht_hg}**
    ###### **Home wins : {ht_hw} | Win rate : {ht_hwp}%**
    ###### **Home draws : {ht_hd} | Draw rate : {ht_hdp}%**
    ###### **Home losses : {ht_hl} | Loss rate : {ht_hlp}%**
    """
    )

    win_loss_grid_ht = dag.AgGrid(
        rowData=ht.to_dict("records"),
        columnSize="responsiveSizeToFit",
        style={"height": 100, "width": "100%"},
        columnDefs=[{"field": i} for i in ht.columns],
    )

    ## team 2 home wins, losses and draws
    at = pivot_dff[pivot_dff["HomeTeam"] == t2]
    at = add_missing_columns(at, all_columns)
    del at["AwayTeam"]

    win_loss_grid_at = dag.AgGrid(
        rowData=at.to_dict("records"),
        columnSize="responsiveSizeToFit",
        style={"height": 100, "width": "100%"},
        columnDefs=[{"field": i} for i in at.columns],
    )

    at_ag = int(at.sum(axis=1))
    at_aw = int(at["H"])
    at_awp = round(at_aw / at_ag * 100, 2)
    at_ad = int(at["D"])
    at_adp = round(at_ad / at_ag * 100, 2)
    at_al = int(at["A"])
    at_alp = round(at_al / at_ag * 100, 2)

    text_at1 = dcc.Markdown(
        f"""
    ###### **Home games : {at_ag}**
    ###### **Home wins : {at_aw} | Win rate : {at_awp}%**
    ###### **Home draws : {at_ad} | Draw rate : {at_adp}%**
    ###### **Home losses : {at_al} | Loss rate : {at_alp}%**
    """
    )

    ## team 1 away wins, losses and draws
    ht_ag = int(at.select_dtypes(include="number").sum(axis=1).sum())
    ht_aw = int(at["A"])
    ht_awp = round(ht_aw / ht_ag * 100, 2)
    ht_ad = int(at["D"])
    ht_adp = round(ht_ad / ht_ag * 100, 2)
    ht_al = int(at["H"])
    ht_alp = round(ht_al / ht_ag * 100, 2)

    text_ht2 = dcc.Markdown(
        f"""
    ###### **Away games : {ht_ag}**
    ###### **Away wins : {ht_aw} | Win rate : {ht_awp}%**
    ###### **Away draws : {ht_ad} | Draw rate : {ht_adp}%**
    ###### **Away losses : {ht_al} | Loss rate : {ht_alp}%**
    """
    )

    ## team 2 away wins, losses and draws
    at_hg = int(ht.select_dtypes(include="number").sum(axis=1).sum())
    at_hw = int(ht["A"])
    at_hwp = round(at_hw / at_hg * 100, 2)
    at_hd = int(ht["D"])
    at_hdp = round(at_hd / at_hg * 100, 2)
    at_hl = int(ht["H"])
    ht_hlp = round(at_hl / at_hg * 100, 2)

    text_at2 = dcc.Markdown(
        f"""
    ###### **Away games : {at_hg}**
    ###### **Away wins : {at_hw} | Win rate : {at_hwp}%**
    ###### **Away draws : {at_hd} | Draw rate : {at_hdp}%**
    ###### **Away losses : {at_hl} | Loss rate : {ht_hlp}%**
    """
    )

    ## team 1 total wins, losses and draws
    t1_tg = int(ht_hg + ht_ag)
    t1_tw = int(ht_hw + ht_aw)
    t1_twp = round(t1_tw / t1_tg * 100, 2)
    t1_td = int(ht_hd + ht_ad)
    t1_tdp = round(t1_td / t1_tg * 100, 2)
    t1_tl = int(ht_hl + ht_al)
    t1_tlp = round(t1_tl / t1_tg * 100, 2)

    text_t1t = dcc.Markdown(
        f"""
    ###### **Total games : {t1_tg}**
    ###### **Total wins : {t1_tw} | Win rate : {t1_twp}%**
    ###### **Total draws : {t1_td} | Draw rate : {t1_tdp}%**
    ###### **Total losses : {t1_tl} | Loss rate : {t1_tlp}%**

    """
    )

    ## team 2 total wins, losses and draws
    t2_tg = int(at_ag + at_hg)
    t2_tw = int(at_aw + at_hw)
    t2_twp = round(t2_tw / t2_tg * 100, 2)
    t2_td = int(at_ad + at_hd)
    t2_tdp = round(t2_td / t2_tg * 100, 2)
    t2_tl = int(at_al + at_hl)
    t2_tlp = round(t2_tl / t2_tg * 100, 2)

    text_t2t = dcc.Markdown(
        f"""
    ###### **Total games : {t2_tg}**
    ###### **Total wins : {t2_tw} | Win rate : {t2_twp}%**
    ###### **Total draws : {t2_td} | Draw rate : {t2_tdp}%**
    ###### **Total losses : {t2_tl} | Loss rate : {t2_tlp}%**

    """
    )

    ## team 1 goal stats
    pt1 = pd.pivot_table(
        dff,
        index=["HomeTeam"],
        values=["FullTimeHomeTeamGoals"],
        aggfunc=["sum", "count"],
    )
    pt1.columns = pt1.columns.droplevel(0)
    pt1 = pt1.reset_index()
    pt1.columns = ["Team", "HomeGoals", "Matches"]

    pt2 = pd.pivot_table(
        dff,
        index=["AwayTeam"],
        values=["FullTimeAwayTeamGoals"],
        aggfunc=["sum", "count"],
    )
    pt2.columns = pt2.columns.droplevel(0)
    pt2 = pt2.reset_index()
    pt2.columns = ["Team", "AwayGoals", "Matches"]

    dft = pd.merge(pt1, pt2, on="Team")
    dft.columns = ["Team", "HGoals", "HGames", "AGoals", "AGames"]

    dft1 = dft[dft["Team"] == t1]
    dft1["TGoals"] = dft1["HGoals"] + dft1["AGoals"]
    dft1["TGames"] = dft1["HGames"] + dft1["AGames"]
    dft1["AveGG"] = round(dft1["TGoals"] / dft1["TGames"], 2)
    dft1["AveHGG"] = round(dft1["HGoals"] / dft1["HGames"], 2)
    dft1["AveAGG"] = round(dft1["AGoals"] / dft1["AGames"], 2)
    dft1_f = dft1[dft1.columns[1:-3]]

    dft1_dag = dag.AgGrid(
        rowData=dft1.to_dict("records"),
        columnSize="responsiveSizeToFit",
        style={"height": 10, "width": "100%"},
        columnDefs=[{"field": i} for i in dft1.columns],
    )

    gg_t1 = dft1.iloc[0, 7]
    hgg_t1 = dft1.iloc[0, 8]
    agg_t1 = dft1.iloc[0, 9]

    text_t1 = dcc.Markdown(
        f"""
    ###### **Average goals per game : {gg_t1}**
    ###### **Average home goals per game : {hgg_t1}**
    ###### **Average away goals per game : {agg_t1}**
    """
    )

    dt1 = dag.AgGrid(
        rowData=dft1_f.to_dict("records"),
        columnSize="responsiveSizeToFit",
        style={"height": 100, "width": "100%"},
        columnDefs=[{"field": i} for i in dft1_f.columns],
    )

    ## team 2 goal stats
    dft2 = dft[dft["Team"] == t2]
    dft2["TGoals"] = dft2["HGoals"] + dft2["AGoals"]
    dft2["TGames"] = dft2["HGames"] + dft2["AGames"]
    dft2["AveGG"] = round(dft2["TGoals"] / dft2["TGames"], 2)
    dft2["AveHGG"] = round(dft2["HGoals"] / dft2["HGames"], 2)
    dft2["AveAGG"] = round(dft2["AGoals"] / dft2["AGames"], 2)
    dft2_f = dft2[dft2.columns[1:-3]]

    gg_t2 = dft2.iloc[0, 7]
    hgg_t2 = dft2.iloc[0, 8]
    agg_t2 = dft2.iloc[0, 9]

    text_t2 = dcc.Markdown(
        f"""
    ###### **Average goals per game : {gg_t2}**
    ###### **Average home goals per game : {hgg_t2}**
    ###### **Average away goals per game : {agg_t2}**
    """
    )

    dt2 = dag.AgGrid(
        rowData=dft2_f.to_dict("records"),
        columnSize="responsiveSizeToFit",
        style={"height": 100, "width": "100%"},
        columnDefs=[{"field": i} for i in dft2_f.columns],
    )

    ## team 1 overall stats
    bpt_df = match_data[
        [
            "Season",
            "HomeTeam",
            "AwayTeam",
            "FullTimeHomeTeamGoals",
            "FullTimeAwayTeamGoals",
        ]
    ]
    bpt_df = bpt_df[bpt_df["Season"].isin(season)]
    bpt_df["GoalDiff"] = (
        bpt_df["FullTimeHomeTeamGoals"] - bpt_df["FullTimeAwayTeamGoals"]
    )
    bpt_df

    htw_conditions = [
        (bpt_df["GoalDiff"] >= 1),
        (bpt_df["GoalDiff"] <= -1),
        (bpt_df["GoalDiff"] == 0),
    ]
    htw_values = [1, 0, 0]

    bpt_df["HomeTeamWin"] = np.select(htw_conditions, htw_values)

    atw_conditions = [
        (bpt_df["GoalDiff"] <= -1),
        (bpt_df["GoalDiff"] >= 1),
        (bpt_df["GoalDiff"] == 0),
    ]
    atw_values = [1, 0, 0]

    bpt_df["HomeTeamWin"] = np.select(htw_conditions, htw_values)
    bpt_df["AwayTeamWin"] = np.select(atw_conditions, atw_values)

    hw_df = pd.pivot_table(
        bpt_df, index=["HomeTeam"], values=["HomeTeamWin"], aggfunc=["sum", "count"]
    )
    print(hw_df)
    hw_df.columns = hw_df.columns.droplevel(0)
    print(hw_df)
    hw_df = hw_df.reset_index()
    hw_df.columns = ["Team", "HomeTeamWins", "Games"]
    hw_df["HomeWinPercentage"] = round(hw_df["HomeTeamWins"] / hw_df["Games"] * 100, 2)
    final_hw_df = hw_df.sort_values(by="HomeTeamWins", ascending=False)

    aw_df = pd.pivot_table(
        bpt_df, index=["AwayTeam"], values=["AwayTeamWin"], aggfunc=["sum", "count"]
    )
    aw_df.columns = aw_df.columns.droplevel(0)
    aw_df = aw_df.reset_index()
    aw_df.columns = ["Team", "AwayTeamWins", "Games"]
    aw_df["AwayWinPercentage"] = round(aw_df["AwayTeamWins"] / aw_df["Games"] * 100, 2)
    final_aw_df = aw_df.sort_values(by="AwayTeamWins", ascending=False)

    hd_df = bpt_df[bpt_df["GoalDiff"] == 0]
    del hd_df["AwayTeam"]
    hd_df_pivot_table = pd.pivot_table(
        hd_df, index=["HomeTeam"], values=["GoalDiff"], aggfunc=["count"]
    )

    hd_df_pivot_table.columns = hd_df_pivot_table.columns.droplevel(0)
    hd_df_pivot_table = hd_df_pivot_table.reset_index()
    hd_df_pivot_table.columns = ["Team", "HomeDraws"]

    hd_final_df = hd_df_pivot_table.sort_values(by="HomeDraws", ascending=False)

    ad_df = bpt_df[bpt_df["GoalDiff"] == 0]
    del ad_df["HomeTeam"]
    ad_df_pivot_table = pd.pivot_table(
        ad_df, index=["AwayTeam"], values=["AwayTeamWin"], aggfunc=["count"]
    )
    ad_df_pivot_table.columns = ad_df_pivot_table.columns.droplevel(0)
    ad_df_pivot_table = ad_df_pivot_table.reset_index()
    ad_df_pivot_table.columns = ["Team", "AwayDraws"]
    ad_final_df = ad_df_pivot_table.sort_values(by="AwayDraws", ascending=False)

    all_dfs = [final_hw_df, final_aw_df, hd_final_df, ad_final_df]
    comdined_df = reduce(
        lambda left, right: pd.merge(left, right, on=["Team"], how="outer"), all_dfs
    )

    comdined_df.columns = [
        "Team",
        "HTWins",
        "HGames",
        "HWinPerc",
        "ATWins",
        "AGames",
        "AWinPerc",
        "HDraws",
        "ADraws",
    ]

    final_comdined_df = comdined_df
    final_comdined_df = final_comdined_df.fillna(0)

    final_comdined_df["TGames"] = (
        final_comdined_df["HGames"] + final_comdined_df["AGames"]
    )
    final_comdined_df["HGLosses"] = (
        final_comdined_df["HGames"]
        - final_comdined_df["HTWins"]
        - final_comdined_df["HDraws"]
    )
    final_comdined_df["AGLosses"] = (
        final_comdined_df["AGames"]
        - final_comdined_df["ATWins"]
        - final_comdined_df["ADraws"]
    )
    final_comdined_df["TWins"] = (
        final_comdined_df["HTWins"] + final_comdined_df["ATWins"]
    )
    final_comdined_df["TWinPerc"] = round(
        final_comdined_df["TWins"] / final_comdined_df["TGames"] * 100, 2
    )

    final_comdined_df["Tdraws"] = (
        final_comdined_df["HDraws"] + final_comdined_df["ADraws"]
    )
    final_comdined_df["TDrawPerc"] = round(
        final_comdined_df["Tdraws"] / final_comdined_df["TGames"] * 100, 2
    )

    final_comdined_df["HDrawPerc"] = round(
        final_comdined_df["HDraws"] / final_comdined_df["HGames"] * 100, 2
    )
    final_comdined_df["ADrawPerc"] = round(
        final_comdined_df["ADraws"] / final_comdined_df["AGames"] * 100, 2
    )

    final_comdined_df["Tlosses"] = (
        final_comdined_df["HGLosses"] + final_comdined_df["AGLosses"]
    )
    final_comdined_df["TLossPerc"] = round(
        final_comdined_df["Tlosses"] / final_comdined_df["TGames"] * 100, 2
    )

    final_comdined_df["HLossPerc"] = round(
        final_comdined_df["HGLosses"] / final_comdined_df["HGames"] * 100, 2
    )
    final_comdined_df["ALossPerc"] = round(
        final_comdined_df["AGLosses"] / final_comdined_df["AGames"] * 100, 2
    )

    final_comdined_df = final_comdined_df[
        [
            "Team",
            "TGames",
            "TWins",
            "TWinPerc",
            "Tdraws",
            "TDrawPerc",
            "Tlosses",
            "TLossPerc",
            "HGames",
            "HTWins",
            "HWinPerc",
            "HDraws",
            "HDrawPerc",
            "HGLosses",
            "HLossPerc",
            "AGames",
            "ATWins",
            "AWinPerc",
            "ADraws",
            "ADrawPerc",
            "AGLosses",
            "ALossPerc",
        ]
    ]

    t1_full_stats = final_comdined_df[final_comdined_df["Team"] == t1]
    t1_total_games = t1_full_stats.iloc[0, 1]
    t1_total_wins = t1_full_stats.iloc[0, 2]
    t1_total_win_perc = t1_full_stats.iloc[0, 3]
    t1_total_draws = t1_full_stats.iloc[0, 4]
    t1_total_draw_perc = t1_full_stats.iloc[0, 5]
    t1_total_losses = t1_full_stats.iloc[0, 6]
    t1_total_loss_perc = t1_full_stats.iloc[0, 7]

    ovr_text_t1 = dcc.Markdown(
        f"""
    ###### **Total games : {t1_total_games}**
    ###### **Wins : {t1_total_wins} | Win rate : {t1_total_win_perc}%**
    ###### **Draws : {t1_total_draws} | Draw rate : {t1_total_draw_perc}%**
    ###### **Losses : {t1_total_losses} | Loss rate : {t1_total_loss_perc}%**
    """
    )

    t1_total_home_games = t1_full_stats.iloc[0, 8]
    t1_total_home_wins = t1_full_stats.iloc[0, 9]
    t1_total_home_win_perc = t1_full_stats.iloc[0, 10]
    t1_total_home_draws = t1_full_stats.iloc[0, 11]
    t1_total_home_draw_perc = t1_full_stats.iloc[0, 12]
    t1_total_home_losses = t1_full_stats.iloc[0, 13]
    t1_total_home_loss_perc = t1_full_stats.iloc[0, 14]

    ovr_hg_text_t1 = dcc.Markdown(
        f"""
    ###### **Home games : {t1_total_home_games}**
    ###### **Wins : {t1_total_home_wins} | Win rate : {t1_total_home_win_perc}%**
    ###### **Draws : {t1_total_home_draws} | Draw rate : {t1_total_home_draw_perc}%**
    ###### **Losses : {t1_total_home_losses} | Loss rate : {t1_total_home_loss_perc}%**
    """
    )

    t1_total_away_games = t1_full_stats.iloc[0, 15]
    t1_total_away_wins = t1_full_stats.iloc[0, 16]
    t1_total_away_win_perc = t1_full_stats.iloc[0, 17]
    t1_total_away_draws = t1_full_stats.iloc[0, 18]
    t1_total_away_draw_perc = t1_full_stats.iloc[0, 19]
    t1_total_away_losses = t1_full_stats.iloc[0, 20]
    t1_total_away_loss_perc = t1_full_stats.iloc[0, 21]

    ovr_ag_text_t1 = dcc.Markdown(
        f"""
    ###### **Away Games : {t1_total_away_games}**
    ###### **Wins : {t1_total_away_wins} | Win rate : {t1_total_away_win_perc}%**
    ###### **Draws : {t1_total_away_draws} | Draw rate : {t1_total_away_draw_perc}%**
    ###### **Losses : {t1_total_away_losses} | Loss rate : {t1_total_away_loss_perc}%**
    """
    )

    ## team 2 overall stats
    t2_full_stats = final_comdined_df[final_comdined_df["Team"] == t2]
    t2_total_games = t2_full_stats.iloc[0, 1]
    t2_total_wins = t2_full_stats.iloc[0, 2]
    t2_total_win_perc = t2_full_stats.iloc[0, 3]
    t2_total_draws = t2_full_stats.iloc[0, 4]
    t2_total_draw_perc = t2_full_stats.iloc[0, 5]
    t2_total_losses = t2_full_stats.iloc[0, 6]
    t2_total_loss_perc = t2_full_stats.iloc[0, 7]

    ovr_text_t2 = dcc.Markdown(
        f"""
    ###### **Total games : {t2_total_games}**
    ###### **Wins : {t2_total_wins} | Win rate : {t2_total_win_perc}%**
    ###### **Draws : {t2_total_draws} | Draw rate : {t2_total_draw_perc}%**
    ###### **Losses : {t2_total_losses} | Loss rate : {t2_total_loss_perc}%**
    """
    )

    t2_total_home_games = t2_full_stats.iloc[0, 8]
    t2_total_home_wins = t2_full_stats.iloc[0, 9]
    t2_total_home_win_perc = t2_full_stats.iloc[0, 10]
    t2_total_home_draws = t2_full_stats.iloc[0, 11]
    t2_total_home_draw_perc = t2_full_stats.iloc[0, 12]
    t2_total_home_losses = t2_full_stats.iloc[0, 13]
    t2_total_home_loss_perc = t2_full_stats.iloc[0, 14]

    ovr_hg_text_t2 = dcc.Markdown(
        f"""
    ###### **Home games : {t2_total_home_games}**
    ###### **Wins : {t2_total_home_wins} | Win rate : {t2_total_home_win_perc}%**
    ###### **Draws : {t2_total_home_draws} | Draw rate : {t2_total_home_draw_perc}%**
    ###### **Losses : {t2_total_home_losses} | Loss rate : {t2_total_home_loss_perc}%**
    """
    )

    t2_total_away_games = t2_full_stats.iloc[0, 15]
    t2_total_away_wins = t2_full_stats.iloc[0, 16]
    t2_total_away_win_perc = t2_full_stats.iloc[0, 17]
    t2_total_away_draws = t2_full_stats.iloc[0, 18]
    t2_total_away_draw_perc = t2_full_stats.iloc[0, 19]
    t2_total_away_losses = t2_full_stats.iloc[0, 20]
    t2_total_away_loss_perc = t2_full_stats.iloc[0, 21]

    ovr_ag_text_t2 = dcc.Markdown(
        f"""
    ###### **Away Games : {t2_total_away_games}**
    ###### **Wins : {t2_total_away_wins} | Win rate : {t2_total_away_win_perc}%**
    ###### **Draws : {t2_total_away_draws} | Draw rate : {t2_total_away_draw_perc}%**
    ###### **Losses : {t2_total_away_losses} | Loss rate : {t2_total_away_loss_perc}%**
    """
    )

    ## team 1 ovr goal stats
    dff_ovr = match_data.copy()
    dff_ovr = dff_ovr[dff_ovr["Season"].isin(season)]
    t1_ovr_goals_df = pd.pivot_table(
        dff_ovr,
        index=["HomeTeam"],
        values=["FullTimeHomeTeamGoals"],
        aggfunc=["sum", "count"],
    )
    t1_ovr_goals_df.columns = t1_ovr_goals_df.columns.droplevel(0)
    t1_ovr_goals_df = t1_ovr_goals_df.reset_index()
    t1_ovr_goals_df.columns = ["Team", "HomeGoals", "Matches"]

    t2_ovr_goals_df = pd.pivot_table(
        dff_ovr,
        index=["AwayTeam"],
        values=["FullTimeAwayTeamGoals"],
        aggfunc=["sum", "count"],
    )
    t2_ovr_goals_df.columns = t2_ovr_goals_df.columns.droplevel(0)
    t2_ovr_goals_df = t2_ovr_goals_df.reset_index()
    t2_ovr_goals_df.columns = ["Team", "AwayGoals", "Matches"]

    dftovr = pd.merge(t1_ovr_goals_df, t2_ovr_goals_df, on="Team")
    dftovr.columns = ["Team", "HGoals", "HGames", "AGoals", "AGames"]

    dft1over = dftovr[dftovr["Team"] == t1]
    dft1over["TGoals"] = dft1over["HGoals"] + dft1over["AGoals"]
    dft1over["TGames"] = dft1over["HGames"] + dft1over["AGames"]
    dft1over["AveGG"] = round(dft1over["TGoals"] / dft1over["TGames"], 2)
    dft1over["AveHGG"] = round(dft1over["HGoals"] / dft1over["HGames"], 2)
    dft1over["AveAGG"] = round(dft1over["AGoals"] / dft1over["AGames"], 2)

    ovr_gg_t1 = dft1over.iloc[0, 7]
    ovr_hgg_t1 = dft1over.iloc[0, 8]
    ovr_agg_t1 = dft1over.iloc[0, 9]

    ovr_t1_goals = dcc.Markdown(
        f"""
    ###### **Average goals per game : {ovr_gg_t1}**
    ###### **Average home goals per game : {ovr_hgg_t1}**
    ###### **Average away goals per game : {ovr_agg_t1}**
    """
    )

    ## team 2 ovr goal stats
    t2_ovr_goals_df = dftovr[dftovr["Team"] == t2]

    dft2over = dftovr[dftovr["Team"] == t2]
    dft2over["TGoals"] = dft2over["HGoals"] + dft2over["AGoals"]
    dft2over["TGames"] = dft2over["HGames"] + dft2over["AGames"]
    dft2over["AveGG"] = round(dft2over["TGoals"] / dft2over["TGames"], 2)
    dft2over["AveHGG"] = round(dft2over["HGoals"] / dft2over["HGames"], 2)
    dft2over["AveAGG"] = round(dft2over["AGoals"] / dft2over["AGames"], 2)

    ovr_gg_t2 = dft2over.iloc[0, 7]
    ovr_hgg_t2 = dft2over.iloc[0, 8]
    ovr_agg_t2 = dft2over.iloc[0, 9]

    ovr_t2_goals = dcc.Markdown(
        f"""
    ###### **Average goals per game : {ovr_gg_t2}**
    ###### **Average home goals per game : {ovr_hgg_t2}**
    ###### **Average away goals per game : {ovr_agg_t2}**
    """
    )

    return (
        h2h_grid,
        text_t1,
        text_ht1,
        text_ht2,
        text_t1t,
        ovr_hg_text_t1,
        ovr_ag_text_t1,
        ovr_text_t1,
        ovr_t1_goals,
        ovr_hg_text_t2,
        ovr_ag_text_t2,
        ovr_text_t2,
        text_at1,
        text_at2,
        text_t2t,
        text_t2,
        ovr_t2_goals,
    )


if __name__ == "__main__":
    app.run(debug=True)
