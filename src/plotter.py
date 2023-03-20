import warnings
import pandas as pd
import plotly.io as pio
import plotly.express as px

warnings.filterwarnings('ignore')

def plot_as_radar_chart(disrupted_pathway_result_as_dct, maxNumberFeature, rawCountTreshold, showCaseInfo, Output):
    pathways, counts = list(),list()
    for k,v in disrupted_pathway_result_as_dct.items():
        if rawCountTreshold is None or v >= rawCountTreshold:
                pathways.append(k)
                counts.append(v)

    be_plot_dbs = pd.DataFrame({"Pathways":pathways,"Disrupted_Score":counts})
    #print(be_plot_dbs)

    try:
        fig = px.line_polar(be_plot_dbs, r = "Disrupted_Score", theta="Pathways",line_close=True )
        fig.update_traces(fill='toself')
        fig.show()

        if Output is not None:
            try:
                pio.write_image(fig,Output)
            except:
                print("It doesn't write !")
    except:
        print("No Data remained to show ! because of harsh filter.")
        pass

    """
    note = "Trial: Trial<br>Trial: Trial"
    fig.update_layout(
        annotations=[
            dict(
                x=1,
                y=1,
                xref="paper",
                yref="paper",
                text=note,
                showarrow=False,
                font=dict(size=14)
            )
        ]
    )
    """
