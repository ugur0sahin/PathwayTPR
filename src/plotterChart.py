import random
import warnings
import pandas as pd
import plotly.io as pio
import plotly.express as px
import plotly.colors as plc
from sklearn.decomposition import PCA
import numpy as np



warnings.filterwarnings('ignore')

def plotRadarChart(disrupted_pathway_result_as_dct, maxNumberFeature, rawCountTreshold, Output):
    pathways, counts = list(),list()
    for k,v in disrupted_pathway_result_as_dct.items():
        if rawCountTreshold is None or int(float(v)) >= rawCountTreshold:
                try:
                    pathways.append(k.replace(" - Homo sapiens (human)","- KEGG"))
                except:
                    pathways.append(k)

                counts.append(float(v))


    be_plot_dbs = pd.DataFrame({"Pathways":pathways,"Disrupted_Score":counts})

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

def plotRadarChart_multipleSamples(disruptedPathways, maxNumberFeature, rawCountTreshold, Output):
    # Merge all pathways
    mergedPathways = list()
    for Sample in disruptedPathways.keys():
        mergedPathways.extend(list(disruptedPathways[Sample].keys()))
    mergedPathways = list(set(mergedPathways))

    # Create a list of lists to store pathway values for each sample
    collectionSample = list()
    for Pathway in mergedPathways:
        collectionPathway_tmp = list()
        for Sample in disruptedPathways.keys():
            try:
                collectionPathway_tmp.append(np.log10(float(disruptedPathways[Sample][Pathway]) + 1))
            except:
                collectionPathway_tmp.append(0)

        if any(disruptionScore > rawCountTreshold for disruptionScore in collectionPathway_tmp):
            collectionPathway_tmp.insert(0, Pathway.replace(" - Homo sapiens (human)", "- KEGG"))
            collectionSample.append(collectionPathway_tmp)
        else:
            pass

    # Create dataframe with pathway values for each sample
    columns = [Sample.replace(".json", "") for Sample in disruptedPathways.keys()]
    columns.insert(0, "Pathway")
    df = pd.DataFrame(collectionSample, columns=columns)

    # Melt dataframe and apply log10 transformation to pathway values
    df_melted = df.melt(id_vars=['Pathway'], var_name='Sample', value_name='Value')

    # Plot pathway disruption comparison using polar line chart
    fig = px.line_polar(df_melted, r='Value', theta='Pathway', color='Sample', line_close=True)
    fig.update_traces(fill='toself')

    fig.update_layout(
        title='Pathway Disruption Comparison - Log10Scale',
        polar=dict(
            radialaxis=dict(visible=True, type='log', tickformat='.2f'),
            angularaxis=dict(visible=True)
        ),
        showlegend=True
    )

    fig.show()


def showPCA_basedPathway(Matrix_of_alteredPathway,plotName="defaultPCA.html"):
    sample_names = list(Matrix_of_alteredPathway.columns)
    groups = [sample.split("|")[0].split("_")[1] for sample in sample_names]

    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(Matrix_of_alteredPathway.T)
    principal_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=sample_names)

    # Sample isimlerini ve grupları DataFrame'e ekleyin
    principal_df['Sample'] = sample_names
    principal_df['Group'] = groups

    # Renkleri rastgele ata
    unique_groups = set(groups)
    colors = {group: "#" + "".join([random.choice("0123456789ABCDEF") for j in range(6)]) for group in unique_groups}

    # Scatter plot oluşturun ve renkleri özelleştirin
    fig = px.scatter(principal_df, x='PC1', y='PC2', hover_name='Sample',
                     color='Group', color_discrete_map=colors)

    # Plot'u gösterin
    fig.update_layout(title=plotName.split(".html")[0],width=1000, height=1000)
    fig.write_html(plotName)
    #fig.show()

def get_value_by_substring(dictionary, substring):
    for key, value in dictionary.items():
        if substring in key:
            return value
    return None