import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def checkPathway_Interactome_compatability(objectedPathway,G):
    subGraph = nx.subgraph(G, objectedPathway) # Replace 'your_main_graph' with the main graph object
    isolatedNodes = list(nx.isolates(subGraph))

    # Create a new graph object with the same nodes and edges as the subGraph
    modifiable_subGraph = nx.Graph(subGraph)

    # Remove isolated nodes from the new modifiable_subGraph
    modifiable_subGraph.remove_nodes_from(isolatedNodes)

    return modifiable_subGraph, isolatedNodes


def grandGraph_build(InteractionDBS, filterConfidence=None):
    if filterConfidence is not None:
        InteractionDBS = InteractionDBS[InteractionDBS["Confidence Value"] >= filterConfidence]

    InteractionList = [(row["Gene Name Interactor A"], row["Gene Name Interactor B"]) for index, row in
                   InteractionDBS.iterrows()]

    return nx.Graph(InteractionList)

def Pathway_parser(PathwayDBS):
    Pathway_defined_HugoSymb = dict()
    for index,row in PathwayDBS.iterrows():
        Pathway_defined_HugoSymb[row["pathway"]] = row["hgnc_symbol_ids"].split(",")
    return Pathway_defined_HugoSymb

HIPPIE_dbs, CPDB_dbs = pd.read_csv("../dbs/HIPPIE-2.2.mitab.txt",delimiter="\t"), \
                       pd.read_csv("../dbs/CPDB_pathways_genes.tsv",delimiter="\t")

def drawGraph(G,Pathway_name,saveGraphML=None, plot_show = None):
    if saveGraphML is not None:
        try:
            nx.write_gexf(G,"../dbs/Graphlets/"+Pathway_name.replace(" ","_")+".gexf")
        except:
            print(Pathway_name)

    if plot_show is not None:
        nx.draw(G, with_labels=True, node_color='lightblue', font_weight='bold', node_size=500)
        plt.show()

if __name__ == '__main__':

    grandGraph_overallPathways = grandGraph_build(HIPPIE_dbs)

    for Pathway_name,HugoSymb_ls in Pathway_parser(CPDB_dbs).items():
        objectedPathway = grandGraph_overallPathways.subgraph(HugoSymb_ls)
        isolatedGraph_of_Pathway, isolatedNodes = checkPathway_Interactome_compatability(objectedPathway,grandGraph_overallPathways)
        drawGraph(isolatedGraph_of_Pathway, Pathway_name,saveGraphML=True)
