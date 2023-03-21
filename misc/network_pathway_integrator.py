import pandas as pd
import networkx as nx

HIPPIE_dbs = pd.read_csv("../dbs/HIPPIE-2.2.mitab.txt",delimiter="\t")
HIPPIE_dbs_filtered = HIPPIE_dbs[HIPPIE_dbs["Confidence Value"] >= 0.5]

InteractionList = [(row["Gene Name Interactor A"], row["Gene Name Interactor B"]) for index,row in HIPPIE_dbs_filtered.iterrows()]
G = nx.Graph(InteractionList)
print(G)

def grandGraphBuild(InteractionDBS):
    pass