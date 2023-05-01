import json
import networkx as nx
import numpy as np
import pandas as pd

def implementMutations_vectorizedPathway(Pathway_definedMutations, Pathway):
    # Create filename for the graphlet
    Graphlet_flname = str(Pathway).replace(" ", "_") + ".gexf"
    # Try to read the graphlet from the file
    try:
        Gpre = nx.read_gexf("dbs/Graphlets/" + Graphlet_flname)
    except:
        # If the graphlet file doesn't exist, return an empty list of alterations and genes
        return nullVector_todefineAlterations.tolist()[0], indexedGenes

    # Create a sorted list of the gene names in the graphlet
    indexedGenes = sorted(list(Gpre.nodes))
    # Create a null vector of alterations for the genes in the graphlet
    nullVector_todefineAlterations = np.zeros((1, len(Gpre.nodes)), dtype=int)

    # Set the value to 1 in the null vector for each mutation in the pathway that is also in the graphlet
    try:
        for Mutation in Pathway_definedMutations:
            try:
                nullVector_todefineAlterations[0, indexedGenes.index(Mutation)] = 1
            except:
                pass
    except:
        pass

    # Return the null vector and indexed genes as lists
    return nullVector_todefineAlterations.tolist()[0], indexedGenes


def buildMatrix_alteredPathway_ofSamples(Mutations_ofSamples_definedPathways):
    Pathways, Samples = list(set([Pathway for Sample in Mutations_ofSamples_definedPathways.keys() for Pathway in Mutations_ofSamples_definedPathways[Sample].keys()])), Mutations_ofSamples_definedPathways.keys()

    AllPathways_alterationMatrices = dict()
    for Pathway in Pathways:
        alterationMatrix_ofPathway = dict()
        for Sample in Samples:
            try:
                Pathway_definedMutations =Mutations_ofSamples_definedPathways[Sample][Pathway]
                alterationMatrix_ofPathway[Sample], indexedGenes = implementMutations_vectorizedPathway(Pathway_definedMutations,
                                                                                          Pathway)

            except:
                alterationMatrix_ofPathway[Sample], indexedGenes = implementMutations_vectorizedPathway(None,
                                                                                          Pathway)

        AllPathways_alterationMatrices[Pathway] = pd.DataFrame(alterationMatrix_ofPathway, index=indexedGenes)

    #json.dump(AllPathways_alterationMatrices, open("trial.json","w"))
    return AllPathways_alterationMatrices


if __name__ == '__main__':
    pass