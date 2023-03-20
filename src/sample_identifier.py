import json
from collections import Counter

def filterMutations(MutationProfile,isDeleterious,isCOSMIChotspot,isTCGAhotspot):

    nested_Mutations = list(MutationProfile.keys())

    if isDeleterious:
        isDeleterious_genes = list()
        for gene, mutations in MutationProfile.items():
            if mutations["isDeleterious"]:
                isDeleterious_genes.append(gene)
        nested_Mutations = set(isDeleterious_genes).intersection(set(nested_Mutations))


    if isCOSMIChotspot:
        isCOSMIChotspot_genes = list()
        for gene, mutations in MutationProfile.items():
            if mutations["isCOSMIChotspot"]:
                isCOSMIChotspot_genes.append(gene)

        nested_Mutations = set(isCOSMIChotspot_genes).intersection(set(nested_Mutations))


    if isTCGAhotspot:
        isTCGAhotspot_genes = list()
        for gene, mutations in MutationProfile.items():
            if mutations["isTCGAhotspot"]:
                isTCGAhotspot_genes.append(gene)

        nested_Mutations = set(isTCGAhotspot_genes).intersection(set(nested_Mutations))


    return nested_Mutations

def affected_pathway_finder(MutationSet,justOncoPaths=False):
    with open("/Users/ugur.sahin/PycharmProjects/TumorCharacterization/dbs/indexed_ConsensusPDB.json", "r") as indexed_ConsensusPBD_fl:
        indexed_ConsensusPBD = json.load(indexed_ConsensusPBD_fl)

    affected_pathways_for_sample = list()
    for mutation in MutationSet:
        gene = mutation.split("_")[0]
        try:
            affected_pathways_for_sample.extend(indexed_ConsensusPBD[gene])
        except:
            #print("There is no Any defined path to this Gene.")
            pass
    return Counter(affected_pathways_for_sample)

"""
if __name__ != '__main__':
    with open("../misc/ACH-000219.json","r") as case_json_fl:
        case_features = json.load(case_json_fl)

    filteredMutations = filterMutations(case_features["MutationProfile"])
    affected_pathway_finder(filteredMutations)
"""