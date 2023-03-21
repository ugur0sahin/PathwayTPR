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

def affected_pathway_finder(MutationSet,justOncoPaths=False,rawCountScore=False):
    with open("/Users/ugur.sahin/PycharmProjects/TumorCharacterization/dbs/indexed_ConsensusPDB.json", "r") as indexed_ConsensusPBD_fl:
        indexed_ConsensusPBD = json.load(indexed_ConsensusPBD_fl)

    affected_pathways_for_sample, definedMutations_toPathways = list(), dict()
    for Mutation in MutationSet:
        Gene = Mutation.split("_")[0]
        try:
            objectedMutation_ofPatway = indexed_ConsensusPBD[Gene]
            affected_pathways_for_sample.extend(objectedMutation_ofPatway)
        except:
            continue

        if not rawCountScore:
            from src.Network_score import main_disruption_rateGraph

            for Pathway in objectedMutation_ofPatway:
                if Pathway not in definedMutations_toPathways:
                    definedMutations_toPathways[Pathway] = []
                definedMutations_toPathways[Pathway].append(Gene)

    if rawCountScore:
        return Counter(affected_pathways_for_sample)
    else:
        return main_disruption_rateGraph(definedMutations_toPathways)


if __name__ != '__main__':
    with open("../example/ACH-000219.json","r") as case_json_fl:
        case_features = json.load(case_json_fl)

    filteredMutations = filterMutations(case_features["MutationProfile"],True,False,False)
    print(affected_pathway_finder(filteredMutations))
