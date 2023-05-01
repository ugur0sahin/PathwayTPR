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

def affectedPathway_finder(MutationSet,justOncoPaths=False,ScoreBaseModel=1, justReturnPathways=False ,differentiationIndicatorModel="Default"):
    """Section1 -> Iterate All Mutations passed from the filtering to find it's related Pathways, add all Pathways into the
    affectedPathways_forSample and do Pathways as keys to define related Mutations as list in the definedMutations_toPathways dictionary."""
    from src.scoreNetwork import main_disruption_rateGraph

    with open("dbs/indexed_ConsensusPDB.json", "r") as indexed_ConsensusPBD_fl:
        indexed_ConsensusPBD = json.load(indexed_ConsensusPBD_fl)

    affectedPathways_forSample, definedMutations_toPathways = list(), dict()
    for Mutation in MutationSet:
        Gene = Mutation.split("_")[0]
        try:
            objectedMutation_ofPatway = indexed_ConsensusPBD[Gene]
            affectedPathways_forSample.extend(objectedMutation_ofPatway)
        except:
            continue

        from src.scoreNetwork import main_disruption_rateGraph

        for Pathway in objectedMutation_ofPatway:
            if Pathway not in definedMutations_toPathways:
                definedMutations_toPathways[Pathway] = []
            definedMutations_toPathways[Pathway].append(Gene)

    """ This is the block which is the pathway provider into another pipelines, if program go in alternative way rather than scoring"""
    if justReturnPathways:
        return definedMutations_toPathways

    """Section2 -> The affectedPathways_forSample and definedMutations_toPathways are constructed 
    ScoreBaseModel 1 uses affectedPathways_forSample just counting the collected Mutations which impacts Pathways, 
    the counts of Mutations belonged Pathways indicates degree of disruption. 
     
     ScoreBaseModel 2 and 3 are us Thoing definedMutations_toPathways dictionary contain <Pathway> : list(Mutations),
     ScoreBaseModel 2 send dictionary into main_disruption_rateGraph(), ScoreBaseModel 3 send calculate_affectMutation_intoPathway()
     
     The expected output ->  Both of Them are <Pathway> : <DisruptionScore>
     
     """
    if ScoreBaseModel == 1: #WORKING!
        return Counter(affectedPathways_forSample)

    elif ScoreBaseModel == 2: #NOT WORKING YET!
        return main_disruption_rateGraph(definedMutations_toPathways, differentiationIndicatorModel=differentiationIndicatorModel)

    elif ScoreBaseModel == 3: #WORKING!
        from src.scoreNetwork import calculate_affectMutation_intoPathway
        Counter_Analog_dict = dict()
        for Pathway, MutationSet in definedMutations_toPathways.items():
            try:
                returnedScore = calculate_affectMutation_intoPathway(Pathway, MutationSet)
                print(Pathway,len(MutationSet) , returnedScore)
                Counter_Analog_dict[Pathway] = returnedScore
            except:
                pass
        return Counter_Analog_dict
"""
if __name__ != '__main__':
    with open("../example/ACH-000219.json","r") as case_json_fl:
        case_features = json.load(case_json_fl)

    filteredMutations = filterMutations(case_features["MutationProfile"],True,False,False)
    print(affectedPathway_finder(filteredMutations))
"""