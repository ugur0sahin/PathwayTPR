import json
import random
import argparse
import pandas as pd
from src.identifierSample import filterMutations, affectedPathway_finder
from src.plotterChart import plotRadarChart, plotRadarChart_multipleSamples, showPCA_basedPathway, get_value_by_substring


parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--Sample',  nargs='+', default=None, help='...')
parser.add_argument('--Output', default=None, help='...')
parser.add_argument('--maxNumberFeature', default=None, type=int, help='...')
parser.add_argument('--rawCountTreshold', default=0, type=float, help='...')
parser.add_argument('--showCaseInfo', default=False, action='store_true', help='...')
parser.add_argument('--isDeleterious', default=False, action='store_true', help='...')
parser.add_argument('--isCOSMIChotspot', default=False, action='store_true', help='...')
parser.add_argument('--isTCGAhotspot', default=False, action='store_true', help='...')
parser.add_argument('--justOncoPaths', default=False, action='store_false', help='...')
parser.add_argument('--differentiationIndicatorModel', default="Default", type=str, help='...')
parser.add_argument('--scoreBaseModel', default=1, type=int, help='...')
parser.add_argument('--doSimilarity', default=False, action='store_false', help='...')
parser.add_argument('--showReducedDimension', nargs='+',default="Default", type=str, help='...')


args = parser.parse_args()

def main():

    collectionSamples, indicatedMutations_ofSamples_of_definedPathways = dict(), dict()
    for SamplePath in args.Sample:
        try:
            with open(SamplePath, "r") as case_json_fl:
                case_features = json.load(case_json_fl)
        except:
            continue

        filteredMutations = filterMutations(case_features["MutationProfile"], args.isDeleterious, args.isCOSMIChotspot,
                                            args.isTCGAhotspot)

        counted_affected_path_dict, definedMutations_toPathways = affectedPathway_finder(filteredMutations,
                                                            ScoreBaseModel=args.scoreBaseModel,
                                                            differentiationIndicatorModel=args.differentiationIndicatorModel), \
                                                                  affectedPathway_finder(filteredMutations, justReturnPathways=True)

        collectionSamples[str(case_features["CCLE_Name"])+"|"+str(case_features["Subtype_Disease"])], \
            indicatedMutations_ofSamples_of_definedPathways[str(case_features["CCLE_Name"])+"|"+str(case_features["Subtype_Disease"])] =\
            counted_affected_path_dict, definedMutations_toPathways

        case_json_fl.close()

    if len(collectionSamples.keys()) == 1:
        plotRadarChart(collectionSamples[list(collectionSamples.keys())[0]], args.maxNumberFeature, args.rawCountTreshold, args.Output)

    else:
        plotRadarChart_multipleSamples(collectionSamples, args.maxNumberFeature, args.rawCountTreshold, args.Output)


    if args.showReducedDimension != "Default":
        from src.architectPathway import buildMatrix_alteredPathway_ofSamples
        alterationMatrices_ofAllPathways_definedSamples = buildMatrix_alteredPathway_ofSamples(
            indicatedMutations_ofSamples_of_definedPathways)

        for Pathway in args.showReducedDimension:
            print(Pathway)
            #Matrix_of_alteredPathway = alterationMatrices_ofAllPathways_definedSamples[Pathway.replace("_"," ")]
            Matrix_of_alteredPathway=get_value_by_substring(alterationMatrices_ofAllPathways_definedSamples, Pathway.replace("_"," "))
            showPCA_basedPathway(Matrix_of_alteredPathway,plotName=Pathway+".PCA.html")


if __name__ == '__main__':
    main()

