import json
import argparse
from src.identifierSample import filterMutations, affected_pathway_finder
from src.plotterChart import plot_as_radar_chart


parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--Sample', default=None, help='...')
parser.add_argument('--Output', default=None, help='...')
parser.add_argument('--maxNumberFeature', default=None, type=int, help='...')
parser.add_argument('--rawCountTreshold', default=None, type=int, help='...')
parser.add_argument('--showCaseInfo', default=False, action='store_true', help='...')
parser.add_argument('--isDeleterious', default=False, action='store_true', help='...')
parser.add_argument('--isCOSMIChotspot', default=False, action='store_true', help='...')
parser.add_argument('--isTCGAhotspot', default=False, action='store_true', help='...')
parser.add_argument('--justOncoPaths', default=False, action='store_false', help='...')
#parser.add_argument('--rawCountScore', default=False, action='store_false', help='...')
parser.add_argument('--differentiationIndicatorModel', default="Default", type=str, help='...')
parser.add_argument('--scoreBaseModel', default=1, type=int, help='...')


args = parser.parse_args()


def main():

    with open(args.Sample,"r") as case_json_fl:
        case_features = json.load(case_json_fl)

    filteredMutations = filterMutations(case_features["MutationProfile"],args.isDeleterious, args.isCOSMIChotspot, args.isTCGAhotspot)

    counted_affected_path_dict = affected_pathway_finder(filteredMutations,
                                                         ScoreBaseModel=args.scoreBaseModel,
                                                         differentiationIndicatorModel=args.differentiationIndicatorModel)

    plot_as_radar_chart(counted_affected_path_dict,args.maxNumberFeature, args.rawCountTreshold, args.showCaseInfo, args.Output)



if __name__ == '__main__':
    main()

