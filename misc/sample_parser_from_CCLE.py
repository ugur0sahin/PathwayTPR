import json
import pandas as pd

CCLE_mutation_dbs = pd.read_csv("../dbs/CCLE_mutations.csv")
DepMap_Cellines_dbs  = pd.read_csv("../dbs/DepMap-2018q3-celllines.csv")

def render_case_row(Profile_dbs, Tumor_Barcode, MutationProfile):

    case_ID_dbs = Profile_dbs[Profile_dbs["Broad_ID"] == Tumor_Barcode].to_dict()

    values = list(case_ID_dbs.values())
    rendered_case_dict =  {k: next(iter(v.values())) for k, v in case_ID_dbs.items()}

    del rendered_case_dict["Broad_ID"]
    rendered_case_dict["MutationProfile"]=MutationProfile

    return rendered_case_dict

if __name__ == '__main__':
    for chosen_sample_ID in list(set(CCLE_mutation_dbs["DepMap_ID"].to_list())):
        try:
            mutation_profile_of_Sample= dict()
            for index, row in CCLE_mutation_dbs[CCLE_mutation_dbs["DepMap_ID"] == chosen_sample_ID].iterrows():
                try:
                    mutation_profile_of_Sample[row["Hugo_Symbol"]+"_"+str(row["Protein_Change"]).split(".")[1]] = {"isDeleterious":row["isDeleterious"],
                                                                                    "isTCGAhotspot":row["isTCGAhotspot"],
                                                                                    "TCGAhsCnt":row["TCGAhsCnt"],
                                                                                    "isCOSMIChotspot":row["isCOSMIChotspot"],
                                                                                    "COSMIChsCnt":row["COSMIChsCnt"]}

                except:
                    mutation_profile_of_Sample[str(row["Hugo_Symbol"])+"_InDel"] = {"isDeleterious":row["isDeleterious"],
                                                                                    "isTCGAhotspot":row["isTCGAhotspot"],
                                                                                    "TCGAhsCnt":row["TCGAhsCnt"],
                                                                                    "isCOSMIChotspot":row["isCOSMIChotspot"],
                                                                                    "COSMIChsCnt":row["COSMIChsCnt"]}

            sample_file = open("../example/"+str(chosen_sample_ID)+".json","w")
            json.dump(render_case_row(DepMap_Cellines_dbs,chosen_sample_ID,mutation_profile_of_Sample),sample_file)
            sample_file.close()
        except:
            print(chosen_sample_ID)



