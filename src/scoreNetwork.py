import json
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

"""
==========================================================================================
Section1A : Defining Score Models to Measure Affect of Mutations into the Mutated Pathways ![Network Based Scoring]!
==========================================================================================

All Mutations introduced into the pathways cause disruptions of the Hubs, decreases Flow and Reachability,
Mutations takes places varied Genes(Nodes) have varied degree impacts to Pathways.

#Model 1 -> Betweeness Centrality for each node in both WTGraph and MTGraph, you can compare their values
to understand the changes in the pathways. A significant change in the betweenness centrality of a gene (node) may indicate
that the mutation has altered the gene's role within the pathway, affecting the overall information flow and functionality of the system.

#Model 2 -> Flow change between WTGraph and MTGraph, by measuring the difference in average shortest path
lengths between nodes in the two graphs.The resulting value, called "flow change," indicates the change in 
how easily nodes can communicate with each other in the mutated graph compared to the healthy graph. A positive 
change suggests that the flow between nodes is slower in the mutated graph, mwhich may indicate a significant impact
of the mutation on the pathway.

#Model 3 -> The kernel_similarity function calculates the similarity between two pathway graphs using graph kernel methods
(Euclidean, Manhattan, RBF). It compares feature vectors based on topological properties and returns similarity scores for
each kernel to assess structural differences.

"""

def Model1_betweenness_centrality(Gpre, Gpost):
     Gpre_betweenness, Gpost_betweenness, change_ratios = nx.betweenness_centrality(Gpre), nx.betweenness_centrality(Gpost), {}
     for Gene in Gpre_betweenness:
         if Gene in Gpost_betweenness:
             WTvalue = Gpre_betweenness[Gene]
             MTvalue = Gpost_betweenness[Gene]
             try:
                 dChange = (MTvalue - WTvalue) / WTvalue
             except:
                 dChange = 0
             change_ratios[Gene] = dChange

     return sum(list(change_ratios.values()))

def Model2_flow_change_w_average_shortest_path(Gpre, Gpost):

    def ASPlen(Graph):
        if nx.is_connected(Graph):
           return nx.average_shortest_path_length((Graph))
        else:
           components, lengths = nx.connected_components(Graph), list()
           for component in components:
               subGraph = Graph.subgraph(component)
               lengths.append(nx.average_shortest_path_length(subGraph))
           return sum(lengths) / len(lengths)

    def flow_change(Gpre, Gpost):
        Gpre_AVGShortestPath, Gpost_AVGShortestPath = ASPlen(Gpre), ASPlen(Gpost)

        return Gpost_AVGShortestPath - Gpre_AVGShortestPath

    return flow_change(Gpre, Gpost)

def Model3_kernel_similarity(Gpre, Gpost, distanceMetric=None):
    def euclidean_distance(vec1, vec2):
        return np.sqrt(np.sum((vec1 - vec2) ** 2))

    def manhattan_distance(vec1, vec2):
        return np.sum(np.abs(vec1 - vec2))

    unique_nodes = sorted(set(Gpre.nodes()).union(set(Gpost.nodes())))

    healthy_vector = np.array([1 if node in Gpre.nodes() else 0 for node in unique_nodes])
    mutated_vector = np.array([1 if node in Gpost.nodes() else 0 for node in unique_nodes])

    if distanceMetric == 'kernelSimilarityEuclidian':
        return euclidean_distance(healthy_vector, mutated_vector)
    elif distanceMetric == 'kernelSimilarityManhattan':
        return manhattan_distance(healthy_vector, mutated_vector)
    else:
        raise ValueError("Invalid distanceMetric. Supported values are 'euclidean' and 'manhattan'.")


"""
==========================================================================================
Section1B : Prepared Mutated Pathways through Knocking-out Mutated Nodes from the Pathway(s) ![Network Based Scoring]! .
==========================================================================================

# main_disruption_rateGraph() is a roof function contains subfunctions denoted below, 
|                                              It is directly called in the Other modules fully spans this script.
|
|----> # callPathwayGraphlet_toremoveNodes() is responsible to finding Pathways of interested, after it removes Mutated Nodes.
|----> # disruption_score_indicator() is responsible to scoring Graphs of Pathways based on chosen measuring Model.

"""

def callPathwayGraphlet_toremoveNodes(Pathway, tobe_removedNode_ls):
    #It will give outputs as preGraph and removed(processed)Graph
    Graphlet_flname = str(Pathway).replace(" ", "_") + ".gexf"

    try:
        Gpre = nx.read_gexf("dbs/Graphlets/"+Graphlet_flname) #call Pathway and define into NetworkX Graph
        Gpost = Gpre.copy()

        # Remove nodes from Graphlet
        Gpost.remove_nodes_from(tobe_removedNode_ls)
    except:
        return None
        pass

    return Gpre,Gpost

def disruption_score_indicator(Gpre,Gpost,differentiationIndicatorModel="Default"):
    if differentiationIndicatorModel == "betweennessCentrality":
        return Model1_betweenness_centrality(Gpre, Gpost)

    elif differentiationIndicatorModel == "flowChangeASP":
        return Model2_flow_change_w_average_shortest_path(Gpre, Gpost)

    elif "kernelSimilarity" in differentiationIndicatorModel:

        if differentiationIndicatorModel == "kernelSimilarityEuclidian":
           return Model3_kernel_similarity(Gpre, Gpost, distanceMetric="kernelSimilarityEuclidian")
        else:
           return Model3_kernel_similarity(Gpre, Gpost, distanceMetric="kernelSimilarityManhattan")


def main_disruption_rateGraph(tobe_removedNodes_definedPathways, differentiationIndicatorModel="kernelSimilarityEuclidian"):
    # Not actually filtered!
    definedMutations_toPathways_filtered = tobe_removedNodes_definedPathways

    Pathway_disruption_rates = dict()
    for Pathway, tobe_removedNode_ls in definedMutations_toPathways_filtered.items():
     try:
            Gpre, Gpost = callPathwayGraphlet_toremoveNodes(Pathway, tobe_removedNode_ls)
            Pathway_disruption_rates[Pathway] = disruption_score_indicator(Gpre,Gpost,
                                                                           differentiationIndicatorModel=differentiationIndicatorModel)
            print(Pathway, disruption_score_indicator(Gpre,Gpost,
                                                                       differentiationIndicatorModel=differentiationIndicatorModel))
     except:
            pass

    return Pathway_disruption_rates

"""
==========================================================================================
Section2 : Calculating Overall Mutation Effects through Collecting of Already prepared Node Importance Database. ![Node Based]!
==========================================================================================

This section aims to identify important nodes in a biological network by comparing various centrality metrics of the given network with those of random graphs.
To accomplish this, random graphs are generated with the same number of nodes and edges as the original graph, and betweenness centrality is calculated for each node in
both the original and random graphs. Subsequently, the Kolmogorov-Smirnov (KS) test is applied to compare the betweenness centrality values of the original network nodes
with those of the random graphs. By performing Monte Carlo simulations, the most optimal KS statistics are determined, which can be utilized to assess the importance of 
nodes in the biological network.

"""

from misc.Node_importance_assigner import buildFeature_analyseMatrix_wMCS_Grand, evaluateNodesImportance_byFeatures, callPathway_tomeasureNode_vitality

def notation_multiply_coeff_value(removeNode, tunedCoefficients, definedNodesImportance):
    score_summation_ls = list()
    # Every iteration indicates FeatureNCoefficient x ChosenNodes_value
    for Feature in tunedCoefficients.keys():
        try:
            score_summation_ls.append(tunedCoefficients[Feature] * definedNodesImportance[Feature][removeNode])

        except:
            print(removeNode + " is not defined into the Pathway.")
            break


    return sum(score_summation_ls)

def calculate_affectMutation_intoPathway(Pathway, tobe_removedNode_ls, already_recordedPathway_dbs=True):
     # Already calculated coefficients w/ KS Analyse for defined Pathway
     # are assigned into tunedCoefficients as Feature:Coefficient dict.

     if already_recordedPathway_dbs:
        tunedCoefficients = json.load(open("dbs/recordCoefficient_ofPathways.json","r"))[Pathway]

     else:
        tunedCoefficients, _ = buildFeature_analyseMatrix_wMCS_Grand(Pathway)
        print(tunedCoefficients)

     #Responsible the call Graph from name of Pathway and
     # find all Features of all Nodes of Graph of Pathway.
     Gbio = callPathway_tomeasureNode_vitality(Pathway)

     if already_recordedPathway_dbs:
         definedNodesImportance = json.load(open("dbs/recordGbio_definedNodeFeatures.json","r"))[Pathway]
     else:
        _, definedNodesImportance = evaluateNodesImportance_byFeatures(Gbio, definedNodeFeatures=True)

     #Every iteration as an evaluation for a Mutation.
     disrupted_score_ofAll_removedNodeS = dict()
     for removeNode in tobe_removedNode_ls:
         nodeScore_byAllFeatures = notation_multiply_coeff_value(removeNode, tunedCoefficients, definedNodesImportance)
         disrupted_score_ofAll_removedNodeS[removeNode] = nodeScore_byAllFeatures
     # Last state of disrupted_score_ofAll_removedNodeS be like -> {"MutA":0.3212, ...., "MutN":0.7543}
     return format(float(sum(list(disrupted_score_ofAll_removedNodeS.values()))),".6g")

"""
==========================================================================================
Section3 : CONTROL w/ main.
==========================================================================================
"""

'''
if __name__ == '__main__':
    tobe_removedNodes_definedPathways = {'Spinal Cord Injury': ['ACAN'], 'Endochondral Ossification': ['ACAN'],
     'Endochondral Ossification with Skeletal Dysplasias': ['ACAN'], 'Keratan sulfate biosynthesis': ['ACAN'],
     'Keratan sulfate/keratin metabolism': ['ACAN'], 'Glycosaminoglycan metabolism': ['ACAN'],
     'Metabolism of carbohydrates': ['ACAN'], 'Integrin': ['ACAN'],
     'Metabolism': ['ACAN', 'AHCY', 'MOCS1', 'PTPN13', 'CYP7A1', 'ACOX1', 'PPP1CA', 'PNPLA6', 'SLC22A5', 'SLC16A3',
                    'PTPN13', 'ACSL1', 'SBF2', 'HPGD'], 'ECM proteoglycans': ['ACAN'],
     'Extracellular matrix organization': ['ACAN', 'CAPN15', 'CAPN2'],
     'Degradation of the extracellular matrix': ['ACAN', 'CAPN15', 'CAPN2'],
     'Signaling events mediated by the Hedgehog family': ['STIL'],
     'Amyotrophic lateral sclerosis - Homo sapiens (human)': ['SETX', 'DNAH5', 'GRIA2', 'RB1CC1'],
     'Interactome of polycomb repressive complex 2 (PRC2)': ['SETX'], 'Iron metabolism in placenta': ['HEPHL1'],
     'FOXA1 transcription factor network': ['NFIC', 'EP300'],
     'Developmental Biology': ['KLK5', 'NCK1', 'DNM1', 'SH3GL2', 'EP300', 'SCN9A', 'NFASC'], 'Keratinization': ['KLK5'],
     'Formation of the cornified envelope': ['KLK5'], 'Dilated cardiomyopathy - Homo sapiens (human)': ['TTN'],
     'Hypertrophic cardiomyopathy - Homo sapiens (human)': ['TTN'], 'Striated Muscle Contraction Pathway': ['TTN'],
     'Platelet degranulation ': ['TTN', 'CHID1'], 'Platelet activation, signaling and aggregation': ['TTN', 'CHID1'],
     'Response to elevated platelet cytosolic Ca2+': ['TTN', 'CHID1'], 'Striated Muscle Contraction': ['TTN'],
     'Muscle contraction': ['TTN', 'ANXA2', 'SCN9A'], 'Hemostasis': ['TTN', 'SLC16A3', 'CHID1'],
     'Structural Pathway of Interleukin 1 (IL-1)': ['SAFB'],
     'Post-translational protein modification': ['SAFB', 'INO80D', 'PPP6C', 'CDKN2A', 'EP300', 'CDKN2A'],
     'Metabolism of proteins': ['SAFB', 'TCP1', 'LONP2', 'INO80D', 'PPP6C', 'CDKN2A', 'EP300', 'MRPS10', 'CDKN2A'],
     'SUMOylation of transcription cofactors': ['SAFB', 'EP300'],
     'SUMO E3 ligases SUMOylate target proteins': ['SAFB', 'CDKN2A', 'EP300', 'CDKN2A'],
     'Validated nuclear estrogen receptor alpha network': ['SAFB', 'EP300'],
     'SUMOylation': ['SAFB', 'CDKN2A', 'EP300', 'CDKN2A'],
     'Generic Transcription Pathway': ['ZNF668', 'PLAGL1', 'ZNF264', 'ZNF257', 'ZNF264', 'GRIA2', 'CDKN2A', 'NOTCH4',
                                       'EP300', 'CNOT7', 'CDKN2A'],
     'Gene expression (Transcription)': ['ZNF668', 'CSTF1', 'PLAGL1', 'ZNF264', 'ZNF257', 'ZNF264', 'GRIA2', 'CDKN2A',
                                         'NOTCH4', 'EP300', 'CNOT7', 'CDKN2A'],
     'RNA Polymerase II Transcription': ['ZNF668', 'CSTF1', 'PLAGL1', 'ZNF264', 'ZNF257', 'ZNF264', 'GRIA2', 'CDKN2A',
                                         'NOTCH4', 'EP300', 'CNOT7', 'CDKN2A'],
     'mRNA surveillance pathway - Homo sapiens (human)': ['CSTF1', 'PPP1CA', 'PPP2R2B'], 'mRNA Processing': ['CSTF1'],
     'Processing of Intronless Pre-mRNAs': ['CSTF1'], 'Processing of Capped Intronless Pre-mRNA': ['CSTF1'],
     'Metabolism of RNA': ['CSTF1', 'CNOT7'], 'mRNA Splicing': ['CSTF1'], 'polyadenylation of mrna': ['CSTF1'],
     'RNA Polymerase II Transcription Termination': ['CSTF1'], 'BARD1 signaling events': ['CSTF1'],
     'mRNA Splicing - Major Pathway': ['CSTF1'], 'mRNA 3,-end processing': ['CSTF1'],
     'Processing of Capped Intron-Containing Pre-mRNA': ['CSTF1'],
     'Aryl Hydrocarbon Receptor Netpath': ['PLAGL1', 'EP300', 'NF1'],
     'TP53 Regulates Transcription of Cell Cycle Genes': ['PLAGL1', 'EP300', 'CNOT7'],
     'TP53 regulates transcription of additional cell cycle genes whose exact role in the p53 pathway remain uncertain': [
         'PLAGL1', 'CNOT7'], 'Transcriptional Regulation by TP53': ['PLAGL1', 'CDKN2A', 'EP300', 'CNOT7', 'CDKN2A'],
     '16p11.2 proximal deletion syndrome': ['TCP1'],
     'Cooperation of PDCL (PhLP1) and TRiC/CCT in G-protein beta folding': ['TCP1'],
     'Association of TriC/CCT with target proteins during biosynthesis': ['TCP1', 'LONP2'],
     'Chaperonin-mediated protein folding': ['TCP1', 'LONP2'], 'Protein folding': ['TCP1', 'LONP2'],
     'Cooperation of Prefoldin and TriC/CCT  in actin and tubulin folding': ['TCP1'],
     'Prefoldin mediated transfer of substrate  to CCT/TriC': ['TCP1'],
     'Formation of tubulin folding intermediates by CCT/TriC': ['TCP1'], 'Folding of actin by CCT/TriC': ['TCP1'],
     'Pathways of neurodegeneration - multiple diseases - Homo sapiens (human)': ['DNAH5', 'GPR37', 'GRIA2', 'RB1CC1',
                                                                                  'CACNA1B', 'CAPN2'],
     'Huntington disease - Homo sapiens (human)': ['DNAH5', 'GRIA2', 'RB1CC1', 'EP300', 'CACNA1B'],
     'lissencephaly gene (lis1) in neuronal migration and development': ['DNAH5'],
     'T cell receptor signaling pathway - Homo sapiens (human)': ['NCK1'],
     'Axon guidance - Homo sapiens (human)': ['NCK1'],
     'Pathogenic Escherichia coli infection - Homo sapiens (human)': ['NCK1', 'NLRP4'],
     'ErbB signaling pathway - Homo sapiens (human)': ['NCK1'], 'Pathogenic Escherichia coli infection': ['NCK1'],
     'B Cell Receptor Signaling Pathway': ['NCK1'],
     'Brain-derived neurotrophic factor (BDNF) signaling pathway': ['NCK1', 'GRIA2'],
     'Primary focal segmental glomerulosclerosis (FSGS)': ['NCK1', 'DNM1'], 'RAC1-PAK1-p38-MMP2 Pathway': ['NCK1'],
     'Association Between Physico-Chemical Features and Toxicity Associated Pathways': ['NCK1'],
     'T-Cell antigen Receptor (TCR) pathway during Staphylococcus aureus infection': ['NCK1'],
     'VEGFA-VEGFR2 Signaling Pathway': ['NCK1', 'NOTCH4', 'PPP1CA', 'CAPN2'],
     'EGF-EGFR signaling pathway': ['NCK1', 'DNM1', 'SH3GL2'], 'Genes controlling nephrogenesis': ['NCK1'],
     'ErbB signaling pathway': ['NCK1'], 'T-cell receptor (TCR) signaling pathway': ['NCK1'],
     'TCR': ['NCK1', 'AHCY', 'ANXA2'], 'VEGFA-VEGFR2 Pathway': ['NCK1'], 'Signaling by VEGF': ['NCK1'],
     'Prolactin': ['NCK1'], 'RHO GTPase Effectors': ['NCK1', 'DIAPH3', 'PPP1R12B'],
     'Signaling by Rho GTPases': ['NCK1', 'DIAPH3', 'PTPN13', 'PPP1R12B', 'PTPN13'],
     'Signaling by Rho GTPases, Miro GTPases and RHOBTB3': ['NCK1', 'DIAPH3', 'PTPN13', 'PPP1R12B', 'PTPN13'],
     'Signaling by the B Cell Receptor (BCR)': ['NCK1'], 'Disease': ['NCK1', 'SH3GL2', 'EP300', 'CAPN2'],
     'BCR': ['NCK1'], 'Antigen activates B Cell Receptor (BCR) leading to generation of second messengers': ['NCK1'],
     'Signaling by Receptor Tyrosine Kinases': ['NCK1', 'DNM1', 'SH3GL2', 'EP300'], 'PDGF': ['NCK1'],
     'BDNF': ['NCK1', 'GRIA2'], 'DCC mediated attractive signaling': ['NCK1'], 'Netrin-1 signaling': ['NCK1'],
     'Activation of RAC1': ['NCK1'], 'Nephrin family interactions': ['NCK1'], 'Cell-Cell communication': ['NCK1'],
     'Signaling by ROBO receptors': ['NCK1'], 'Fcgamma receptor (FCGR) dependent phagocytosis': ['NCK1'],
     'Regulation of actin dynamics for phagocytic cup formation': ['NCK1'],
     'Axon guidance': ['NCK1', 'DNM1', 'SH3GL2', 'SCN9A', 'NFASC'],
     'Stabilization and expansion of the E-cadherin adherens junction': ['NCK1'], 'Insulin Pathway': ['NCK1'],
     'Netrin-mediated signaling events': ['NCK1'],
     'Signaling events mediated by focal adhesion kinase': ['NCK1', 'CAPN2'], 'EPHB forward signaling': ['NCK1'],
     'Arf6 signaling events': ['NCK1'], 'Signaling events regulated by Ret tyrosine kinase': ['NCK1'],
     'PDGFR-beta signaling pathway': ['NCK1', 'PPP2R2B'],
     'Nervous system development': ['NCK1', 'DNM1', 'SH3GL2', 'SCN9A', 'NFASC'],
     'Signaling events mediated by VEGFR1 and VEGFR2': ['NCK1', 'MYOF'],
     'TCR signaling in na&#xef;ve CD4+ T cells': ['NCK1'], 'Nephrin/Neph1 signaling in the kidney podocyte': ['NCK1'],
     'Innate Immune System': ['NCK1', 'NLRP4', 'LAIR1', 'PTPRN2', 'ANXA2', 'CFHR1', 'EP300', 'DEFB115', 'NFASC'],
     'VEGFR1 specific signals': ['NCK1'],
     'Immune System': ['NCK1', 'NLRP4', 'LAIR1', 'PTPRN2', 'DNM1', 'MRC2', 'ANXA2', 'CFHR1', 'SH3GL2', 'EP300',
                       'DEFB115', 'NFASC'], 'FCGR3A-mediated phagocytosis': ['NCK1'],
     'Leishmania phagocytosis': ['NCK1'], 'Parasite infection': ['NCK1'], 'Leishmania infection': ['NCK1'],
     'EGFR1': ['NCK1', 'DNM1', 'ANXA2', 'SH3GL2'], 'Generation of second messenger molecules': ['NCK1'],
     'y branching of actin filaments': ['NCK1'], 'TCR signaling': ['NCK1'],
     'Signaling events mediated by Hepatocyte Growth Factor Receptor (c-Met)': ['NCK1', 'SH3GL2'],
     'Angiopoietin receptor Tie2-mediated signaling': ['NCK1'], 'RHOU GTPase cycle': ['NCK1'],
     'RHOV GTPase cycle': ['NCK1'], 'RHO GTPase cycle': ['NCK1', 'DIAPH3', 'PTPN13', 'PTPN13'],
     'Infectious disease': ['NCK1', 'SH3GL2'], 'Downstream signal transduction': ['NCK1'],
     'Signaling by PDGF': ['NCK1'],
     'Signal Transduction': ['NCK1', 'GPR37', 'KDM4B', 'DIAPH3', 'PTPN13', 'DNM1', 'NOTCH4', 'PPP1R12B', 'PPP1CA',
                             'SH3GL2', 'EP300', 'NF1', 'PTPN13', 'TACR3', 'S1PR5'],
     'RHO GTPases Activate WASPs and WAVEs': ['NCK1'],
     'Adaptive Immune System': ['NCK1', 'LAIR1', 'DNM1', 'MRC2', 'SH3GL2'],
     'ECM-receptor interaction - Homo sapiens (human)': ['FREM1'],
     'Development of ureteric collection system': ['FREM1'], 'Parkinson disease - Homo sapiens (human)': ['GPR37'],
     'Parkin-Ubiquitin Proteasomal System pathway': ['GPR37'], 'Parkinson,s disease pathway': ['GPR37'],
     'GPCRs, Class A Rhodopsin-like': ['GPR37'], 'Signaling by GPCR': ['GPR37', 'PPP1CA', 'TACR3', 'S1PR5'],
     'G alpha (i) signalling events': ['GPR37', 'PPP1CA', 'S1PR5'],
     'GPCR downstream signalling': ['GPR37', 'PPP1CA', 'TACR3', 'S1PR5'],
     'role of parkin in ubiquitin-proteasomal pathway': ['GPR37'],
     'Peptide ligand-binding receptors': ['GPR37', 'TACR3'],
     'Class A/1 (Rhodopsin-like receptors)': ['GPR37', 'TACR3', 'S1PR5'],
     'GPCR ligand binding': ['GPR37', 'TACR3', 'S1PR5'],
     'Olfactory transduction - Homo sapiens (human)': ['OR6C75', 'OR6C3'],
     'Olfactory Signaling Pathway': ['OR6C75', 'OR6C3'], 'Sensory Perception': ['OR6C75', 'RGS9BP', 'OR6C3', 'RPE65'],
     'IL-18 signaling pathway': ['LONP2'], 'Protein localization': ['LONP2', 'ACOX1'],
     'Peroxisomal protein import': ['LONP2', 'ACOX1'],
     'Cysteine and methionine metabolism - Homo sapiens (human)': ['AHCY'],
     'S-Adenosylhomocysteine (SAH) Hydrolase Deficiency': ['AHCY'], 'Methionine Metabolism': ['AHCY'],
     'Methionine Adenosyltransferase Deficiency': ['AHCY'], 'Glycine N-methyltransferase Deficiency': ['AHCY'],
     'Hypermethioninemia': ['AHCY'], 'Methylenetetrahydrofolate Reductase Deficiency (MTHFRD)': ['AHCY'],
     'Betaine Metabolism': ['AHCY'],
     'Homocystinuria-megaloblastic anemia due to defect in cobalamin metabolism, cblG complementation type': ['AHCY'],
     'Selenoamino Acid Metabolism': ['AHCY'], 'Cystathionine Beta-Synthase Deficiency': ['AHCY'],
     'Folate Metabolism': ['AHCY'], 'Trans-sulfuration pathway': ['AHCY'], 'One-carbon metabolism': ['AHCY'],
     'Trans-sulfuration and one-carbon metabolism': ['AHCY'], 'Methionine De Novo and Salvage Pathway': ['AHCY'],
     'Ethanol effects on histone modifications': ['AHCY'], 'methionine degradation': ['AHCY'],
     'cysteine biosynthesis': ['AHCY'], 'Methionine Cysteine metabolism': ['AHCY'],
     'superpathway of methionine degradation': ['AHCY'], 'Methionine and cysteine metabolism': ['AHCY'],
     'Metabolism of amino acids and derivatives': ['AHCY'], 'Biological oxidations': ['AHCY', 'CYP7A1'],
     'Selenoamino acid metabolism': ['AHCY'], 'Methylation': ['AHCY'], 'Phase II - Conjugation of compounds': ['AHCY'],
     'Metabolism of ingested SeMet, Sec, MeSec into H2Se': ['AHCY'], 'Sulfur amino acid metabolism': ['AHCY'],
     'DNA Double Strand Break Response': ['KDM4B'], 'DNA Double-Strand Break Repair': ['KDM4B'],
     'DNA Repair': ['KDM4B', 'EP300'], 'HDMs demethylate histones': ['KDM4B'],
     'Chromatin modifying enzymes': ['KDM4B', 'EP300'],
     'Recruitment and ATM-mediated phosphorylation of repair and signaling proteins at DNA double strand breaks': [
         'KDM4B'], 'ESR-mediated signaling': ['KDM4B', 'EP300'], 'Chromatin organization': ['KDM4B', 'EP300'],
     'Estrogen-dependent gene expression': ['KDM4B', 'EP300'], 'Signaling by Nuclear Receptors': ['KDM4B', 'EP300'],
     'Condensation of Prometaphase Chromosomes': ['NCAPD2'],
     'akap95 role in mitosis and chromosome dynamics': ['NCAPD2'], 'Mitotic Prometaphase': ['NCAPD2', 'HAUS3'],
     'Aurora B signaling': ['NCAPD2'], 'M Phase': ['NCAPD2', 'HAUS3'],
     'Cell Cycle, Mitotic': ['NCAPD2', 'CDKN2A', 'PPP1R12B', 'HAUS3', 'EP300', 'CDKN2A'],
     'Cell Cycle': ['NCAPD2', 'PPP6C', 'CDKN2A', 'PPP1R12B', 'HAUS3', 'EP300', 'CDKN2A'],
     'Nucleotide-binding Oligomerization Domain (NOD) pathway': ['NLRP4'], 'TNFalpha': ['NLRP4', 'PPP6C', 'EP300'],
     'Regulation of innate immune responses to cytosolic DNA': ['NLRP4'],
     'Cytosolic sensors of pathogen-associated DNA ': ['NLRP4', 'EP300'],
     'mTOR signaling pathway - Homo sapiens (human)': ['FNIP1'],
     'Neutrophil degranulation': ['LAIR1', 'PTPRN2', 'ANXA2', 'NFASC'],
     'Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell': ['LAIR1'],
     'Salmonella infection - Homo sapiens (human)': ['AHNAK2', 'ANXA2'], 'Gastric Cancer Network 1': ['INO80D'],
     'UCH proteinases': ['INO80D'], 'Deubiquitination': ['INO80D', 'EP300'],
     'Folate biosynthesis - Homo sapiens (human)': ['MOCS1'], 'Molybdenum cofactor (Moco) biosynthesis': ['MOCS1'],
     'Molybdenum cofactor biosynthesis': ['MOCS1'], 'Metabolism of water-soluble vitamins and cofactors': ['MOCS1'],
     'Metabolism of vitamins and cofactors': ['MOCS1'],
     'Regulation of actin cytoskeleton - Homo sapiens (human)': ['DIAPH3', 'PPP1R12B', 'PPP1CA'],
     'Regulation of Actin Cytoskeleton': ['DIAPH3'], 'RHO GTPases Activate Formins': ['DIAPH3'],
     'RHOD GTPase cycle': ['DIAPH3'], 'ErbB1 downstream signaling': ['DIAPH3', 'CAPN2'],
     'CDC42 signaling events': ['DIAPH3'], 'RHOA GTPase cycle': ['DIAPH3'], 'RHOB GTPase cycle': ['DIAPH3'],
     'RHOC GTPase cycle': ['DIAPH3'], 'RAC1 GTPase cycle': ['DIAPH3'], 'RAC2 GTPase cycle': ['DIAPH3'],
     'RAC3 GTPase cycle': ['DIAPH3'], 'RHOG GTPase cycle': ['DIAPH3'], 'CDC42 GTPase cycle': ['DIAPH3'],
     'RHOJ GTPase cycle': ['DIAPH3'], 'RHOQ GTPase cycle': ['DIAPH3'], 'RHOF GTPase cycle': ['DIAPH3'],
     'Herpes simplex virus 1 infection - Homo sapiens (human)': ['ZNF257', 'PPP1CA'],
     'Vitamin D Receptor Pathway': ['ZNF257', 'CYP7A1', 'CDKN2A', 'CDKN2A'],
     'Apoptosis - Homo sapiens (human)': ['PTPN13', 'PTPN13', 'CAPN2'],
     'Apoptosis Modulation and Signaling': ['PTPN13', 'CDKN2A', 'PTPN13', 'CDKN2A'],
     'Ectoderm Differentiation': ['PTPN13', 'PTPN13'],
     'Metabolism of lipids': ['PTPN13', 'CYP7A1', 'ACOX1', 'PPP1CA', 'PNPLA6', 'SLC22A5', 'PTPN13', 'ACSL1', 'SBF2',
                              'HPGD'], 'fas signaling pathway (cd95)': ['PTPN13', 'PTPN13'],
     'Fc-epsilon receptor I signaling in mast cells': ['PTPN13', 'PTPN13'],
     'Ephrin B reverse signaling': ['PTPN13', 'DNM1', 'PTPN13'],
     'Synthesis of PIPs at the plasma membrane': ['PTPN13', 'PTPN13', 'SBF2'],
     'PI Metabolism': ['PTPN13', 'PNPLA6', 'PTPN13', 'SBF2'],
     'Phospholipid metabolism': ['PTPN13', 'PNPLA6', 'PTPN13', 'SBF2'], 'RND1 GTPase cycle': ['PTPN13', 'PTPN13'],
     'RND2 GTPase cycle': ['PTPN13', 'PTPN13'], 'RND3 GTPase cycle': ['PTPN13', 'PTPN13'],
     'Type I diabetes mellitus - Homo sapiens (human)': ['PTPRN2'],
     'Steroid hormone biosynthesis - Homo sapiens (human)': ['CYP7A1'],
     'Bile secretion - Homo sapiens (human)': ['CYP7A1'], 'Cholesterol metabolism - Homo sapiens (human)': ['CYP7A1'],
     'Primary bile acid biosynthesis - Homo sapiens (human)': ['CYP7A1'],
     'PPAR signaling pathway - Homo sapiens (human)': ['CYP7A1', 'ACOX1', 'ACSL1'],
     'Bile Acid Biosynthesis': ['CYP7A1'], '27-Hydroxylase Deficiency': ['CYP7A1'],
     'Congenital Bile Acid Synthesis Defect Type II': ['CYP7A1'], 'Cerebrotendinous Xanthomatosis (CTX)': ['CYP7A1'],
     'Zellweger Syndrome': ['CYP7A1'], 'Familial Hypercholanemia (FHCA)': ['CYP7A1'],
     'Congenital Bile Acid Synthesis Defect Type III': ['CYP7A1'], 'Drug Induction of Bile Acid Pathway': ['CYP7A1'],
     'Liver X receptor pathway': ['CYP7A1'], 'PPAR-alpha pathway': ['CYP7A1'],
     'Farnesoid X receptor pathway': ['CYP7A1'],
     'Nuclear Receptors Meta-Pathway': ['CYP7A1', 'MYOF', 'ACOX1', 'EP300', 'TSC22D3'],
     'Nuclear Receptors in Lipid Metabolism and Toxicity': ['CYP7A1'],
     'Angiopoietin Like Protein 8 Regulatory Pathway': ['CYP7A1'],
     'PPAR signaling pathway': ['CYP7A1', 'ACOX1', 'ACSL1'], 'Statin inhibition of cholesterol production': ['CYP7A1'],
     'Oxidation by Cytochrome P450': ['CYP7A1'], 'Bile acids synthesis and enterohepatic circulation': ['CYP7A1'],
     '7-oxo-C and 7beta-HC pathways': ['CYP7A1'], 'Metapathway biotransformation Phase I and II': ['CYP7A1'],
     'bile acid biosynthesis, neutral pathway': ['CYP7A1'], 'Phase I - Functionalization of compounds': ['CYP7A1'],
     'Bile acid biosynthesis': ['CYP7A1', 'ACOX1'], 'Steroids metabolism': ['CYP7A1'], 'Endogenous sterols': ['CYP7A1'],
     'Synthesis of bile acids and bile salts via 7alpha-hydroxycholesterol': ['CYP7A1'],
     'Synthesis of bile acids and bile salts via 27-hydroxycholesterol': ['CYP7A1'],
     'Synthesis of bile acids and bile salts': ['CYP7A1'], 'Cytochrome P450 - arranged by substrate type': ['CYP7A1'],
     'Bile acid and bile salt metabolism': ['CYP7A1'], 'Metabolism of steroids': ['CYP7A1'],
     'Bacterial invasion of epithelial cells - Homo sapiens (human)': ['DNM1'],
     'Endocrine and other factor-regulated calcium reabsorption - Homo sapiens (human)': ['DNM1'],
     'Phospholipase D signaling pathway - Homo sapiens (human)': ['DNM1'],
     'Endocytosis - Homo sapiens (human)': ['DNM1', 'SH3GL2'],
     'Synaptic vesicle cycle - Homo sapiens (human)': ['DNM1', 'CACNA1B'],
     'Synaptic Vesicle Pathway': ['DNM1', 'CACNA1B'], 'One-carbon metabolism and related pathways': ['DNM1'],
     'Gap junction trafficking and regulation': ['DNM1'], 'Retrograde neurotrophin signalling': ['DNM1', 'SH3GL2'],
     'Signaling by NTRK1 (TRKA)': ['DNM1', 'SH3GL2', 'EP300'], 'Signaling by NTRKs': ['DNM1', 'SH3GL2', 'EP300'],
     'IL8- and CXCR2-mediated signaling events': ['DNM1'], 'MHC class II antigen presentation': ['DNM1', 'SH3GL2'],
     'Formation of annular gap junctions': ['DNM1'], 'Gap junction degradation': ['DNM1'],
     'Gap junction trafficking': ['DNM1'], 'Purine metabolism': ['DNM1'], 'Recycling pathway of L1': ['DNM1', 'SH3GL2'],
     'Clathrin-mediated endocytosis': ['DNM1', 'SH3GL2'], 'L1CAM interactions': ['DNM1', 'SH3GL2', 'SCN9A', 'NFASC'],
     'Membrane Trafficking': ['DNM1', 'PPP6C', 'SH3GL2', 'SBF2'], 'PAR1-mediated thrombin signaling events': ['DNM1'],
     'Neurotrophic factor-mediated Trk receptor signaling': ['DNM1'],
     'IL8- and CXCR1-mediated signaling events': ['DNM1'], 'Internalization of ErbB1': ['DNM1', 'SH3GL2'],
     'Thromboxane A2 receptor signaling': ['DNM1'], 'CXCR3-mediated signaling events': ['DNM1'],
     'Notch signaling pathway': ['DNM1', 'NOTCH4', 'EP300'],
     'Vesicle-mediated transport': ['DNM1', 'PPP6C', 'SH3GL2', 'SBF2'],
     'gamma-aminobutyric acid receptor life cycle pathway': ['DNM1'], '-arrestins in gpcr desensitization': ['DNM1'],
     'roles of  arrestin dependent recruitment of src kinases in gpcr signaling': ['DNM1'],
     'endocytotic role of ndk phosphins and dynamin': ['DNM1'],
     'role of -arrestins in the activation and targeting of map kinases': ['DNM1'],
     'CXCR4-mediated signaling events': ['DNM1'], 'IL1': ['PPP6C'], 'COPII-mediated vesicle transport': ['PPP6C'],
     'ER to Golgi Anterograde Transport': ['PPP6C'], 'Transport to the Golgi and subsequent modification': ['PPP6C'],
     'Asparagine N-linked glycosylation': ['PPP6C'], 'Telomere Extension By Telomerase': ['PPP6C'],
     'Extension of Telomeres': ['PPP6C'], 'Telomere Maintenance': ['PPP6C'], 'Chromosome Maintenance': ['PPP6C'],
     'Cocaine addiction - Homo sapiens (human)': ['GRIA2'],
     'Amphetamine addiction - Homo sapiens (human)': ['GRIA2', 'PPP1CA'],
     'Dopaminergic synapse - Homo sapiens (human)': ['GRIA2', 'PPP1CA', 'PPP2R2B', 'CACNA1B'],
     'Glutamatergic synapse - Homo sapiens (human)': ['GRIA2'],
     'Long-term depression - Homo sapiens (human)': ['GRIA2'],
     'Nicotine addiction - Homo sapiens (human)': ['GRIA2', 'CACNA1B'],
     'Spinocerebellar ataxia - Homo sapiens (human)': ['GRIA2', 'RB1CC1'],
     'Long-term potentiation - Homo sapiens (human)': ['GRIA2', 'PPP1CA', 'EP300'],
     'Circadian entrainment - Homo sapiens (human)': ['GRIA2'],
     'cAMP signaling pathway - Homo sapiens (human)': ['GRIA2', 'ACOX1', 'PPP1CA', 'EP300'],
     'Neuroactive ligand-receptor interaction - Homo sapiens (human)': ['GRIA2', 'TACR3', 'S1PR5'],
     'Retrograde endocannabinoid signaling - Homo sapiens (human)': ['GRIA2', 'CACNA1B'],
     'Common Pathways Underlying Drug Addiction': ['GRIA2', 'PPP1CA'],
     'Fragile X Syndrome': ['GRIA2', 'CDKN2A', 'PPP1CA', 'NF1', 'CDKN2A'],
     'PKC-gamma calcium signaling pathway in ataxia': ['GRIA2'],
     'MECP2 regulates neuronal receptors and channels': ['GRIA2'], 'Transcriptional Regulation by MECP2': ['GRIA2'],
     'Trafficking of AMPA receptors': ['GRIA2'], 'Activation of AMPA receptors': ['GRIA2'],
     'Glutamate binding, activation of AMPA receptors and synaptic plasticity': ['GRIA2'],
     'Neurotransmitter receptors and postsynaptic signal transmission': ['GRIA2'],
     'Transmission across Chemical Synapses': ['GRIA2', 'CACNA1B'], 'Neuronal System': ['GRIA2', 'CACNA1B'],
     'N-cadherin signaling events': ['GRIA2'], 'Trafficking of GluR2-containing AMPA receptors': ['GRIA2'],
     'fosb gene expression and drug abuse': ['GRIA2'], 'Long-term potentiation': ['GRIA2'],
     'Post NMDA receptor activation events': ['GRIA2'],
     'Activation of NMDA receptors and postsynaptic events': ['GRIA2'],
     'Unblocking of NMDA receptors, glutamate binding and activation': ['GRIA2'],
     'Non-small cell lung cancer - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'Chronic myeloid leukemia - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'Melanoma - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'Cushing syndrome - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'Bladder cancer - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'Cell cycle - Homo sapiens (human)': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Human T-cell leukemia virus 1 infection - Homo sapiens (human)': ['CDKN2A', 'EP300', 'CDKN2A'],
     'p53 signaling pathway - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'Hepatocellular carcinoma - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'Glioma - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'MicroRNAs in cancer - Homo sapiens (human)': ['CDKN2A', 'NOTCH4', 'EP300', 'CDKN2A'],
     'Pathways in cancer - Homo sapiens (human)': ['CDKN2A', 'NOTCH4', 'EP300', 'CDKN2A'],
     'Viral carcinogenesis - Homo sapiens (human)': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Cellular senescence - Homo sapiens (human)': ['CDKN2A', 'PPP1CA', 'CAPN2', 'CDKN2A'],
     'Human cytomegalovirus infection - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'],
     'Pancreatic cancer - Homo sapiens (human)': ['CDKN2A', 'CDKN2A'], 'Cell cycle': ['CDKN2A', 'EP300', 'CDKN2A'],
     'TP53 network': ['CDKN2A', 'CDKN2A'], 'Glioblastoma signaling pathways': ['CDKN2A', 'EP300', 'NF1', 'CDKN2A'],
     'Apoptosis': ['CDKN2A', 'CDKN2A', 'CDKN2A', 'CDKN2A'], 'Bladder cancer': ['CDKN2A', 'CDKN2A'],
     'Tumor suppressor activity of SMARCB1': ['CDKN2A', 'CDKN2A'], 'Non-small cell lung cancer': ['CDKN2A', 'CDKN2A'],
     'Pancreatic adenocarcinoma pathway': ['CDKN2A', 'CDKN2A'], 'G1 to S cell cycle control': ['CDKN2A', 'CDKN2A'],
     'Gastrin signaling pathway': ['CDKN2A', 'ANXA2', 'CDKN2A'],
     'Head and Neck Squamous Cell Carcinoma': ['CDKN2A', 'CDKN2A'], 'Melanoma': ['CDKN2A', 'NF1', 'CDKN2A'],
     'Senescence and Autophagy in Cancer': ['CDKN2A', 'RB1CC1', 'CDKN2A'],
     'DNA damage response (only ATM dependent)': ['CDKN2A', 'CDKN2A'], 'Hedgehog': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Oncogene Induced Senescence': ['CDKN2A', 'CDKN2A'], 'Oxidative Stress Induced Senescence': ['CDKN2A', 'CDKN2A'],
     'Cellular Senescence': ['CDKN2A', 'CDKN2A'], 'RUNX3 regulates p14-ARF': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Transcriptional regulation by RUNX3': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Transcriptional Regulation by VENTX': ['CDKN2A', 'CDKN2A'],
     'C-MYB transcription factor network': ['CDKN2A', 'EP300', 'CDKN2A'],
     'SUMOylation of DNA damage response and repair proteins': ['CDKN2A', 'CDKN2A'],
     'SUMOylation of transcription factors': ['CDKN2A', 'CDKN2A'],
     'Cellular responses to stress': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Cyclin D associated events in G1': ['CDKN2A', 'CDKN2A'], 'G1 Phase': ['CDKN2A', 'CDKN2A'],
     'Mitotic G1 phase and G1/S transition': ['CDKN2A', 'CDKN2A'],
     'Apoptotic factor-mediated response': ['CDKN2A', 'CDKN2A'],
     'Intrinsic Pathway for Apoptosis': ['CDKN2A', 'CDKN2A'],
     'Regulation of nuclear beta catenin signaling and target gene transcription': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Coregulation of Androgen receptor activity': ['CDKN2A', 'CDKN2A'],
     'AP-1 transcription factor network': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Hypoxic and oxygen homeostasis regulation of HIF-1-alpha': ['CDKN2A', 'CDKN2A'],
     'FOXM1 transcription factor network': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Validated transcriptional targets of AP1 family members Fra1 and Fra2': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Regulation of retinoblastoma protein': ['CDKN2A', 'EP300', 'CDKN2A'], 'C-MYC pathway': ['CDKN2A', 'CDKN2A'],
     'E2F transcription factor network': ['CDKN2A', 'EP300', 'CDKN2A'], 'Programmed Cell Death': ['CDKN2A', 'CDKN2A'],
     'cyclins and cell cycle regulation': ['CDKN2A', 'CDKN2A'],
     'Senescence-Associated Secretory Phenotype (SASP)': ['CDKN2A', 'CDKN2A'],
     'Regulation of TP53 Degradation': ['CDKN2A', 'CDKN2A'],
     'Regulation of TP53 Expression and Degradation': ['CDKN2A', 'CDKN2A'],
     'Validated transcriptional targets of TAp63 isoforms': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Validated transcriptional targets of deltaNp63 isoforms': ['CDKN2A', 'CDKN2A'],
     'tumor suppressor arf inhibits ribosomal biogenesis': ['CDKN2A', 'CDKN2A'],
     'cell cycle: g1/s check point': ['CDKN2A', 'CDKN2A'],
     'Cellular responses to external stimuli': ['CDKN2A', 'EP300', 'CDKN2A'],
     'Regulation of TP53 Activity': ['CDKN2A', 'EP300', 'CDKN2A'], 'p53 pathway': ['CDKN2A', 'EP300', 'CDKN2A'],
     'ctcf: first multivalent nuclear factor': ['CDKN2A', 'CDKN2A'], 'Phagosome - Homo sapiens (human)': ['MRC2'],
     'Tuberculosis - Homo sapiens (human)': ['MRC2', 'EP300'],
     'Cross-presentation of soluble exogenous antigens (endosomes)': ['MRC2'],
     'Antigen processing-Cross presentation': ['MRC2'],
     'Class I MHC mediated antigen processing & presentation': ['MRC2'],
     'Prostaglandin Synthesis and Regulation': ['ANXA2', 'HPGD'], 'Gastrin': ['ANXA2'],
     'Smooth Muscle Contraction': ['ANXA2'], 'Breast cancer - Homo sapiens (human)': ['NOTCH4'],
     'Notch signaling pathway - Homo sapiens (human)': ['NOTCH4', 'EP300'],
     'Thyroid hormone signaling pathway - Homo sapiens (human)': ['NOTCH4', 'EP300'],
     'Human papillomavirus infection - Homo sapiens (human)': ['NOTCH4', 'EP300', 'PPP2R2B'], 'NOTCH-Ncore': ['NOTCH4'],
     'Neural Crest Differentiation': ['NOTCH4'], 'Notch Signaling': ['NOTCH4'],
     'Canonical and non-canonical Notch signaling': ['NOTCH4'],
     'Role of Osx and miRNAs in tooth development': ['NOTCH4'],
     'Epithelial to mesenchymal transition in colorectal cancer': ['NOTCH4'], 'Breast cancer pathway': ['NOTCH4'],
     'Osteoblast differentiation': ['NOTCH4'], 'Notch Signaling Pathway Netpath': ['NOTCH4', 'EP300'],
     'Pre-NOTCH Processing in the Endoplasmic Reticulum': ['NOTCH4'],
     'Pre-NOTCH Expression and Processing': ['NOTCH4', 'EP300'], 'Pre-NOTCH Processing in Golgi': ['NOTCH4'],
     'NOTCH4 Activation and Transmission of Signal to the Nucleus': ['NOTCH4'],
     'NOTCH4 Intracellular Domain Regulates Transcription': ['NOTCH4'],
     'Negative regulation of NOTCH4 signaling': ['NOTCH4'], 'Signaling by NOTCH4': ['NOTCH4'],
     'Signaling by NOTCH': ['NOTCH4', 'EP300'], 'Notch': ['NOTCH4', 'NOTCH4', 'EP300', 'EP300'],
     'Notch-HLH transcription pathway': ['NOTCH4'], 'Visual signal transduction: Rods': ['RGS9BP', 'RPE65'],
     'Visual signal transduction: Cones': ['RGS9BP', 'RPE65'],
     'Inactivation, recovery and regulation of the phototransduction cascade': ['RGS9BP'],
     'The phototransduction cascade': ['RGS9BP'], 'Visual phototransduction': ['RGS9BP', 'RPE65'],
     'Autophagy - animal - Homo sapiens (human)': ['RB1CC1'],
     'Alzheimer disease - Homo sapiens (human)': ['RB1CC1', 'CAPN2'],
     'Longevity regulating pathway - Homo sapiens (human)': ['RB1CC1'],
     'PI3K-AKT-mTOR signaling pathway and therapeutic opportunities': ['RB1CC1'],
     'Neurodegeneration with brain iron accumulation (NBIA) subtypes pathway': ['RB1CC1'],
     'Autophagy': ['RB1CC1', 'RB1CC1'], 'mTOR signaling pathway': ['RB1CC1'], 'Macroautophagy': ['RB1CC1'],
     'Aryl Hydrocarbon Receptor Pathway': ['MYOF', 'EP300'], 'Hepatitis C and Hepatocellular Carcinoma': ['MYOF'],
     'Focal adhesion - Homo sapiens (human)': ['PPP1R12B', 'PPP1CA', 'CAPN2'],
     'Oxytocin signaling pathway - Homo sapiens (human)': ['PPP1R12B', 'PPP1CA'],
     'Vascular smooth muscle contraction - Homo sapiens (human)': ['PPP1R12B', 'PPP1CA'],
     'Proteoglycans in cancer - Homo sapiens (human)': ['PPP1R12B', 'PPP1CA'],
     'Focal Adhesion': ['PPP1R12B', 'PPP1CA', 'CAPN2'], 'RHO GTPases Activate ROCKs': ['PPP1R12B'],
     'G2/M Transition': ['PPP1R12B', 'HAUS3', 'EP300'], 'Mitotic G2-G2/M phases': ['PPP1R12B', 'HAUS3', 'EP300'],
     'Regulation of PLK1 Activity at G2/M Transition': ['PPP1R12B', 'HAUS3'],
     'ccr3 signaling in eosinophils': ['PPP1R12B'], 'integrin signaling pathway': ['PPP1R12B'],
     'pkc-catalyzed phosphorylation of inhibitory phosphoprotein of myosin phosphatase': ['PPP1R12B'],
     'rho cell motility signaling pathway': ['PPP1R12B'],
     'thrombin signaling and protease-activated receptors': ['PPP1R12B'],
     'rac1 cell motility signaling pathway': ['PPP1R12B'], 'RHO GTPases activate PAKs': ['PPP1R12B'],
     'RHO GTPases activate PKNs': ['PPP1R12B'], 'Loss of Nlp from mitotic centrosomes': ['HAUS3'],
     'Centrosome maturation': ['HAUS3'], 'Anchoring of the basal body to the plasma membrane': ['HAUS3'],
     'Cilium Assembly': ['HAUS3'], 'Organelle biogenesis and maintenance': ['HAUS3'],
     'Recruitment of mitotic centrosome proteins and complexes': ['HAUS3'],
     'Loss of proteins required for interphase microtubule organization from the centrosome': ['HAUS3'],
     'Recruitment of NuMA to mitotic centrosomes': ['HAUS3'], 'AURKA Activation by TPX2': ['HAUS3'],
     'beta-Alanine metabolism - Homo sapiens (human)': ['ACOX1'],
     'Propanoate metabolism - Homo sapiens (human)': ['ACOX1'], 'Peroxisome - Homo sapiens (human)': ['ACOX1', 'ACSL1'],
     'Biosynthesis of unsaturated fatty acids - Homo sapiens (human)': ['ACOX1'],
     'Fatty acid degradation - Homo sapiens (human)': ['ACOX1', 'ACSL1'],
     'alpha-Linolenic acid metabolism - Homo sapiens (human)': ['ACOX1'], 'Estrogen Receptor Pathway': ['ACOX1'],
     'Eicosanoid metabolism via cyclooxygenases (COX)': ['ACOX1', 'HPGD'],
     'Eicosanoid metabolism via lipooxygenases (LOX)': ['ACOX1', 'HPGD'],
     'Omega-3-Omega-6 FA synthesis': ['ACOX1', 'ACSL1'], 'Tyrosine metabolism': ['ACOX1'],
     'Di-unsaturated fatty acid beta-oxidation': ['ACOX1', 'ACSL1'], 'Leukotriene metabolism': ['ACOX1', 'ACSL1'],
     'Mono-unsaturated fatty acid beta-oxidation': ['ACOX1', 'ACSL1'],
     'Omega-3 fatty acid metabolism': ['ACOX1', 'ACSL1'], 'Omega-6 fatty acid metabolism': ['ACOX1', 'ACSL1'],
     'Saturated fatty acids beta-oxidation': ['ACOX1'], 'Vitamin E metabolism': ['ACOX1'],
     'alpha-linolenic acid (ALA) metabolism': ['ACOX1', 'ACSL1'],
     'alpha-linolenic (omega3) and linoleic (omega6) acid metabolism': ['ACOX1', 'ACSL1'],
     'Beta-oxidation of very long chain fatty acids': ['ACOX1'], 'Peroxisomal lipid metabolism': ['ACOX1'],
     'Fatty acid metabolism': ['ACOX1', 'SLC22A5', 'ACSL1', 'HPGD'],
     'Trihydroxycoprostanoyl-CoA beta-oxidation': ['ACOX1'], 'Phytanic acid peroxisomal oxidation': ['ACOX1', 'ACSL1'],
     '3-oxo-10R-octadecatrienoate beta-oxidation': ['ACOX1'],
     'mechanism of gene regulation by peroxisome proliferators via ppara': ['ACOX1', 'EP300'],
     'fatty acid &beta;-oxidation (peroxisome)': ['ACOX1', 'ACSL1'], 'TYSND1 cleaves peroxisomal proteins': ['ACOX1'],
     'Complement and coagulation cascades - Homo sapiens (human)': ['CFHR1'],
     'Regulation of Complement cascade': ['CFHR1'], 'Complement cascade': ['CFHR1'],
     'Inflammatory mediator regulation of TRP channels - Homo sapiens (human)': ['PPP1CA'],
     'Platelet activation - Homo sapiens (human)': ['PPP1CA'], 'Alcoholism - Homo sapiens (human)': ['PPP1CA'],
     'Insulin resistance - Homo sapiens (human)': ['PPP1CA'], 'Oocyte meiosis - Homo sapiens (human)': ['PPP1CA'],
     'Diabetic cardiomyopathy - Homo sapiens (human)': ['PPP1CA'],
     'Adrenergic signaling in cardiomyocytes - Homo sapiens (human)': ['PPP1CA', 'PPP2R2B'],
     'Hippo signaling pathway - Homo sapiens (human)': ['PPP1CA', 'PPP2R2B'],
     'Insulin signaling pathway - Homo sapiens (human)': ['PPP1CA'],
     'cGMP-PKG signaling pathway - Homo sapiens (human)': ['PPP1CA'],
     'Excitatory Neural Signalling Through 5-HTR 4 and Serotonin': ['PPP1CA'],
     'Excitatory Neural Signalling Through 5-HTR 7 and  Serotonin ': ['PPP1CA'],
     'Excitatory Neural Signalling Through 5-HTR 6 and Serotonin ': ['PPP1CA'],
     'Intracellular Signalling Through LHCGR Receptor and Luteinizing Hormone/Choriogonadotropin': ['PPP1CA'],
     'Intracellular Signalling Through FSH Receptor and Follicle Stimulating Hormone': ['PPP1CA'],
     'Intracellular Signalling Through Histamine H2 Receptor and Histamine': ['PPP1CA'], 'TGF-Ncore': ['PPP1CA'],
     'Sphingolipid pathway': ['PPP1CA'], 'Nicotine Activity on Dopaminergic Neurons': ['PPP1CA'],
     'Hippo-Merlin Signaling Dysregulation': ['PPP1CA'], 'Netrin-UNC5B signaling pathway': ['PPP1CA'],
     'Signaling by TGF-beta Receptor Complex': ['PPP1CA'], 'Signaling by TGFB family members': ['PPP1CA'],
     'Triglyceride catabolism': ['PPP1CA'], 'Downregulation of TGF-beta receptor signaling': ['PPP1CA'],
     'TGF-beta receptor signaling activates SMADs': ['PPP1CA'], 'Triglyceride metabolism': ['PPP1CA'],
     'AndrogenReceptor': ['PPP1CA', 'EP300'], 'ALK1 signaling events': ['PPP1CA'],
     'TGF-beta receptor signaling': ['PPP1CA'], 'BMP receptor signaling': ['PPP1CA'], 'DARPP-32 events': ['PPP1CA'],
     'Opioid Signalling': ['PPP1CA'], 'regulation of ck1/cdk5 by type 1 glutamate receptors': ['PPP1CA'],
     'protein kinase a at the centrosome': ['PPP1CA'],
     'Glycerophospholipid metabolism - Homo sapiens (human)': ['PNPLA6'], 'Glycerophospholipid catabolism': ['PNPLA6'],
     'EGFR downregulation': ['SH3GL2'], 'Signaling by EGFR': ['SH3GL2'], 'Lysosome Vesicle Biogenesis': ['SH3GL2'],
     'Negative regulation of MET activity': ['SH3GL2'], 'Signaling by MET': ['SH3GL2'],
     'Golgi Associated Vesicle Biogenesis': ['SH3GL2'], 'trans-Golgi Network Vesicle Budding': ['SH3GL2'],
     'Cargo recognition for clathrin-mediated endocytosis': ['SH3GL2'], 'RhoA signaling pathway': ['SH3GL2'],
     'InlB-mediated entry of Listeria monocytogenes into host cell': ['SH3GL2'],
     'Listeria monocytogenes entry into host cells': ['SH3GL2'],
     'TGF-beta signaling pathway - Homo sapiens (human)': ['EP300'],
     'Adherens junction - Homo sapiens (human)': ['EP300'],
     'Growth hormone synthesis, secretion and action - Homo sapiens (human)': ['EP300'],
     'HIF-1 signaling pathway - Homo sapiens (human)': ['EP300'],
     'JAK-STAT signaling pathway - Homo sapiens (human)': ['EP300'],
     'Renal cell carcinoma - Homo sapiens (human)': ['EP300'],
     'Wnt signaling pathway - Homo sapiens (human)': ['EP300'],
     'FoxO signaling pathway - Homo sapiens (human)': ['EP300'], 'Prostate cancer - Homo sapiens (human)': ['EP300'],
     'Glucagon signaling pathway - Homo sapiens (human)': ['EP300'], 'Melanogenesis - Homo sapiens (human)': ['EP300'],
     'Influenza A - Homo sapiens (human)': ['EP300'], 'Hepatitis B - Homo sapiens (human)': ['EP300'],
     'Kaposi sarcoma-associated herpesvirus infection - Homo sapiens (human)': ['EP300'],
     'Androgen receptor signaling pathway': ['EP300'], 'Energy Metabolism': ['EP300'],
     'Integrated breast cancer pathway': ['EP300'], 'Kit receptor signaling pathway': ['EP300'],
     'Initiation of transcription and translation elongation at the HIV-1 LTR': ['EP300'],
     'Pathways affected in adenoid cystic carcinoma': ['EP300'],
     'Hematopoietic Stem Cell Gene Regulation by GABP alpha-beta Complex': ['EP300'],
     'TGF-beta Signaling Pathway': ['EP300'], 'IL-4 signaling pathway': ['EP300'],
     'Wnt signaling pathway and pluripotency': ['EP300', 'PPP2R2B'], 'Prion disease pathway': ['EP300'],
     'Pathways in clear cell renal cell carcinoma': ['EP300'], 'Ebola Virus Pathway on Host': ['EP300'],
     'Type 2 papillary renal cell carcinoma': ['EP300'], 'Hepatitis B infection': ['EP300'],
     'Male infertility': ['EP300'], 'TGF-beta receptor signaling in skeletal dysplasias': ['EP300'],
     'CAMKK2 Pathway': ['EP300'], 'SARS-CoV-2 innate immunity evasion and cell-specific immune response': ['EP300'],
     'TGF-beta Receptor Signaling': ['EP300'], 'Sudden Infant Death Syndrome (SIDS) Susceptibility Pathways': ['EP300'],
     'IL6': ['EP300'], 'Attenuation phase': ['EP300'], 'Oncostatin_M': ['EP300'],
     'RUNX1 interacts with co-factors whose precise effect on RUNX1 targets is not known': ['EP300'],
     'Transcriptional regulation by RUNX1': ['EP300'],
     'Formation of the beta-catenin:TCF transactivating complex': ['EP300'],
     'TCF dependent signaling in response to WNT': ['EP300'], 'Cellular response to heat stress': ['EP300'],
     'B-WICH complex positively regulates rRNA expression': ['EP300'],
     'Positive epigenetic regulation of rRNA expression': ['EP300'],
     'Transcriptional regulation of granulopoiesis': ['EP300'], 'KitReceptor': ['EP300'],
     'Epigenetic regulation of gene expression': ['EP300'], 'HSF1-dependent transactivation': ['EP300'],
     'RUNX3 regulates NOTCH signaling': ['EP300'], 'Constitutive Signaling by NOTCH1 PEST Domain Mutants': ['EP300'],
     'Signaling by NOTCH1 PEST Domain Mutants in Cancer': ['EP300'],
     'Signaling by NOTCH1 HD+PEST Domain Mutants in Cancer': ['EP300'], 'Signaling by NOTCH1 in Cancer': ['EP300'],
     'TRAF3-dependent IRF activation pathway': ['EP300'], 'NGF-stimulated transcription': ['EP300'],
     'Nuclear Events (kinase and transcription factor activation)': ['EP300'],
     'NOTCH2 intracellular domain regulates transcription': ['EP300'], 'Signaling by NOTCH2': ['EP300'],
     'Formation of TC-NER Pre-Incision Complex': ['EP300'],
     'Constitutive Signaling by NOTCH1 HD+PEST Domain Mutants': ['EP300'],
     'Pre-NOTCH Transcription and Translation': ['EP300'], 'Signaling by NOTCH1': ['EP300'],
     'Dual incision in TC-NER': ['EP300'], 'role of erbb2 in signal transduction and oncology': ['EP300'],
     'Gap-filling DNA repair synthesis and ligation in TC-NER': ['EP300'],
     'Transcription-Coupled Nucleotide Excision Repair (TC-NER)': ['EP300'], 'Nucleotide Excision Repair': ['EP300'],
     'TGF_beta_Receptor': ['EP300', 'SAMD3'], 'NOTCH3 Intracellular Domain Regulates Transcription': ['EP300'],
     'Signaling by NOTCH3': ['EP300'], 'Metalloprotease DUBs': ['EP300'],
     'Regulation of FOXO transcriptional activity by acetylation': ['EP300'],
     'Diseases of signal transduction by growth factor receptors and second messengers': ['EP300'],
     'FOXO-mediated transcription of cell death genes': ['EP300'], 'CD209 (DC-SIGN) signaling': ['EP300'],
     'C-type lectin receptors (CLRs)': ['EP300'], 'FOXO-mediated transcription': ['EP300'], 'LIF signaling': ['EP300'],
     'IL4': ['EP300'], 'TGF-beta super family signaling pathway canonical': ['EP300'],
     'RORA activates gene expression': ['EP300'], 'TRAF6 mediated IRF7 activation': ['EP300'],
     'DDX58/IFIH1-mediated induction of interferon-alpha/beta': ['EP300'], 'BMP signaling Dro': ['EP300'],
     'LRR FLII-interacting protein 1 (LRRFIP1) activates type I IFN production': ['EP300'],
     'BMP2 signaling TGF-beta MV': ['EP300'], 'hypoxia and p53 in the cardiovascular system': ['EP300'],
     'Wnt Canonical': ['EP300', 'DKK3'], 'Notch-mediated HES/HEY network': ['EP300'],
     'Glucocorticoid receptor regulatory network': ['EP300'],
     'Validated targets of C-MYC transcriptional repression': ['EP300'],
     'HIF-2-alpha transcription factor network': ['EP300'], 'HIF-1-alpha transcription factor network': ['EP300'],
     'Regulation of Androgen receptor activity': ['EP300'], 'p73 transcription factor network': ['EP300'],
     'Signaling events mediated by HDAC Class III': ['EP300'], 'ATF-2 transcription factor network': ['EP300', 'NF1'],
     'Role of Calcineurin-dependent NFAT signaling in lymphocytes': ['EP300'], 'IFN-gamma pathway': ['EP300'],
     'Direct p53 effectors': ['EP300'],
     'TP53 Regulates Transcription of Genes Involved in G2 Cell Cycle Arrest': ['EP300'],
     'role of mef2d in t-cell apoptosis': ['EP300', 'CAPN2'], 'Polo-like kinase mediated events': ['EP300'],
     'FoxO family signaling': ['EP300'], 'Transcriptional regulation of white adipocyte differentiation': ['EP300'],
     'nfkb activation by nontypeable hemophilus influenzae': ['EP300'],
     'hypoxia-inducible factor in the cardivascular system': ['EP300'],
     'Regulation of TP53 Activity through Methylation': ['EP300'],
     'Retinoic acid receptors-mediated signaling': ['EP300'], 'cell cycle: g2/m checkpoint': ['EP300'],
     'erythropoietin mediated neuroprotection through nf-kb': ['EP300'],
     'melanocyte development and pigmentation pathway': ['EP300'], 'il-7 signal transduction': ['EP300'],
     'pelp1 modulation of estrogen receptor activity': ['EP300'],
     'acetylation and deacetylation of rela in nucleus': ['EP300'],
     'transcription regulation by methyltransferase of carm1': ['EP300'], 'tgf beta signaling pathway': ['EP300'],
     'multi-step regulation of transcription by pitx2': ['EP300'],
     'carm1 and regulation of the estrogen receptor': ['EP300'], 'Wnt Mammals': ['EP300', 'DKK3'],
     'PI5P Regulates TP53 Acetylation': ['EP300'], 'Regulation of TP53 Activity through Acetylation': ['EP300'],
     'Activation of the TFAP2 (AP-2) family of transcription factors': ['EP300'],
     'Signaling events mediated by HDAC Class I': ['EP300'],
     'RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function': ['EP300'],
     'Validated targets of C-MYC transcriptional activation': ['EP300'],
     'Regulation of nuclear SMAD2/3 signaling': ['EP300'],
     'Transcriptional regulation by the AP-2 (TFAP2) family of transcription factors': ['EP300'],
     'NOTCH1 Intracellular Domain Regulates Transcription': ['EP300'], 'NR1H2 and NR1H3-mediated signaling': ['EP300'],
     'NR1H3 & NR1H2 regulate gene expression linked to cholesterol transport and efflux': ['EP300'],
     'Signaling by WNT': ['EP300'], 'Circadian Clock': ['EP300'], 'HATs acetylate histones': ['EP300'],
     'Choline metabolism in cancer - Homo sapiens (human)': ['SLC22A5'],
     'Glycine, serine, alanine and threonine metabolism': ['SLC22A5', 'SLC16A3'], 'Lipoate metabolism': ['SLC22A5'],
     'Carnitine metabolism': ['SLC22A5'], 'Organic cation transport': ['SLC22A5'],
     'Glycerophospholipid metabolism': ['SLC22A5'], 'SLC-mediated transmembrane transport': ['SLC22A5', 'SLC16A3'],
     'Transport of small molecules': ['SLC22A5', 'SLC16A3', 'CLCN1', 'TSC22D3'],
     'Organic cation/anion/zwitterion transport': ['SLC22A5'],
     'Transport of bile salts and organic acids, metal ions and amine compounds': ['SLC22A5', 'SLC16A3'],
     'PI3K-Akt signaling pathway - Homo sapiens (human)': ['PPP2R2B'],
     'Tight junction - Homo sapiens (human)': ['PPP2R2B'], 'AMPK signaling pathway - Homo sapiens (human)': ['PPP2R2B'],
     'Sphingolipid signaling pathway - Homo sapiens (human)': ['PPP2R2B', 'S1PR5'],
     'Chagas disease - Homo sapiens (human)': ['PPP2R2B'], 'Hepatitis C - Homo sapiens (human)': ['PPP2R2B'],
     'Focal Adhesion-PI3K-Akt-mTOR-signaling pathway': ['PPP2R2B'], 'PI3K-Akt signaling pathway': ['PPP2R2B'],
     'Glycogen Synthesis and Degradation': ['PPP2R2B'], 'ATR signaling pathway': ['PPP2R2B'],
     'GPCR Dopamine D1like receptor': ['PPP2R2B'], 'insulin': ['PPP2R2B'], 'insulin Mam': ['PPP2R2B'],
     'Central carbon metabolism in cancer - Homo sapiens (human)': ['SLC16A3'],
     'Metabolic reprogramming in colon cancer': ['SLC16A3'],
     'Pyruvate metabolism and Citric Acid (TCA) cycle': ['SLC16A3'], 'Glycolysis and Gluconeogenesis': ['SLC16A3'],
     'Butanoate metabolism': ['SLC16A3'], 'Proton-coupled monocarboxylate transport': ['SLC16A3'],
     'Basigin interactions': ['SLC16A3'], 'Cell surface interactions at the vascular wall': ['SLC16A3'],
     'Pyruvate metabolism': ['SLC16A3'], 'The citric acid (TCA) cycle and respiratory electron transport': ['SLC16A3'],
     'Ras signaling pathway - Homo sapiens (human)': ['NF1'],
     'MAPK signaling pathway - Homo sapiens (human)': ['NF1', 'CACNA1B'], 'Pilocytic astrocytoma': ['NF1'],
     'MECP2 and Associated Rett Syndrome': ['NF1'], 'MAPK Signaling Pathway': ['NF1', 'CACNA1B'],
     'Ras signaling': ['NF1'], 'Transcription co-factors SKI and SKIL protein partners': ['NF1'],
     'Synaptic signaling pathways associated with autism spectrum disorder': ['NF1'],
     'EGFR Tyrosine Kinase Inhibitor Resistance': ['NF1'],
     'chromatin remodeling by hswi/snf atp-dependent complexes': ['NF1'], 'Regulation of RAS by GAPs': ['NF1'],
     'RAF/MAP kinase cascade': ['NF1'], 'MAPK1/MAPK3 signaling': ['NF1'], 'MAPK family signaling cascades': ['NF1'],
     'Regulation of Ras family activation': ['NF1'], 'FOXA2 and FOXA3 transcription factor networks': ['NF1'],
     'Syndecan-2-mediated signaling events': ['NF1'], 'Ion channel transport': ['CLCN1', 'TSC22D3'],
     'Stimuli-sensing channels': ['CLCN1', 'TSC22D3'], 'Prion disease - Homo sapiens (human)': ['CACNA1B'],
     'Morphine addiction - Homo sapiens (human)': ['CACNA1B'],
     'Serotonergic synapse - Homo sapiens (human)': ['CACNA1B'],
     'Type II diabetes mellitus - Homo sapiens (human)': ['CACNA1B'],
     'GABAergic synapse - Homo sapiens (human)': ['CACNA1B'],
     'Calcium signaling pathway - Homo sapiens (human)': ['CACNA1B', 'TACR3'],
     'Cholinergic synapse - Homo sapiens (human)': ['CACNA1B'], 'Calcium Regulation in the Cardiac Cell': ['CACNA1B'],
     'Presynaptic depolarization and calcium channel opening': ['CACNA1B'],
     'Protein processing in endoplasmic reticulum - Homo sapiens (human)': ['CAPN2'],
     'Necroptosis - Homo sapiens (human)': ['CAPN2'], 'Shigellosis - Homo sapiens (human)': ['CAPN2'],
     'Integrin-mediated Cell Adhesion': ['CAPN2'], 'Alzheimer,s disease': ['CAPN2'],
     'Deregulated CDK5 triggers multiple neurodegenerative pathways in Alzheimer,s disease models': ['CAPN2'],
     'Neurodegenerative Diseases': ['CAPN2'], 'Diseases of programmed cell death': ['CAPN2'],
     'mcalpain and friends in cell motility': ['CAPN2'],
     'Regulation of Wnt-B-catenin Signaling by Small Molecule Compounds': ['DKK3'],
     'Adipocytokine signaling pathway - Homo sapiens (human)': ['ACSL1'],
     'Ferroptosis - Homo sapiens (human)': ['ACSL1'], 'Thermogenesis - Homo sapiens (human)': ['ACSL1'],
     'Fatty acid biosynthesis - Homo sapiens (human)': ['ACSL1'],
     'Long chain acyl-CoA dehydrogenase deficiency (LCAD)': ['ACSL1'], 'Trifunctional protein deficiency': ['ACSL1'],
     'Carnitine palmitoyl transferase deficiency (II)': ['ACSL1'],
     'Very-long-chain acyl coa dehydrogenase deficiency (VLCAD)': ['ACSL1'],
     'Medium chain acyl-coa dehydrogenase deficiency (MCAD)': ['ACSL1'],
     'Beta Oxidation of Very Long Chain Fatty Acids': ['ACSL1'],
     'Short Chain Acyl CoA Dehydrogenase Deficiency (SCAD Deficiency)': ['ACSL1'], 'Fatty acid Metabolism': ['ACSL1'],
     'Oxidation of Branched Chain Fatty Acids': ['ACSL1'], 'Glutaric Aciduria Type I': ['ACSL1'],
     'Ethylmalonic Encephalopathy': ['ACSL1'],
     'Mitochondrial Beta-Oxidation of Long Chain Saturated Fatty Acids': ['ACSL1'],
     'Adrenoleukodystrophy, X-linked': ['ACSL1'], 'Carnitine-acylcarnitine translocase deficiency': ['ACSL1'],
     'Carnitine palmitoyl transferase deficiency (I)': ['ACSL1'], 'Fatty acid beta-oxidation': ['ACSL1'],
     'Fatty Acid Biosynthesis': ['ACSL1'], 'Mitochondrial LC-Fatty Acid Beta-Oxidation': ['ACSL1'],
     'Ferroptosis': ['ACSL1'], 'Thermogenesis': ['ACSL1'],
     'Cholesterol metabolism (includes both Bloch and Kandutsch-Russell pathways)': ['ACSL1'],
     'Omega-9 FA synthesis': ['ACSL1'], 'Fatty acid transporters': ['ACSL1'], 'fatty acid &alpha;-oxidation': ['ACSL1'],
     'Synthesis of very long-chain fatty acyl-CoAs': ['ACSL1'], 'Fatty acyl-CoA biosynthesis': ['ACSL1'],
     'Linoleic acid (LA) metabolism': ['ACSL1'], 'stearate biosynthesis': ['ACSL1'],
     'eicosapentaenoate biosynthesis': ['ACSL1'], 'fatty acid &beta;-oxidation': ['ACSL1'],
     '&gamma;-linolenate biosynthesis': ['ACSL1'], 'fatty acid activation': ['ACSL1'],
     'Ribosome - Homo sapiens (human)': ['MRPS10'], 'Mitochondrial translation initiation': ['MRPS10'],
     'Translation': ['MRPS10'], 'Mitochondrial translation elongation': ['MRPS10'],
     'Mitochondrial translation termination': ['MRPS10'], 'Mitochondrial translation': ['MRPS10'],
     'Peptide GPCRs': ['TACR3'], 'G alpha (q) signalling events': ['TACR3'],
     'Tachykinin receptors bind tachykinins': ['TACR3'],
     'Intracellular trafficking proteins involved in CMT neuropathy': ['SBF2'],
     'RAB GEFs exchange GTP for GDP on RABs': ['SBF2'], 'Rab regulation of trafficking': ['SBF2'],
     'Taste transduction - Homo sapiens (human)': ['SCN9A'], 'Interaction between L1 and Ankyrins': ['SCN9A', 'NFASC'],
     'Phase 0 - rapid depolarisation': ['SCN9A'], 'Cardiac conduction': ['SCN9A'],
     'Retinol metabolism - Homo sapiens (human)': ['RPE65'], 'Vitamin A Deficiency': ['RPE65'],
     'Retinol Metabolism': ['RPE65'], 'Vitamin A and carotenoid metabolism': ['RPE65'],
     'the visual cycle I (vertebrates)': ['RPE65'], 'The canonical retinoid cycle in rods (twilight vision)': ['RPE65'],
     'Glucocorticoid Receptor Pathway': ['TSC22D3'], 'RNA degradation - Homo sapiens (human)': ['CNOT7'],
     'Deadenylation of mRNA': ['CNOT7'], 'Deadenylation-dependent mRNA decay': ['CNOT7'],
     'Signal Transduction of S1P Receptor': ['S1PR5'], 'S1P4 pathway': ['S1PR5'], 'S1P5 pathway': ['S1PR5'],
     'Sphingosine 1-phosphate (S1P) pathway': ['S1PR5'], 'Lysosphingolipid and LPA receptors': ['S1PR5'],
     'Antimicrobial peptides': ['DEFB115'], 'Beta defensins': ['DEFB115'], 'Defensins': ['DEFB115'],
     'Cell adhesion molecules - Homo sapiens (human)': ['NFASC'], 'Neurofascin interactions': ['NFASC'],
     'Transcriptional misregulation in cancer - Homo sapiens (human)': ['HPGD'],
     'Biosynthesis of EPA-derived SPMs': ['HPGD'], 'Biosynthesis of E-series 18(S)-resolvins': ['HPGD'],
     'Biosynthesis of specialized proresolving mediators (SPMs)': ['HPGD'],
     'Synthesis of Prostaglandins (PG) and Thromboxanes (TX)': ['HPGD'],
     'Arachidonic acid metabolism': ['HPGD', 'HPGD'], 'Prostaglandin formation from arachidonate': ['HPGD'],
     'Prostaglandin formation from dihomo gama-linoleic acid': ['HPGD'], 'Synthesis of Lipoxins (LX)': ['HPGD'],
     'Biosynthesis of D-series resolvins': ['HPGD'], 'Biosynthesis of DHA-derived SPMs': ['HPGD']}

    #### Choose Process ####
    ScoreBaseModel = 2 #1 corresponds to Graph based disruption score;
                       #2 corresponds to Node based disruption score.

    if ScoreBaseModel == 1:
       print(main_disruption_rateGraph(tobe_removedNodes_definedPathways))

    else:
       Pathway = "TCR signaling in na&#xef;ve CD4+ T cells"
       #tobe_removedNode_ls = ['ACAN','CAPN15','CAPN2'] #Note -> There are Mutations defined into the pathway before,
                                               # although not present in the .gexf Graph (for this example CAPN15) SOLVE IT !

       # or
       # Pathway = "Extracellular matrix organization"
       tobe_removedNode_ls = tobe_removedNodes_definedPathways[Pathway]

       #Just Involved Pathway and Removed (Mutated) Nodes
       print(Pathway ,calculate_affectMutation_intoPathway(Pathway, tobe_removedNode_ls))
'''