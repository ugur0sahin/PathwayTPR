import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import networkx as nx
from scipy.linalg import expm
from scipy.stats import ks_2samp


def run_KS_test(columnVector_forFeature_Gbio, columnVector_forFeature_Grand):
	GbioValues,GrandValues = list(columnVector_forFeature_Gbio), list(columnVector_forFeature_Grand)
	ks_statistic, p_value = ks_2samp(GbioValues, GrandValues)

	return ks_statistic, p_value

def callPathway_tomeasureNode_vitality(Pathway):
	# It will give outputs as preGraph and removed(processed)Graph
	Graphlet_flname = str(Pathway).replace(" ", "_") + ".gexf"  #Sometimes There is no any Pathway cause the bugs!!
	return nx.read_gexf("../dbs/Graphlets/" + Graphlet_flname)  # call Pathway and define into NetworkX Graph

def evaluateNodesImportance_byFeatures(G, NonObjectedFeatures = tuple(), definedNodeFeatures = False):
	#NonObjectedFeatures is a tuple which contain not consider indicative Features
	#If there is one or more indicator in this parameter it won't be consider in the Importance scoring of Nodes
	#! Bug ! -> Inf results are problematic.
	ObjectedFeatures_return, definedNodeFeature_values = dict(), dict()
	if 'Degree Centrality' not in NonObjectedFeatures:
		degree_centralityNodes =  nx.degree_centrality(G)
		ObjectedFeatures_return["Degree Centrality"] = list(degree_centralityNodes.values())

		if definedNodeFeatures is True:
			definedNodeFeature_values["Degree Centrality"] = degree_centralityNodes

	if 'Closeness Centrality' not in NonObjectedFeatures:
		closeness_centralityNodes = nx.closeness_centrality(G)
		ObjectedFeatures_return['Closeness Centrality'] = list(closeness_centralityNodes.values())

		if definedNodeFeatures is True:
			definedNodeFeature_values["Closeness Centrality"] = closeness_centralityNodes

	if 'Betweenness Centrality' not in NonObjectedFeatures:
		betweenness_centrality = nx.betweenness_centrality(G)
		ObjectedFeatures_return['Betweenness Centrality'] = list(betweenness_centrality.values())

		if definedNodeFeatures is True:
			definedNodeFeature_values["Betweenness Centrality"] = degree_centralityNodes

	if 'Eigenvector Centrality' not in NonObjectedFeatures:
		eigenvector_centrality = nx.eigenvector_centrality(G)
		ObjectedFeatures_return['Eigenvector Centrality'] = list(eigenvector_centrality.values())

		if definedNodeFeatures is True:
			definedNodeFeature_values["Eigenvector Centrality"] = degree_centralityNodes

	if 'Eccentricity Centrality' not in NonObjectedFeatures:
		try:
			eccentricity = nx.eccentricity(G)
			max_eccentricity = max(eccentricity.values())
			eccentricity_centrality = {node: 1 / max_eccentricity for node in G.nodes}
			ObjectedFeatures_return['Eccentricity Centrality'] = list(eccentricity_centrality.values())

			if definedNodeFeatures is True:
				definedNodeFeature_values["Eccentricity Centrality"] = degree_centralityNodes

		except:
			pass

	if 'Subgraph Centrality' not in NonObjectedFeatures:
		adj_matrix = nx.adjacency_matrix(G).todense().astype(np.float64)
		exp_adj_matrix = expm(adj_matrix)
		subgraph_centrality = {n: exp_adj_matrix[i, i] for i, n in enumerate(G.nodes)}
		ObjectedFeatures_return['Subgraph Centrality'] = list(subgraph_centrality.values())

		if definedNodeFeatures is True:
			definedNodeFeature_values["Subgraph Centrality"] = degree_centralityNodes

	return ObjectedFeatures_return, definedNodeFeature_values


def buildFeature_analyseMatrix_wMCS_Grand(Pathway,numberSimulation=10):
		Gbio,ks_statistics_forFeatures = callPathway_tomeasureNode_vitality(Pathway), list()
		"""Gbio.remove_nodes_from(list(nx.isolates(Gbio)))"""
		FeaturesGbio_ofNodes, _ = evaluateNodesImportance_byFeatures(Gbio)

		for _ in range(numberSimulation):
			# Generate a random barabasi Graph with the same number of Nodes and Edges as Gbio
			average_degree = sum([degree for node, degree in Gbio.degree()]) / len(Gbio.nodes())
			Grand = nx.barabasi_albert_graph(Gbio.number_of_nodes(), int(average_degree))
			"""Grand.remove_nodes_from(list(nx.isolates(Grand)))"""
			# Calculate overall Nodes' Features (like betweenness centrality, etc) for Gbio and Grand, output is Feature defined lists.
			FeaturesGrand_ofNodes, _ = evaluateNodesImportance_byFeatures(Grand)

			# Run KS test for every feature retrieved from the FeaturesGbio_ofNodes,
			# FeaturesGrand_ofNodes which is the Nodes state dictionary for one Gbio,Grand Graphs.

			stepFeaturesKS_tmp = list()
			for Feature in FeaturesGbio_ofNodes.keys():
				ks_stat_forFeature, _ = run_KS_test(list(FeaturesGbio_ofNodes[Feature]), list(FeaturesGrand_ofNodes[Feature]))
				stepFeaturesKS_tmp.append(ks_stat_forFeature)
			ks_statistics_forFeatures.append(stepFeaturesKS_tmp)
		# The last form of ks_statistics_forFeatures is [1stIteration [ x1, x2 ,x3 ,x4 .... ,xN],
		# 									 2ndIteration [ x1, x2 ,x3 ,x4 .... ,xN],
		#									 ...
		#									 Nth Iteration[ x1, x2 ,x3 ,x4 .... ,xN]]
		# Every column corresponds a Feature's simulation results. Mean of them gives the optimum distance.


		# Calculate the average KS statistic from iteration list as last output of function.
		MC_simulationMatrix, KS_coefficients = pd.DataFrame(np.array(ks_statistics_forFeatures), columns=list(FeaturesGbio_ofNodes.keys())), dict()
		for Feature in list(FeaturesGbio_ofNodes.keys()):
			KS_coefficients[Feature] = np.mean(MC_simulationMatrix[Feature])

		return KS_coefficients


if __name__ != '__main__':
    #print(buildFeature_analyseMatrix_wMCS_Grand('Generic Transcription Pathway'))
	pass