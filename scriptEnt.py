# coding: utf-8
# Powered by Python 3.5

# To cancel the modifications performed by the script
# on the current graph, click on the undo button.

# Some useful keyboards shortcuts : 
#   * Ctrl + D : comment selected lines.
#   * Ctrl + Shift + D  : uncomment selected lines.
#   * Ctrl + I : indent selected lines.
#   * Ctrl + Shift + I  : unindent selected lines.
#   * Ctrl + Return  : run script.
#   * Ctrl + F  : find selected text.
#   * Ctrl + R  : replace selected text.
#   * Ctrl + Space  : show auto-completion dialog.

from tulip import tlp
import random
import numpy as np
from scipy import spatial
from math import sqrt
from numpy import median, percentile


# The updateVisualization(centerViews = True) function can be called
# during script execution to update the opened views

# The pauseScript() function can be called to pause the script execution.
# To resume the script execution, you will have to click on the "Run script " button.

# The runGraphScript(scriptFile, graph) function can be called to launch
# another edited script on a tlp.Graph object.
# The scriptFile parameter defines the script name to call (in the form [a-zA-Z0-9_]+.py)

# The main(graph) function must be defined 
# to run the script on the current graph

def label(graph, Locus, viewLabel):  
  for n in graph.getNodes():
    #print(n)
    #print( Locus[n])
    viewLabel[n] = Locus[n]
    
def setSize(graphe, viewSize):
  
  baseSize = tlp.Size(2,0.5,1)

  for n in graph.getNodes():
    viewSize[n] = baseSize #* (graph.deg(n) + 1)
    
def changeEdge(graph, Negative, Positive,viewColor):
  
  blue = tlp.Color(58,95,205)
  green = tlp.Color(69,139,0)
  
  for e in graph.getEdges():
    if Negative[e] is True:
      viewColor[e] = blue
    else:
      viewColor[e] = green
     
def NodePosition(graph, viewLayout):
  # Apply an FM^3 layout on it
  fm3pParams = tlp.getDefaultPluginParameters("FM^3 (OGDF)", graph)
#  fm3pParams["Unit edge length"] = 200
  graph.applyLayoutAlgorithm("FM^3 (OGDF)", viewLayout, fm3pParams)
 
  #for n in graph.getNodes():
   # x = random.random() * size
   # y = random.random() * size
   # viewLayout[n] = tlp.Coord(x, y, 0)


def geneClustering(graph, all_tp, similarityFct, seuil):
  
  # creation ou reinitialisation du graph de correlation 
  graph_correlation_name = "correlation_{}".format(similarityFct.__name__)
  correlation_graph = reinitializeSubgraph(graph, graph_correlation_name)
  
  
  # Supression de toutes les arretes
  print("Supression des aretes du graphe ",graph_correlation_name)
  correlation_graph.delEdges(correlation_graph.getEdges())
  
  # creation de la doubleProperty "correlation" permettant d'enregistrer 
  # la correlation entre chaque gene au niveau des arretes
  
  
  similarity_prop = computeCorrelation(correlation_graph, all_tp, similarityFct, seuil)
  print(similarity_prop)
  
  
  clusteringMCL(correlation_graph, similarity_prop )
  
  return correlation_graph



def computeSimilarity(similarity_graph, all_tp, similarityFct):
  """
  parametres:
    similarity_graph: graph qui va contenir toutes les aretes entre chaque neouds portant une valeur de similarité pour leurs niveau d'exression
    all_tp :dictionnaire des propriété des niveaux d'expression à chaque temps
    similarityFct: function calculant une similarité entre deux liste de valeur. Peut correspondre à la fonction cosineSimilarity, jackknife ou pearson 
  """

#  correlation_properties = {}
#  for s in seuils:
#    correlation_properties[s] =  graph.getDoubleProperty("correlation_{}".format(s))
     
  correlation_value = similarity_graph.getDoubleProperty("correlation_{}".format(similarityFct.__name__))
  
  
  nannodes = []
  treated_node = []
  cpt = 0.0
  step = 0
  nodes = list(similarity_graph.getNodes())
  nb_node = len(nodes)
  print("nb node", nb_node)
  for nodeA in similarity_graph.getNodes():
    try:
      nodes.remove(nodeA)
    except ValueError:
      continue
   
    exprA = []
    for tp in all_tp:
      exprA.append(tp[nodeA])
    
    if sum(exprA) == 0:
      nannodes.append(nodeA)
      similarity_graph.delNode(nodeA)
      continue 
    


    for nodeB in nodes:
      exprB = []

      for tp in all_tp:
        exprB.append(tp[nodeB])
    
      if sum(exprB) == 0:
        nannodes.append(nodeB)
        similarity_graph.delNode(nodeB)
        nodes.remove(nodeB)
        continue
       
       
      result = similarityFct(exprA, exprB)
#      print(result)
#      print('cosine scipy', 1-spatial.distance.cosine(exprA, exprB))
#      result = np.corrcoef(exprA, exprB)[0][1]
#      pearson_corr = pearson(exprA, exprB)
#      cosine = cosine_similarity(exprA, exprB)
#      jk = jackknife(exprA, exprB)
#      print('cosine', cosine)
#      print('cosine scipy', spatial.distance.cosine(exprA, exprB)) 
#      print('pearson scipy', result)
#      print('pearson homemade',pearson_corr)
#      print('jk', jk)
      
      #Add edge between node A and node B
      edge = similarity_graph.addEdge(nodeA, nodeB)
      correlation_value[edge] = result
      
    
#    print(float(cpt)/nb_node)
    if cpt/nb_node >= step:
      print("{}% de la mesure de similarity".format(int(step)*100))
      step += 0.05
    cpt += 1
      
#      if result > seuil:
#        edge = similarity_graph.addEdge(nodeA, nodeB)
#        correlation_value[edge] = result
        
        
#      for s in correlation_properties:
#        if result > s:
#          edge = graph.addEdge(nodeA, nodeB)
#          correlation_value[edge] = result
#          
#          
#          correlation_properties[s][edge] = result
#          continue
#        correlation_properties[s][edge] = result

#  return {seuil:correlation_value}
  return correlation_value



def cosineSimilarity(v1,v2):
    
  sumxx, sumxy, sumyy = 0, 0, 0
  for i in range(len(v1)):
    x = v1[i]; y = v2[i]
    sumxx += x*x
    sumyy += y*y
    sumxy += x*y
  return sumxy/sqrt(sumxx*sumyy)

def jackknife(v1, v2):
  
  sumx_all, sumy_all, sumxx_all, sumyy_all, sumxy_all = pearson(v1, v2, jk_flag = True)
  size =  len(v1) -1
  result = []
  
  for i in range(len(v1)):
    x = v1[i]
    y = v2[i]
    
    sumxx = sumxx_all - x*x
    sumyy = sumyy_all - y*y
    sumxy = sumxy_all - x*y
    sumx = sumx_all - x
    sumy = sumy_all - y
    numerator = sumxy - ((sumx * sumy)/size)
    denominator = sqrt((sumxx - (sumx*sumx) / size) * (sumyy - (sumy*sumy)/size))
    
    if denominator == 0:
      continue
    result.append(numerator/denominator)
    
  return min(result)
    
    
def pearson(v1, v2, jk_flag = False):
  
  sumx = sum(v1)
  sumy = sum(v2)
  size = len(v1)
  
  sumxx, sumxy, sumyy = 0, 0, 0
  for i in range(len(v1)):
 
    x = v1[i]
    y = v2[i]
    
    sumxx += x*x
    sumyy += y*y
    sumxy += x*y

  numerator = sumxy - ((sumx * sumy)/size)
#  print((sumxx - (sumx*sumx) / size) * (sumyy - (sumy*sumy)/size))
  
  denominator = sqrt((sumxx - (sumx*sumx) / size) * (sumyy - (sumy*sumy)/size))
    
    
  if jk_flag:
    return sumx, sumy, sumxx, sumyy, sumxy
  if denominator == 0:
    return 0
 
  return numerator/denominator


def clusteringMCL(graph, weight_prop):
  params = tlp.getDefaultPluginParameters('MCL Clustering', graph)
  
  params['weights'] = weight_prop

  
  resultMetric = graph.getDoubleProperty('cluster')
  print("MCL clustering start")
  success = graph.applyDoubleAlgorithm('MCL Clustering', resultMetric, params)
  print("MCL clustering is over", success)
  return resultMetric
  
  
def writeClusters(graph, grp_metric, prop_to_write, nb_gene_threshold = 20):
  clusters = getClusters(graph, grp_metric)
  print(clusters.keys())
  for c in clusters:
    if len(clusters[c]) <  nb_gene_threshold:
      continue
    with open("Cluster_{}.txt".format(c), "w") as fl:
      for n in clusters[c]:
        fl.write(prop_to_write[n]+"\n")
        

def getClusters(graph, metric):
  clusters = {}
  for n in graph.getNodes():
    # le try except peut être vu comme maladroit mais python apprécie particulierement cette structure
    try:
      clusters[metric[n]].append(n)
    except:
      clusters[metric[n]] = [n]
  return clusters
  
def buildHeatMap(graph, heatmap, all_tp, cluster_metric):
  
  viewSize = graph.getSizeProperty("viewSize")
  viewColor = graph.getColorProperty("viewColor")
  viewLayout = graph.getLayoutProperty("viewLayout")
  viewBorderColor = graph.getColorProperty("viewBorderColor")
  
  clusters = getClusters(graph, cluster_metric)
  print(clusters)
  # Recherche de la valeur d'expression min et max pour l'intensitÃ© de couleur du heatmap
  # On met graph en parametre pour bien prendre en compte seulement les niveaux d'expression du graph pour lequel on fait le heatmap
  # Si jamais on a filtrÃ© ce graph lÃ , le max ou min peut diffÃ©rer du graph original
  minimum = all_tp[0].getNodeMin(graph) 
  maximum = all_tp[0].getNodeMax(graph) 
  
  for tp in all_tp:
    minimum = all_tp[0].getNodeMin(graph) if all_tp[0].getNodeMin(graph) < minimum else minimum
    maximum = all_tp[0].geNodetMax(graph) if all_tp[0].getNodeMax(graph) < maximum else maximum 
 
  print(minimum)
  print(maximum) 
  heatmap.clear() # Remove all nodes, edges and sub-graphs from the graph 
  heatmap_color = tlp.ColorScale([tlp.Color.Green, tlp.Color.Black, tlp.Color.Red])
  y_step = 2
  x_step = 100
  
  print("start")
  yi = 0
  
  for c in clusters:
    yi += 20
    for n in clusters[c]:
      yi += y_step
      
      xi =0
      for tp in all_tp:
        xi += x_step
        new_n = heatmap.addNode()
        viewSize[new_n] = tlp.Size(x_step,y_step,1)
        viewLayout[new_n] = tlp.Coord(xi, yi, 0)
  #     
        pos = (tp[n] - minimum )/ (maximum - minimum)
   
        viewColor[new_n] = heatmap_color.getColorAtPos(pos)
        viewBorderColor[new_n] = heatmap_color.getColorAtPos(pos)
  
  # Ajout de la legende tp dans le heatmap
  viewLabel = heatmap.getStringProperty("viewLabel")
  viewLabelPosition = heatmap.getIntegerProperty("viewLabelPosition")
  viewFontSize = heatmap.getIntegerProperty("viewFontSize")
  xi = 0
  yi = - y_step*40
  print(yi)
  
  # display number of TP at the bottom of the heatmap
  for i in range(len(all_tp)):
    xi += x_step
    new_n = heatmap.addNode()
    viewLabel[new_n] = str(i+1)
    viewLayout[new_n] = tlp.Coord(xi, yi, 1)
    viewFontSize[new_n] = 400
    viewSize[new_n] = tlp.Size(x_step,y_step*40,1)
    viewColor[new_n] = tlp.Color.White
    viewBorderColor[new_n] = tlp.Color.White
     
def reinitializeSubgraph(graph, name, clone_name = "clone"):
  
  if graph.getSubGraph(name) is None:
    print("creation of {} subGraph".format(name))
    sub = graph.addCloneSubGraph(name)
  else:
    sub = graph.getSubGraph(name)
 
  sub.clear()
  clone = graph.getSubGraph(clone_name)
  tlp.copyToGraph(sub, clone)
  return sub

def createCleanCopy(graph, graph_to_copy, copy_name):
  
  copy = graph.addCloneSubGraph(copy_name)
  copy.clear()

  tlp.copyToGraph(copy, graph_to_copy)
  return copy


  


def multiLevelClustering(multi_scale_graph, similarity_graph, nodes_cluster, corr_metric, nb_gene_min=20, level_max=3, level=1, group=0):
  
  print("level {} groupe {}".format(level, group))
  
  #Creation d'un sous graph qui va contenir tous les neouds du cluster
  similarity_graph.inducedSubGraph(nodes_cluster, multi_scale_graph, name="similarity_level{}_grp{}".format(level, group))
  graph_lvl = multi_scale_graph.getSubGraph("similarity_level{}_grp{}".format(level, group))
  
  # Calcule du seuil pour filtrer les aretes (ici la médiane des aretes de correlation)
  correlations  = [corr_metric[e] for e in graph_lvl.getEdges()]
  seuil = median(correlations)
  print("seuil mediane", seuil)
  
  percentile_75 = percentile(correlations, 75)
  print("percentile 75", percentile_75)
  
  
  # Filtrage des aretes du graph_lvl selon le seuil calculé
  for e in graph_lvl.getEdges():
    if corr_metric[e] < seuil:
      graph_lvl.delEdge(e)
  
  # Etape de Partionnement
  metric = clusteringMCL(graph_lvl, weight_prop=corr_metric)
  
  # Recupération des groupes définit par MCL
  clusters = getClusters(graph_lvl, metric)
  
  # Si la clusteriastion n'a renvoyé qu'un seul cluster on stope la clusterisation
  if len(clusters) == 1 or level >= level_max:
    clusters["len"] = graph_lvl.numberOfNodes()
    return clusters
  
  #Recursion de la fonction sur chaque nouveau groupe si le nombre de genes est suffisant
  for c in clusters:
    if len(clusters[c]) >= nb_gene_min:
      clusters[c] = multiLevelClustering(graph_lvl, similarity_graph, nodes_cluster=clusters[c], level=level+1, group=int(c), corr_metric=corr_metric)
    
   
 
  print(clusters)
  
  #A la fin on retourne clusters qui correspond à un dictionnaire multi echel
  clusters["len"] = graph_lvl.numberOfNodes()
  return clusters
  
  

def main(graph): 
  Locus = graph.getStringProperty("Locus")
  Negative = graph.getBooleanProperty("Negative")
  Positive = graph.getBooleanProperty("Positive")
  
  
  all_tp = [] 
  
  # Creation d'une liste comprenant tous les doubleProperty tp
  for i in range(1, 18):
    k = "tp{} s".format(i)
    tp = graph.getDoubleProperty("tp{} s".format(i))
    all_tp.append(tp)
  
  
  
  
  viewBorderColor = graph.getColorProperty("viewBorderColor")
  viewBorderWidth = graph.getDoubleProperty("viewBorderWidth")
  viewColor = graph.getColorProperty("viewColor")
  viewFont = graph.getStringProperty("viewFont")
  viewFontSize = graph.getIntegerProperty("viewFontSize")
  viewIcon = graph.getStringProperty("viewIcon")
  viewLabel = graph.getStringProperty("viewLabel")
  viewLabelBorderColor = graph.getColorProperty("viewLabelBorderColor")
  viewLabelBorderWidth = graph.getDoubleProperty("viewLabelBorderWidth")
  viewLabelColor = graph.getColorProperty("viewLabelColor")
  viewLabelPosition = graph.getIntegerProperty("viewLabelPosition")
  viewLayout = graph.getLayoutProperty("viewLayout")
  viewMetric = graph.getDoubleProperty("viewMetric")
  viewRotation = graph.getDoubleProperty("viewRotation")
  viewSelection = graph.getBooleanProperty("viewSelection")
  viewShape = graph.getIntegerProperty("viewShape")
  viewSize = graph.getSizeProperty("viewSize")
  viewSrcAnchorShape = graph.getIntegerProperty("viewSrcAnchorShape")
  viewSrcAnchorSize = graph.getSizeProperty("viewSrcAnchorSize")
  viewTexture = graph.getStringProperty("viewTexture")
  viewTgtAnchorShape = graph.getIntegerProperty("viewTgtAnchorShape")
  viewTgtAnchorSize = graph.getSizeProperty("viewTgtAnchorSize")

  # Creation of clone subgraph if no other graph have been created yet, otherwise "clone" already exist
  if graph.numberOfDescendantGraphs() == 1:
    print("creation of graph clone")
    graph.addCloneSubGraph("clone")
  

  #

#  correlation = graph.addCloneSubGraph("correlation")
  
#  working = reinitializeSubgraph(graph, "working")
#  heatmap = reinitializeSubgraph(working, "heatmap")
  working = graph.getSubGraph("working")
  
  


  #graph_working = graph.delSubGraph(graph_working)
  #working = graph.addCloneSubGraph("working")
  #correlation = graph.getSubGraph("correlation")
  
 
  
  label(working, Locus, viewLabel)
  setSize(working, viewSize)
  changeEdge(working, Negative, Positive, viewColor)
  NodePosition(working, viewLayout)
    
  #  similarity = reinitializeSubgraph(graph, "similarity")
  
  # Utilisation de la fonction jackknife pour le calcule de correlation
  similarityFct=jackknife
  
  if graph.getSubGraph("similarity") is None:
    
    print("creation of similarity subGraph")
    
    similarity = graph.addCloneSubGraph("similarity")
    
    print("Computation of similarity among genes")
    computeSimilarity(similarity, all_tp, similarityFct)
    
  else:
    similarity = graph.getSubGraph('similarity')
  
  if graph.getSubGraph("similarity_multi_scale") is None:
    multi_scale = graph.addCloneSubGraph('similarity_multi_scale')
    multi_scale.clear() # netoyage de multi_scale graph car on a pas besoin de ses nodes et edges 
  
  else:
    multi_scale  = graph.getSubGraph("similarity_multi_scale")
  
  all_nodes = list(similarity.getNodes()) #pour le premier niveau de partionnement tout les neouds sont requis
  correlation_prop = similarity.getDoubleProperty("correlation_{}".format(similarityFct.__name__))
  
  multiLevelClustering(multi_scale, similarity, nodes_cluster=all_nodes, corr_metric=correlation_prop)
  
  

  
  
  
  
    
#  cosineScipy = spatial.distance.cosine

#  correlation_graph = geneClustering(graph, all_tp, similarityFct=jackknife)
  
#  correlation_graph = geneClustering(graph, all_tp, similarityFct=cosineSimilarity,  seuils = [0, 0.5, 0.8, 0.95, 0.99, 0.98, 0.995])
  
#  correlation_graph = geneClustering(graph, all_tp, similarityFct=cosineSimilarity,  seuil=0.8)
  
  
  
#  correlation_graph = graph.getSubGraph("correlation_jackknife")

#  
#  seuils = [0.8]  #[0, 0.5, 0.8, 0.95, 0.99, 0.98, 0.995]
#  similarity_prop = {}
#  for s in seuils:
#    similarity_prop[s] = correlation_graph.getDoubleProperty("correlation_{}".format(s))
    
#  similarity_prop = {0.995:correlation_graph.getDoubleProperty("correlation_O.995")}

#  clusteringMCL(graph, similarity_prop)

#  cluster_metric = correlation_graph.getDoubleProperty("resultMetric_0.8")

#  buildHeatMap(correlation_graph, heatmap, all_tp, cluster_metric)
  
  
  
#  writeClusters(correlation_graph, cluster_metric, Locus)
  
  print('end')
  #working = graph.addCloneSubGraph("working")
  #clone = graph.addCloneSubGraph("clone")
  #correlation = graph.addCloneSubGraph("correlation")
