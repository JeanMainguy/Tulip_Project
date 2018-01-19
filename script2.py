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
  fm3pParams["Unit edge length"] = 200
  graph.applyLayoutAlgorithm("FM^3 (OGDF)", viewLayout, fm3pParams)
 
  size = 200
  #for n in graph.getNodes():
   # x = random.random() * size
   # y = random.random() * size
   # viewLayout[n] = tlp.Coord(x, y, 0)


def geneClustering(graph, all_tp, seuil=0.95):
  
  # creation ou reinitialisation du graph de correlation portant le seuil dans son nom
  graph_correlation_name = "correlation_{}".format(seuil)
  correlation_graph = reinitializeSubgraph(graph, graph_correlation_name)
  
  # Supression de toutes les arretes
  correlation_graph.delEdges(correlation_graph.getEdges())
  
  # creation de la doubleProperty "correlation" permettant d'enregistrer 
  # la correlation entre chaque gene au niveau des arretes
  correlation_value = correlation_graph.getDoubleProperty("correlation")
  
  computeCorrelation(correlation_graph, all_tp, correlation_value, seuil)
  
  clusteringMCL(correlation_graph, correlation_value)
  
  return correlation_graph

def computeCorrelation(graph, all_tp, correlation_value, seuil):
  
  nannodes = []
  treated_node = []
  cpt = 0
  nodes = list(graph.getNodes())
  for nodeA in graph.getNodes():
    try:
      nodes.remove(nodeA)
    except ValueError:
      continue
   
    exprA = []
    for tp in all_tp:
      exprA.append(tp[nodeA])
    
    if sum(exprA) == 0:
      nannodes.append(nodeA)
      graph.delNode(nodeA)
      continue 
    
    cpt += 1
    print(cpt)

    for nodeB in nodes:
      exprB = []

      for tp in all_tp:
        exprB.append(tp[nodeB])
    
      if sum(exprB) == 0:
        nannodes.append(nodeB)
        graph.delNode(nodeB)
        nodes.remove(nodeB)
        continue
       
      result = np.corrcoef(exprA, exprB)[0][1]
      
      #Add edge between node A and node B
      if result < 0:
        continue
  
      if result < seuil:
        continue 
      
      edge = graph.addEdge(nodeA, nodeB)
      correlation_value[edge] = result


def clusteringMCL(graph, correlation_value):
  params = tlp.getDefaultPluginParameters('MCL Clustering', graph)
  
#  params['inflate'] = .
  params['weights'] = correlation_value
#  params['pruning'] = 5
  
  resultMetric = graph.getDoubleProperty('resultMetric')
  success = graph.applyDoubleAlgorithm('MCL Clustering', resultMetric, params)
  
  
  

def buildHeatmap(graph, heatmap, all_tp, extrem = {"min":0, "max":15}):
  
  viewSize = graph.getSizeProperty("viewSize")
  viewColor = graph.getColorProperty("viewColor")
  viewLayout = graph.getLayoutProperty("viewLayout")
  viewBorderColor = graph.getColorProperty("viewBorderColor")
  

  # Recherche de la valeur d'expression min et max pour l'intensité de couleur du heatmap
  # On met graph en parametre pour bien prendre en compte seulement les niveaux d'expression du graph pour lequel on fait le heatmap
  # Si jamais on a filtré ce graph là, le max ou min peut différer du graph original
  minimum = all_tp[0].getNodeMin(graph) 
  maximum = all_tp[0].getNodeMax(graph) 
  
  for tp in all_tp:
    minimum = all_tp[0].getNodeMin(graph) if all_tp[0].getNodeMin(graph) < minimum else minimum
    maximum = all_tp[0].geNodetMax(graph) if all_tp[0].getNodeMax(graph) < maximum else maximum 
 
  print(minimum)
  print(maximum) 
  heatmap.clear() # Remove all nodes, edges and sub-graphs from the graph 
  heatMap = tlp.ColorScale([tlp.Color.Green, tlp.Color.Black, tlp.Color.Red])
  y_step = 2
  x_step = 100
  
  print("start")
  yi = 0
  for n in graph.getNodes():
    yi += y_step
    
    xi =0
    for tp in all_tp:
      xi += x_step
      new_n = heatmap.addNode()
      viewSize[new_n] = tlp.Size(x_step,y_step,1)
      viewLayout[new_n] = tlp.Coord(xi, yi, 0)
#     
      pos = (tp[n] - extrem['min'] )/ (extrem['max'] - extrem['min'])
 
      viewColor[new_n] = heatMap.getColorAtPos(pos)
      viewBorderColor[new_n] = heatMap.getColorAtPos(pos)
  
  # Ajout de la legende tp dans le heatmap
  viewLabel = heatmap.getStringProperty("viewLabel")
  viewLabelPosition = heatmap.getIntegerProperty("viewLabelPosition")
  viewFontSize = heatmap.getIntegerProperty("viewFontSize")
  xi = 0
  yi = - y_step*40
  print(yi)
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

def createSafecopy(graph, copyname, nametocopy):
  

  copied_graph = graph.getSubGraph(nametocopy)
  
  safecopy = graph.addCloneSubGraph(copyname)
  safecopy.clear()

  tlp.copyToGraph(safecopy, copied_graph)

  

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
  
  working = reinitializeSubgraph(graph, "working")
  heatmap = reinitializeSubgraph(graph, "heatmap")
  #graph_working = graph.delSubGraph(graph_working)
  #working = graph.addCloneSubGraph("working")
  #correlation = graph.getSubGraph("correlation")
  
  
  
  label(graph, Locus, viewLabel)
  setSize(graph, viewSize)
  changeEdge(graph, Negative, Positive, viewColor)
  NodePosition(graph, viewLayout)


#  correlation_graph = geneClustering(graph, all_tp, seuil=0.95)
  correlation_graph = graph.getSubGraph("correlation_0.95")
  buildHeatmap(correlation_graph, heatmap, all_tp)
  
  
  print('end')
  #working = graph.addCloneSubGraph("working")
  #clone = graph.addCloneSubGraph("clone")
  #correlation = graph.addCloneSubGraph("correlation")
