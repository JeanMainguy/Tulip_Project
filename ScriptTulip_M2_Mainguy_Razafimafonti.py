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

def preTraitementVisualisation(graph):
  """
  Pré-traitement & première visualisation du graph
  """


  # Partie 1

  viewLabel = graph.getStringProperty("viewLabel")
  Locus = graph.getStringProperty("Locus")
  viewSize = graph.getSizeProperty("viewSize")
  viewFontSize = graph.getIntegerProperty("viewFontSize")
  viewColor = graph.getColorProperty("viewColor")
  viewTgtAnchorShape = graph.getIntegerProperty("viewTgtAnchorShape")
  viewTgtAnchorSize = graph.getSizeProperty("viewTgtAnchorSize")


  node_size = tlp.Size(100000,30000,1)


  for n in graph.getNodes():

    #Ajoute des labels issu de Locus aux neouds du graph
    viewLabel[n] = Locus[n]
    viewFontSize[n] = 1

    #Affection d'une taille non null à chaque sommet afin de pouvoir visualiser correctement les étiquettes
    viewSize[n] = node_size

  # Affectation des couleurs des arretes et la forme des extremités selon le type de régulation (positive ou negative)
  Negative = graph.getBooleanProperty("Negative")
  blue = tlp.Color(58,95,205)
  square = tlp.EdgeExtremityShape.Square

  Positive = graph.getBooleanProperty("Positive")
  green = tlp.Color(69,139,0)
  arrow = tlp.EdgeExtremityShape.Arrow

  for e in graph.getEdges():

    if Negative[e] is True:
      #Regulation negatif
      viewColor[e] = blue
      viewTgtAnchorShape[e] = square

    elif Positive[e] is True:
      #regulation poistive
      viewColor[e] = green
      viewTgtAnchorShape[e] = arrow

      viewTgtAnchorSize[e] = node_size


  #Affectation des positions aux sommets du graphe par l'application un model de force "FM^3 (OGDF)"
  viewLayout = graph.getLayoutProperty("viewLayout")

  fm3pParams = tlp.getDefaultPluginParameters("FM^3 (OGDF)", graph)
  fm3pParams["Unit edge length"] = 400
  graph.applyLayoutAlgorithm("FM^3 (OGDF)", viewLayout, fm3pParams)



def geneClustering(graph, all_tp, similarityFct, seuil):
  """
  Fonction pas multiechelle pour le clustering
  Donc fonction n'est plus utlisé
  """

  # creation ou reinitialisation du graph de correlation
  graph_correlation_name = "correlation_{}".format(similarityFct.__name__)
  correlation_graph = reinitializeSubgraph(graph, graph_correlation_name)


  # Supression de toutes les arretes
  print("Supression des aretes du graphe ",graph_correlation_name)
  correlation_graph.delEdges(correlation_graph.getEdges())

  # creation de la doubleProperty "correlation" permettant d'enregistrer
  # la correlation entre chaque gene au niveau des arretes


  similarity_prop = computeCorrelation(correlation_graph, all_tp, similarityFct, seuil)

  clusteringMCL(correlation_graph, similarity_prop )

  return correlation_graph



def computeSimilarity(similarity_graph, all_tp, similarityFct):
  """
  Va calculer la similiraité entre chaque gène.
  Filtration sur le gènes qui présentent un niveau d'expression de 0 à tous les temps

  parametres:
    similarity_graph: graph qui va contenir toutes les aretes entre chaque neouds portant une valeur de similarite pour leurs niveau d'exression
    all_tp :liste des propriete des niveaux d'expression a chaque temps
    similarityFct: function calculant une similarite entre deux liste de valeur. Peut correspondre a la fonction cosineSimilarity, jackknife ou pearson
  """


  correlation_value = similarity_graph.getDoubleProperty("correlation_{}".format(similarityFct.__name__))


  nannodes = []
  treated_node = []
  cpt = 0.0
  step = 0
  nodes = list(similarity_graph.getNodes()) # creation du liste contenant tous les noeuds qui va ^etre ensuite deplété pour pas passer deux fois sur la m^eme pair de neouds
  nb_node = len(nodes)
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
      if result < 0:
        continue

      #Add edge between node A and node B
      edge = similarity_graph.addEdge(nodeA, nodeB)
      correlation_value[edge] = result


    cpt += 1
    if cpt%50 == 0:
      print("{} gene already process..".format(int(cpt)))

  return correlation_value



def cosineSimilarity(v1,v2):
  """
  Calcule de similarité de Cosine etre v1 et v2
  """

  sumxx, sumxy, sumyy = 0, 0, 0

  for i in range(len(v1)):
    x = v1[i]; y = v2[i]
    sumxx += x*x
    sumyy += y*y
    sumxy += x*y
  return sumxy/sqrt(sumxx*sumyy)

def jackknife(v1, v2):
  """
  Calcule de similarité de Jackknife etre v1 et v2
  """

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

  """
  Calcule de la correlation de pearson
  Fonction également utilsé par JK
  """

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


  denominator = sqrt((sumxx - (sumx*sumx) / size) * (sumyy - (sumy*sumy)/size))


  if jk_flag:
    return sumx, sumy, sumxx, sumyy, sumxy
  if denominator == 0:
    return 0

  return numerator/denominator



def clusteringMCL(graph, weight_prop, result_name):
  """
  Effectue un partionnement race au plugin MCL clustering
  """
  params = tlp.getDefaultPluginParameters('MCL Clustering', graph)

  params['weights'] = weight_prop


  resultMetric = graph.getDoubleProperty(result_name)
  print("MCL clustering start")
  success = graph.applyDoubleAlgorithm('MCL Clustering', resultMetric, params)
  print("MCL clustering is over", success)
  return resultMetric


def extractClusters(cluster,all_tp,  niveau, nb_gene_threshold, prop_to_write):
  """
  Fonction récursive pour faire un heatmap pour chaque groupe > nb_gene_thrsholde
  Et pour chacun de ces groupes écrit dan un fichier la liste des
  """

  #le cluster un un nombre de gene > au threshold
  print("cluster type", type(cluster))
  for c in cluster:
    if type(cluster[c]) == dict:
      if cluster[c]['len'] >=  nb_gene_threshold:
        niveau_c = niveau+'_'+str(int(c))

        #Creation d'une heatmap pour ce cluster
        print(cluster[c]['len'])
        print("creation de la heatmap heatmap{} à partir de l ancestre".format(niveau))
        new_heatmap = 'heatmap_{}'.format(niveau_c)
        createHeatmap(graph, heatmap_name=new_heatmap, all_tp=all_tp, lvls_clusters=cluster[c])

        #ecriture du groupe dans un fichier pour l'enrichissmeent
        fl = open("Groupe/groupe_{}".format(niveau_c), 'w')
        writeClusters(cluster, fl, prop_to_write)
        fl.close()

        #Exploration des autres clusters en aval
        extractClusters(cluster[c], all_tp, niveau+'_'+str(int(c)), nb_gene_threshold, prop_to_write)

    if type(cluster[c]) == list:
      if len(cluster[c]) > nb_gene_threshold:
        #Creation d'une heatmap pour ce cluster mais ce sera le dernier pour cette branch
#        print("nb gene dans le groupe ", len(cluster[c]))
#        print("creation de la heatmap heatmap{}_{} à partir de l ancestre".format(niveau, int(c)))
#        print("pas de groupe en aval")
        niveau_c = niveau+'_'+str(int(c))
        new_heatmap = 'heatmap_{}_{}'.format(niveau, str(int(c)))

        #reconstruction d'un dico à un seul groupe pour le donner à createHeatmap
        dico_cluster = {c:cluster[c], 'len':len(cluster[c])}

        createHeatmap(graph, heatmap_name=new_heatmap, all_tp=all_tp, lvls_clusters=dico_cluster)
        #ecriture du cluster
        fl = open("Groupe/groupe_{}".format(niveau_c), 'w')
        writeClusters(cluster, fl, prop_to_write)
        fl.close()



def writeClusters(clusters, fl, prop):
  for c in clusters:
    if type(clusters[c]) == list:
      for n in clusters[c]:
        fl.write(prop[n]+"\n")

    elif type(clusters[c]) == dict:
      writeClusters(clusters[c], fl, prop)


def getClusters(graph, metric):
  """
  Permet de récupéré les clusters dans un dico à prtir de la metric donner par MCL et un graph
  """
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



def getMultiScaleClusters(graph, grp_prop_name, lvl=1, grp=0):
  """
  Permet de reconstruire la structure de clusters obtenu dans multiLevelClustering à partir des graph
  """
  print(grp_prop_name+str(lvl))
  grp_prop = graph.getDoubleProperty(grp_prop_name+str(lvl))
  clusters = getClusters(graph, grp_prop)
  lvl += 1
  for c in clusters:
    subgraph_name = "similarity_level{}_grp{}".format(lvl, int(c))
    if graph.getSubGraph(subgraph_name) is not None:
      next_lvl_graph = graph.getSubGraph(subgraph_name)
      clusters[c] = getMultiScaleClusters(next_lvl_graph, grp_prop_name, lvl=lvl, grp=int(c))

  clusters['len'] =  graph.numberOfNodes()
  return clusters

def getNumberOfLevel(lvls_clusters, lvl):
  for c in lvls_clusters:
    if type(c) == str:
      continue

    if type(lvls_clusters[c]) == dict:
      return getNumberOfLevel(lvls_clusters[c], lvl= lvl+1)
  return lvl


def createHeatmap(graph, heatmap_name, all_tp, lvls_clusters):
  """
  creatHeatmap
  Prépare les paremetre pour lancer la fonction récursive de création de heatmap
  """

  # heatmap_root graph va donner le grpah parent du heatmap

  viewLayout = graph.getLayoutProperty("viewLayout")

  heatmap_root = graph.getSubGraph("heatmap_root")
  if heatmap_root is None:
    print('Creation de heatmap_root')
    heatmap_root = graph.addCloneSubGraph("heatmap_root")


  heatmap = heatmap_root.getSubGraph(heatmap_name)
  if heatmap is None:
    print('Creation de {}'.format(heatmap_name))
    heatmap = heatmap_root.addCloneSubGraph(heatmap_name)
  heatmap.clear() ;#Pour construire le heatmap on part d'un graph vierge

  #Parametre de la carte de chaleur
  nb_gene = lvls_clusters['len'] # nombre d'element en Y
  nb_tp = len(all_tp) # nombre d'element en X
  stepx = 10 # valeur en dure arbitraire
  ratio = 0.5 # ratio largeur/hauteur du heatmap
  stepy = ((nb_tp*stepx)/(ratio*nb_gene)) #determination de stepy en fonction du reste des parametre afin de respecter le ratio
  nb_gene_min = 1 # nombre minimal de gene dans un cluster que la heatmap doit afficher
#  lvl_max=float('inf') #nombre de niveau de partionnement maximal qu'on souhaite afficher


  #Calcul de la distance en x des nouds de l'arbre
  nb_lvl = getNumberOfLevel(lvls_clusters, lvl=1)
#  if lvl_max < nb_lvl:
#    nb_lvl = lvl_max
  print('Il y a {} niveau de partionnement'.format(nb_lvl))
  # step branch correspond a la longuer des branche de l'arbre entre deux niveau partionnemnts
  step_branch = (nb_tp*stepx)/(nb_lvl+1) #largeur du heatmap / nombre de pas (donc ici largeur de l'arbre == largeur heatmap)

  yi = 0
  minimum = all_tp[0].getNodeMin(graph)
  maximum = all_tp[0].getNodeMax(graph)

  for tp in all_tp:
    minimum = tp.getNodeMin(graph) if tp.getNodeMin(graph) < minimum else minimum
    maximum = tp.getNodeMax(graph) if tp.getNodeMax(graph) > maximum else maximum

  ynode = lvls_clusters["len"]*stepy/2
  xnode = -step_branch*(nb_lvl+1)
  ancestor_node = heatmap.addNode() #noeud correspondant à l'origine de l'arbre
  viewLayout[ancestor_node] = tlp.Coord(xnode, ynode, 0) # distant du heatmap d'une largeur de heatmap en X et se situant au milieu du heatmap en Y

  buildMultiScaleHeatmap(heatmap, lvls_clusters, all_tp, lvl=1, ancestor_node=ancestor_node, step_branch=step_branch, yi=yi, maximum=maximum, minimum=minimum, stepx=stepx, stepy=stepy, nb_gene_min=nb_gene_min)


def buildMultiScaleHeatmap(heatmap, clusters, all_tp, lvl, ancestor_node, step_branch, yi, maximum, minimum, stepx, stepy, nb_gene_min=1):
  """
  Recusrsivité
  lancer pour chaque groupe la creation du heatmap ou bien  si il est composé de sous groupe  se relance sur chaque sous groupe
  Gere egalement la création des noeuds et arretes de l'arbre
  """
  viewLayout = graph.getLayoutProperty("viewLayout")

  distance = len(all_tp)*stepx - step_branch*lvl #distance en X des neouds de niveau de partionnement courant

#  if lvl >= lvl_max: #si on a atteint le niveau de partionnement maximal on transforme cluster en list pour pas avoir de niveau superieur
#    flat_dico = {}
#    flat_dico[0] = getFlatList(clusters)
#    flat_dico['len'] =  clusters['len']
#    clusters = flat_dico
#    print(clusters[c])

  for c in clusters:

    #Passer la clé 'len'
    if c == 'len':
      continue


    #Si le cluster est une liste de neouds et non un nouveau dico
    if type(clusters[c]) == list:
      #Si le nombre de gene est inferieur au threshold on passe à la clusteriation suivante
      if nb_gene_min >  len(clusters[c]):
        continue


      yinext = yi + len(clusters[c]) * stepy
#      print(lvl*'    '+"affichage dans le heatmap du cluster groupe"+ str(c)+ " level "+ str(lvl)+"  yi "+str(yi)+"  nb gene" +str(len(clusters[c])) )
      # dans ce cas il n'y pas de niveau supplémentaire on affiche donc les neouds de c dans le heatmap
      displayClusterToHeatmap(heatmap, clusters[c], all_tp, yi,  maximum, minimum, stepx, stepy)

      #creer neouds au niveau de yi et yi +nbneouds * stepy   correspondant au feuilles extrem de chaque groupes
      first_leaf = heatmap.addNode()
      last_leaf = heatmap.addNode()

      viewLayout[first_leaf] = tlp.Coord(0-stepx/2, yi, 0)
      viewLayout[last_leaf] = tlp.Coord(0-stepx/2, yinext-stepy, 0)

      #creer un neouds à distance d de xi et à distance y0 + yi/2 # et donner le label correspondant au numéro du cluster

      current_node = createTreeNode(heatmap, distance, yi, yinext, str(int(c)))

      #relier les neouds precedent à ancestor link
      heatmap.addEdge(ancestor_node, current_node)
      heatmap.addEdge(current_node, first_leaf)
      heatmap.addEdge(current_node, last_leaf)


#
      yi = yinext
    else: #il y a un nouveau niveau à afficher
      yinext = yi + clusters[c]['len'] * stepy
      current_node = createTreeNode(heatmap, distance, yi, yinext, str(int(c)))
      heatmap.addEdge(ancestor_node, current_node)

      print(lvl*'    '+"nouveau niveau de cluster groupe"+ str(c)+ " level "+ str(lvl))
      #creation de ancestor link
      buildMultiScaleHeatmap(heatmap, clusters[c], all_tp, lvl+1, current_node, step_branch, yi, maximum, minimum, stepx, stepy, nb_gene_min)
      yi = yinext




def createTreeNode(heatmap, dx, yi, yinext, label):
  """
  Gere la creation d'un neouds du cluster
  """
  viewLayout = graph.getLayoutProperty("viewLayout")
  viewLabel = graph.getStringProperty("viewLabel")
  viewColor = heatmap.getColorProperty("viewColor")
  viewSize = heatmap.getSizeProperty("viewSize")

  ynode = yi + (yinext-yi)/2.0
  xnode = -dx
  node = heatmap.addNode()
  viewLayout[node] = tlp.Coord(xnode, ynode, 0)
  viewLabel[node] = label
  viewColor[node] = tlp.Color.White
  viewSize[node] = tlp.Size(1,1,0)
  return node



def displayClusterToHeatmap(heatmap, nodes, all_tp, yi, maximum, minimum, stepx, stepy, display_gene_name=False):
  """
  Fonction d'écriture du heatmap en lui même
  Pour chaque gene du groupe crée une ligne avec ces valeur d'expression
  """

  viewSize = heatmap.getSizeProperty("viewSize")
  viewColor = heatmap.getColorProperty("viewColor")
  viewLayout = heatmap.getLayoutProperty("viewLayout")
  viewBorderColor = heatmap.getColorProperty("viewBorderColor")
  viewLabel = graph.getStringProperty("viewLabel")
  Locus = graph.getStringProperty("Locus")
  viewFontSize = graph.getIntegerProperty("viewFontSize")

  heatmap_color = tlp.ColorScale([tlp.Color.Green, tlp.Color.Black, tlp.Color.Red])

  for n in nodes:

    xi =0
    for tpi, tp in enumerate(all_tp):

      new_n = heatmap.addNode()
      viewSize[new_n] = tlp.Size(stepx,stepy,1)
      viewLayout[new_n] = tlp.Coord(xi, yi, 0)

      pos = (tp[n] - minimum )/ (maximum - minimum)

      viewColor[new_n] = heatmap_color.getColorAtPos(pos)
      viewBorderColor[new_n] = heatmap_color.getColorAtPos(pos)
      xi += stepx


    ## ne marche pas encore
    if display_gene_name:
      new_n = heatmap.addNode()
      viewSize[new_n] = tlp.Size(1,1,1)
      viewLayout[new_n] = tlp.Coord(xi, yi, 0)

      viewLabel[new_n] = Locus[n]
      viewFontSize[new_n] = stepy
      viewColor[new_n] = tlp.Color.White
      viewBorderColor[new_n] = tlp.Color.White

    yi += stepy





def multiLevelClustering(multi_scale_graph, similarity_graph, nodes_cluster, corr_metric, nb_gene_min=25, level_max=float('inf'), level=1, group=0):
  
  """
  Fonction récursive pour faire le partionnement sur plusieur niveaux.
  Parametre :
    multi_scale_graph : graph parent des graph créer à un niveau donné
    similarity_graph : graph contenant les proproété avec les valeurs de similiraté
    nodes_cluster : les neouds qu'il faut partionner à ce tour
    corr_metric : propriété de similarty graph donnnant la mesure de similarité
    nb_gene_min : nb minimum de gene sur lequel on réalise un clustering
    level_max=float('inf') : nombre de niveau max
    level : niveau courant de partionnement 1 correspond au partionnement sur tous les neouds neouds
    group : numéro de group dans le niveau va e 0 au nombre de groupe que MCL a trouvé

  """

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
      
#  
#  # copy du graph level pour garder les info à chaque niveau et pour chaque groupe - permet de reconstruire le dico imbriqué à posteriori 
#  graph_lvl_copy = multi_scale_graph.getSubGraph("similarity_level{}_grp{}_copy".format(level, group))
  
  # Etape de Partionnement
  metric = clusteringMCL(graph_lvl, weight_prop=corr_metric, result_name = 'cluster_level{}'.format(level))
  
  # Recupération des groupes définit par MCL
  clusters = getClusters(graph_lvl, metric)
  
  # Si la clusteriastion n'a renvoyé qu'un seul cluster on stope la clusterisation
  if len(clusters) == 1 or level >= level_max:
    clusters["len"] = graph_lvl.numberOfNodes()
    return clusters
  
  #Recursion de la fonction sur chaque nouveau groupe si le nombre de genes est suffisant
  for c in clusters:
    if len(clusters[c]) >= nb_gene_min:
      clusters[c] = multiLevelClustering(graph_lvl, similarity_graph, nodes_cluster=clusters[c], level=level+1, group=int(c), corr_metric=corr_metric, level_max=level_max)
    
   
 
#  print(clusters)
  
  #A la fin on retourne clusters qui correspond à un dictionnaire multi echel
  clusters["len"] = graph_lvl.numberOfNodes()
  return clusters



def getGeneName(graph, filename="correpondance_table.txt"):

  Locus = graph.getStringProperty("Locus")

  try:
    fl = open(filename, 'r')
  except:
    #le fichier n'est pas là
    print("{} n'est pas dans le repertoir courant... On utlisera les ECK pour l'ecriture")
    return Locus

  GeneId = graph.getStringProperty("GeneName")
  dico_correspondance = {}
  for l in fl:
    champs = l.split('\t')
    dico_correspondance[champs[0]] = champs[1]

  cpt = 0
  for n in graph.getNodes():

    if Locus[n] in dico_correspondance:

      GeneId[n] = dico_correspondance[Locus[n]]
    else:
      cpt +=1
#      print(Locus[n]+" not in correspondance table")
      GeneId[n] = Locus[n] #danc ce cas là le Locus est ajouté

  return GeneId


def main(graph):
  Locus = graph.getStringProperty("Locus")


  all_tp = []

  # Creation d'une liste comprenant tous les doubleProperty tp
  for i in range(1, 18):
    k = "tp{} s".format(i)
    tp = graph.getDoubleProperty("tp{} s".format(i))
    all_tp.append(tp)


  # Creation of clone subgraph if no other graph have been created yet, otherwise "clone" already exist
  if graph.numberOfDescendantGraphs() == 1:
    print("creation of graph clone")
    graph.addCloneSubGraph("clone")

  #Ecrasment et/ou réation du graph working  à partir du graph clone
  working = reinitializeSubgraph(graph, "working")

  #Pretraitrement et Premiere Visualisation
  preTraitementVisualisation(working)

  
  ## MESURE DE SIMILARITE ENTRE LES GENES##

  # Utilisation de la fonction jackknife, cosineSimilarity ou pearsonle pour le calcule de correlation:

  similarityFct=cosineSimilarity
  similarityFct=jackknife
  similarityFct=pearson
  
  #Création ou recupération du graph similarity qui va contenit des aretes entre chaque genes et un poid associé correspondant à la similarité entre les deux gènes liés
  similarity = graph.getSubGraph('similarity')
  if similarity is None:

    print("creation of similarity subGraph")
    similarity = reinitializeSubgraph(graph, "similarity")

    print("Computation of similarity among genes")

    computeSimilarity(similarity, all_tp, similarityFct)   #Mesure de la similiratité entre tous les noeuds

  else:
    # Si le graph a deja été déterminé on le récupère sans avoir à recommencer les mesures
    similarity = graph.getSubGraph('similarity')

  
  #PARTIONNEMENT A PLUSIEUR ECHELLE##

#  similarity_multi_sert a être le graph parent des différent sous graph déterminé par le partionnement
  if graph.getSubGraph("similarity_multi_scale") is None:
    multi_scale = graph.addCloneSubGraph('similarity_multi_scale')
    multi_scale.clear() # netoyage de multi_scale graph car on a pas besoin de ses nodes et edges

  else:
    multi_scale  = graph.getSubGraph("similarity_multi_scale")

  all_nodes = list(similarity.getNodes()) #pour le premier niveau de partionnement tout les neouds sont requis

  # les mesures de similarité pour la fonction étudié. SImiliraty peut avoir les valeurs pour COsine et Jackknife dans deux proproeté différentes
  similarity_property = similarity.getDoubleProperty("correlation_{}".format(similarityFct.__name__))


  if multi_scale.getSubGraph("similarity_level1_grp0") is None:
    print("multiLevelClustering in progress..")
    #Fait le partionnement multi echelle en retournant un dictionnaire imbriqué et en créant des sous graph 
    lvls_clusters = multiLevelClustering(multi_scale, similarity, nodes_cluster=all_nodes, corr_metric=similarity_property)

  else:
    #Pour ne pas se coltiner l'étape de clustering et juste recupéré le dico imbriqué des différent niveau de partionnement
    print('getMultiScaleClsuters in progress... retrieval of the  multi level clusters from the graphs')
    graph_lvl1_grp0 = multi_scale.getSubGraph("similarity_level1_grp0")
    cluster_propriety = "cluster_level"
    lvls_clusters = getMultiScaleClusters(graph_lvl1_grp0, cluster_propriety, lvl=1, grp=0)



  createHeatmap(graph, heatmap_name="heatmap", all_tp=all_tp, lvls_clusters=lvls_clusters)


  #Partie pas completement stable mais marche quand même bien
  #Recuperation des noms des gènes pour l'enrichissment
  GeneId = getGeneName(multi_scale, filename="correpondance_table.txt")



#  Creation de heatmap en cascade pour les groupe ayant plus de gene que le seuil et écriture des gènes
  niveau = ''
  nb_gene_threshold = 20
  prop_to_write = Locus
  heatmap_ancestor = 'heatmap_0'

#
#  extractClusters(lvls_clusters, all_tp, niveau, nb_gene_threshold, prop_to_write=GeneId)
#




  print('end of the anlysis')
