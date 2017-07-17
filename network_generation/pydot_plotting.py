# pydot plotting 

import networkx as nx
import pydot
import Image
 
G = nx.DiGraph()
 
... some code that adds edges to G
 
Gp = nx.to_pydot(G)   #converts G to Gp in pydot format
 
... some code that makes the Gp pydot graph prettier
# the default for a DiGraph will have arrows
# i will often do the following
# edges = Gp.get_edge_list()
# for e in edges:
#      source = e.get_source()
#      target = e.get_destination()
#      if source == 'some string that is the name of the source with edges that you want to make look different' and target == 'some string2':
#        e.set_style('bold')
#        e.set_color('red')
#        e.set_arrowhead('tee')  # use tee instead of arrow
#        e.set_label('label for edge')
# can do similar stuff with nodes
 
# now output your graph to a file and display it
outstem = 'bounded'
Gp.write_png(outstem + '_dot.png', prog='dot')  # writes Gp to png file #use prog='neato' or 'fdp' for undirected graph (no default arrows for this)
# the next 2 lines open and display the png file
im = Image.open(outstem + '_dot.png') 
im.show()
 