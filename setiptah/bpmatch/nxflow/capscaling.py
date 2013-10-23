
#from bpmatch_roads.bm_roadnet import costWrapper

from setiptah.basic_graph.mygraph import mygraph
from setiptah.nxopt.cvxcostflow import MinConvexCostFlow

import setiptah.bpmatch.roadmaps as roadbm


class negativeWrapper :
    def __init__(self, func ) :
        self.func = func
        
    def __call__(self, z ) :
        return self.func( -z )




def SOLVER( roadnet, surplus, objectives ) :
    network = mygraph()
    capacity = {}
    supply = { i : 0. for i in roadnet.nodes() }
    cost = {}   # functions
    
    for i,j, road in roadnet.edges_iter( keys=True ) :
        supply[j] += surplus[road]
        
        cc = roadbm.costWrapper( objectives[road] )
        ncc = negativeWrapper( cc )     # don't really have to worry about the C(0) offset
        
        network.add_edge( (road,+1), i, j )
        cost[ (road,+1) ] = cc
        
        network.add_edge( (road,-1), j, i )
        cost[ (road,-1) ] = ncc
    
    f = MinConvexCostFlow( network, {}, supply, cost )
    
    flow = {}
    for i, j, road in roadnet.edges_iter( keys=True ) :
        flow[road] = f[(road,+1)] - f[(road,-1)]
        flow[road] = int( flow[road] )
    
    #print flow
    return flow





