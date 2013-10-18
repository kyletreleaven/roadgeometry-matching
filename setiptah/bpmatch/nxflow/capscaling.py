
#from bpmatch_roads.bm_roadnet import costWrapper

from setiptah.mygraph import mygraph
from setiptah.nxopt.cvxcostflow import MinConvexCostFlow


class costWrapper :
    def __init__(self, lines ) :
        self.lines = lines
        
    def __call__(self, z ) :
        """ this is an O(log n) query function, can be reduced to O(1) by random access with saturation """
        _, line = self.lines.floor_item( -z )
        return line( z )

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
        
        cc = costWrapper( objectives[road] )
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



