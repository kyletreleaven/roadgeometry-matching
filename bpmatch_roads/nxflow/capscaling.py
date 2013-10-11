
#from bpmatch_roads.bm_roadnet import costWrapper

from bpmatch_roads.mygraph import mygraph
from bpmatch_roads.cvxcostflow import MinConvexCostFlow


class costWrapper :
    def __init__(self, lines ) :
        self.lines = lines
        
    def __call__(self, z ) :
        """ this is an O(log n) query function, can be reduced to O(1) by random access with saturation """
        _, line = self.lines.floor_item( -z )
        return line( z )





def SOLVER( roadnet, surplus, objectives ) :
    network = mygraph()
    capacity = {}
    supply = { i : 0. for i in roadnet.nodes() }
    cost = {}   # functions
    
    for i,j, road in roadnet.edges_iter( keys=True ) :
        network.add_edge( road, i, j )
        
        # construct cost function for the road
        cost[road] = costWrapper( objectives[road] )
        
        supply[j] += surplus[road]
        
    flow = MinConvexCostFlow( network, {}, supply, cost )
    for e in flow : flow[e] = int( flow[e] )
    
    #print flow
    return flow



