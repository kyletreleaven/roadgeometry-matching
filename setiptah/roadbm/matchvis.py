
import itertools

import numpy as np
import bintrees

import networkx as nx

""" my dependencies """
import setiptah.roadgeometry.roadmap_basic as ROAD
import setiptah.roadgeometry.astar_basic as ASTAR

import setiptah.roadbm.bm as roadbm

import matplotlib.pyplot as plt



""" convenience functions """

def position( address, roadmap, pos, length_attr='length' ) :
    """
    get the Euclidean position of a street address,
    given roadmap and dictionary of vertex positions
    """
    if isinstance( address, ROAD.RoadAddress ) :
        road = address.road
        coord = address.coord
    else :
        road, coord = address
    coord = float( coord )
        
    u,v, key = ROAD.obtain_edge( roadmap, road )
    assert key == road
    data = ROAD.get_road_data( road, roadmap )
    width = data.get( length_attr, 1 )
    
    #ROAD.get_edge_data( )
    x = pos[u]
    vec = pos[v] - x
    return x + vec * coord / width

def pointsToXY( points ) :
    """ split a list of (x,y) coordinates into X and Y; usually for plotting """
    X = [ x for x,y in points ]
    Y = [ y for x,y in points ]
    return X, Y


ZNODES = 1
ZLABELS = 2
ZEDGES = 3
ZTRAILS = 4
ZPOINTS = 5

def drawRoadmap( roadmap, pos, ax=None, **kwargs ) :
    if ax is None : ax = plt.gca()
    
    # draw the skeleton (undirected)
    skeleton = nx.convert.convert_to_undirected( roadmap )
    nx.draw_networkx_nodes( skeleton, pos, ax=ax, zorder=ZNODES )
    nx.draw_networkx_edges( skeleton, pos=pos, ax=ax, zorder=ZEDGES, **kwargs )
    
    road_labels = { (u,v) : road + '\n'     # the endline is to raise the label
                   for u,v,road in roadmap.edges_iter( keys=True ) }
    nx.draw_networkx_edge_labels( skeleton, pos=pos, ax=ax, 
                                  edge_labels=road_labels, zorder=ZLABELS )





def SHOWMATCH( match, S, T, roadmap, pos, length_attr='length', ax=None,
               **kwargs ) :
    """
    visualize a matching on a roadmap:
    imagine depositing one uniform trail of ink,
    for each match in the matching,
    on the shortest path between the endpoints of the match;
    segments of the network more often covered will obtain more ink
    """
    # labels for three kinds of graph nodes
    VERTEX = 'v'
    POINT_IN_S = 'S'
    POINT_IN_T = 'T'
    
    # draw the roadmap
    if ax is None : ax = plt.gca()
    options = { 'edge_color' : 'g', 'alpha' : .15 }     # lightly, though...
    options.update( kwargs )                            # but let overrides
    drawRoadmap( roadmap, pos, ax=ax, **options )
    ax.set_aspect('equal')  # i just like equal aspect...
    
    """ The hard part is getting the edges with proper thickness """
    # sort points onto segments
    segments = roadbm.SEGMENTS( S, T, roadmap )
    
    # make a path graph
    graph = nx.Graph()
    
    for u, v, road, data in roadmap.edges_iter( keys=True, data=True ) :
        width = data.get( length_attr, 1 )
        
        def traverse() :
            yield 0., VERTEX, u     # location, type, label
            for y, queue in segments[road].iter_items() :
                for s in queue.P : yield y, POINT_IN_S, s
                for t in queue.Q : yield y, POINT_IN_T, t
            yield width, VERTEX, v
            
        ITER = traverse()
        next = ITER.next()
        for y2, type2, label2 in ITER :
            y1, type1, label1 = next
            graph.add_edge( (type1,label1), (type2,label2), weight=y2-y1, score=0 )
            next = y2, type2, label2
            
    # add unit weight to shortest paths
    for i, j in match :
        path = nx.shortest_path( graph, (POINT_IN_S,i), (POINT_IN_T,j),
                                 weight='weight' )
        
        for ii, jj in zip( path[:-1], path[1:] ) :
            data = graph.get_edge_data( ii, jj )
            data['score'] += 1
            
    def vertpos(u) : return pos[u]
    def pos_from_S(u) : return position( S[u], roadmap, pos )
    def pos_from_T(u) : return position( T[u], roadmap, pos )
    switch = { VERTEX : vertpos, 
              POINT_IN_S : pos_from_S, 
              POINT_IN_T : pos_from_T }
    
    other_pos = {}
    for uu in graph.nodes_iter() :
        typeu, labelu = uu
        other_pos[uu] = switch[typeu]( labelu )
    
    if False :
        # figure out how to do this?
        colors = [ data['score'] for _,__,data in graph.edges_iter( data=True ) ]
        nx.draw_networkx_edges( graph, pos=other_pos, edge_color=colors )
        plt.colorbar()  #?
        #nx.draw(G,pos,node_color='#A0CBE2',edge_color=colors,width=4,edge_cmap=plt.cm.Blues,with_labels=False)
    else :
        # plot edges in graph with variable thickness? or some other visual cue
        for uu, vv, data in graph.edges_iter( data=True ) :
            score = data['score']
            if score <= 0 : continue    # would just waste effort
            
            posu = other_pos[uu]
            posv = other_pos[vv]
            
            xu, yu = posu
            xv, yv = posv
            options = { 'color' : 'k',
                       'alpha' : .6,
                        }
            ax.plot( [xu,xv], [yu,yv], linewidth=score, solid_capstyle='butt',
                     # butt style prevents awkward overlap of segments
                     zorder=ZTRAILS,
                     **options )
        
    # plot the points on top, so visible; this isn't working
    # show S points in red
    positions = [ position(addr, roadmap, pos) for addr in S ]
    options = { 'marker' : 'x' }
    X, Y = pointsToXY( positions )
    ax.scatter( X, Y, color='r', zorder=ZPOINTS, **options )
    # show T points in blue
    positions = [ position(addr, roadmap, pos).tolist() for addr in T ]
    X, Y = pointsToXY( positions )
    ax.scatter( X, Y, color='b', zorder=ZPOINTS, **options )





if __name__ == '__main__' :
    plt.close('all')
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--number', type=int, default=10 )
    args = parser.parse_args()
    
    N = args.number
    interchanges = [ np.random.rand(2) for i in xrange(N) ]
    
    import scipy.spatial as spatial
    """ Connect points using a Delaunay triangulation """
    tri = spatial.Delaunay( interchanges )
    
    import networkx as nx
    graph = nx.Graph()
    # find the edges in the triangulation
    indices, seq = tri.vertex_neighbor_vertices
    for i in xrange(N) :
        for j in seq[ indices[i]:indices[i+1] ] :
            graph.add_edge(i,j)
    
    # construct the roadmap
    roadmap = nx.MultiDiGraph()
    for ridx, (u,v) in enumerate( graph.edges() ) :
        x, y = [ tri.points[k] for k in (u,v) ]
        roadmap.add_edge(u,v, 'road %d' % ridx, length=np.linalg.norm(y-x) )
        
    # and positions
    pos = { k : point for k, point in enumerate( tri.points ) }
    
    
    
    # now, draw two sets of points
    import setiptah.roadgeometry.probability as roadprob
    uniform = roadprob.UniformDist( roadmap )
    
    
    
    if False :
        addresses = [ uniform.sample() for i in xrange(M) ]
        positions = [ position(addr, roadmap, pos) for addr in addresses ]
        
        X = [ x for x,y in positions ]
        Y = [ y for x,y in positions ]
        
        nx.draw(graph, pos=pos)
        ax = plt.gca()
        ax.scatter(X,Y)
        ax.set_aspect('equal')
    
    if True :
        # requires bintrees, too.
        #import setiptah.roadgeometry.roadmap_basic as ROAD
        
        M = 500
        unpack = lambda addr : ( addr.road, addr.coord )
        SS = [ unpack( uniform.sample() ) for i in xrange(M) ]
        TT = [ unpack( uniform.sample() ) for i in xrange(M) ]
        
        import setiptah.roadbm.bm as roadbm
        match = roadbm.ROADSBIPARTITEMATCH( SS, TT, roadmap )
        #match = [ (i,i) for i in xrange(M) ]
        SHOWMATCH( match[:], SS, TT, roadmap, pos=pos )
        
        
        
        
        
        
