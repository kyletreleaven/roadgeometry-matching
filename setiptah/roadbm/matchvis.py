
import itertools

import numpy as np
import bintrees

import networkx as nx

""" my dependencies """
import setiptah.roadgeometry.roadmap_basic as ROAD
import setiptah.roadgeometry.astar_basic as ASTAR

import setiptah.roadbm.bm as roadbm

import matplotlib.pyplot as plt













def heightFunctionTex( S, T, ymin, ymax, z=None, steps=None ) :
    if z is None :
        zplus = 0
    else :
        zplus = z
        
    #str = "\\begin{tikzpicture}\n"
    str = ''
    
    # compute necessary data
    segment = roadbm.ONESEGMENT( S, T )
    intervals = roadbm.INTERVALS( segment )
    
    # draw the axis
    # horizontal
    #str += "\\draw [->] (%f,0) -- (%f,0) node [right] {$\\coordvar$} ;\n" % ( YMIN, YMAX+.5 )
    str += texhline( 0., ymin, ymax + .5, style='->' )
    #
    str += "\\draw [thick] (%(ymax)f,-.15) -- (%(ymax)f,.15) " % { 'ymax' : ymax }
    str += "node [above right] {$\\roadlen_\\roadvar$} ;\n"
    # vertical
    fmin, fmax = min( intervals ), max( intervals )
    hmin, hmax = min( fmin + zplus, 0 ), max( fmax + zplus, 0 )
    hmin, hmax = np.floor( hmin ), np.ceil( hmax )
    #np.floor( min( intervals ) )
    #fmax = np.ceil( max( intervals ) )
    data = { 'ymin' : ymin, 'fmin' : hmin-.25, 'fmax' : hmax+.25 }
    str += "\\draw [->] (%(ymin)f,%(fmin)f) -- (%(ymin)f,%(fmax)f) " % data
    str += "node [above] {$\\postcumarcs(y) = \\cumarcs(y) + \\numarcs$} ;\n"
    # vertical ticks
    for tick in np.arange(hmin,hmax+1,1) :
        data = { 'tR' : ymin+.1, 'tL' : ymin-.1, 'tick' : tick }
        str += "\\draw (%(tR)f,%(tick)d) -- (%(tL)f,%(tick)d) " % data
        str += "node [left=5] {$%(tick)d$} ;\n" % data
        
    # place the X's and O's
    str += texhline( zplus, ymin, ymax, style='dashed' )
    for y, item in segment.items() :
        phases = [ (item.P, '${\\color{red}\\times}$' ), (item.Q, '${\\color{blue}\\circ}$') ]
        for Y, mark in phases :
            for i in Y : str += "\\draw (%f,%f) node {%s} ;\n" % ( y, zplus, mark )
            
    # place the levels
    arr = bintrees.RBTree()
    for f, ranges in intervals.items() :
        for a,b in ranges :
            lb = a
            if lb == '-' : lb = -np.inf
            arr[lb] = ( (a,b), f )
            
    #print arr
    #print [ v for v in arr.values() ]
    
    iter = arr.values()
    if steps is not None : iter = itertools.islice( iter, steps )
    for (a,b), f in iter :
    #for f, ranges in intervals.items() :
        h = f + zplus
        label = ''
        if a == '-' :
            a = ymin
            label = "node [%s] {$\\numarcs_\\roadvar$}"
            if h >= 0 :       # zr >= 0?
                label = label % "above"
            else :
                label = label % "below"
        if b == '+' :
            b = ymax
            str += "\\draw (%f,%f) node [right] {$\\numarcs_\\roadvar + \\surplus_\\roadvar$} ;\n" % (b,h)
                        
        levelstr = "\\draw [thick] (%(yl)f,%(z)f) -- %(extra)s (%(yr)f,%(z)f) ;\n"
        data = { 'yl' : a, 'yr' : b, 'z' : h, 'extra' : label }
        str += levelstr % data
        # place the shades
        str += "\\path [fill=black,opacity=.2] (%(yl)f,0) rectangle (%(yr)f,%(z)f) ;\n" % data
        
    #str += "\\end{tikzpicture}\n"
    return str




def discreteCostTex( S, T, ymin, ymax, zmin, zmax, zplus=None ) :
    segment = roadbm.ONESEGMENT( S, T )
    meas = roadbm.MEASURE( segment, ymin, ymax )
    obj = roadbm.OBJECTIVE( meas )
    Cf = roadbm.costWrapper( obj )
    
    Z = range(zmin,zmax+1)
    C = [ Cf(z) for z in Z ]
    cmin, cmax = min(C), max(C)
    
    # obtain the 'star'
    ZMIN, ZMAX = -max( meas ), -min( meas )
    CC = [ ( Cf(z), z ) for z in range(ZMIN,ZMAX+1) ]
    COPT, ZOPT = min(CC)
    
    
    #str = "\\begin{tikzpicture}[x=1cm,y=.1cm]\n"
    str = ''
    # draw the axis
    # horizontal
    horz_level = np.floor(cmin)
    data = { 'zmin' : zmin, 'zmax' : zmax, 'horz' : horz_level }
    fmt = "\\draw [->] (%(zmin)d,%(horz)d) -- (%(zmax)f,%(horz)d) node [right] {$\\coordvar$} ;\n" 
    str +=  fmt % data
    # vertical
    str += "\\draw [->] (0,%(cmin)f-1) -- (0,%(cmax)f) " % { 'cmin' : cmin, 'cmax' : cmax }
    str += "node [above] {$\\cost(\\numarcs)$} ;\n"
    # horizontal ticks
    for tick in np.arange(zmin,zmax+1,1) :
        data = { 'tick' : tick, 'horz' : horz_level }
        str += "\\draw (%(tick)d,%(horz)f-.1) -- (%(tick)f,%(horz)f+.1) " % data
        str += "node [below=5] {$%d$} ;\n" % tick
        
    # scatter
    for z, c in zip( Z, C ) :
        if z == ZOPT :
            str += '\\draw [fill] (%f,%f) node {$\\star$} ;\n' % (z,c)       # {$\\circ$}'
        else :
            str += '\\draw [fill] (%f,%f) circle (.05cm) ;\n' % (z,c)       # {$\\circ$}'
        
    if zplus is not None and zmin <= zplus and zplus <= zmax :
        str += '\\draw [] (%f,%f) circle (.2cm) ;\n' % (zplus,Cf(zplus))       # {$\\circ$}'
        
    #str += "\\end{tikzpicture}\n"
    return str
    
    
    




def costFunctionTex( S, T, ymin, ymax, z=None ) :
    if z is None :
        zplus = 0
    else :
        zplus = z
        
    """ C to tikz """
    segment = roadbm.ONESEGMENT( S, T )
    meas = roadbm.MEASURE( segment, ymin, ymax )
    obj = roadbm.OBJECTIVE( meas )
    Cf = roadbm.costWrapper( obj )
    
    str = "\\begin{tikzpicture}[x=2cm,y=.1cm]\n"
    
    # draw the axis
    ZMIN, ZMAX = -max( meas ), -min( meas )
    WIDTH = ZMAX - ZMIN
    
    C = [ ( Cf(z), z ) for z in range(ZMIN,ZMAX+1) ]
    CMIN, ZOPT = min(C)
    CMAX, _ = max(C)
    COPT = Cf(ZOPT)
    
    ZMIN = int( min( ZMIN, np.floor( zplus ) ) )
    ZMAX = int( max( ZMAX, np.ceil( zplus ) ) )
    
    # horizontal
    horz_level = np.floor(CMIN)
    data = { 'zmin' : ZMIN, 'zmax' : ZMAX, 'horz' : horz_level }
    fmt = "\\draw [->] (%(zmin)d,%(horz)d) -- (%(zmax)f,%(horz)d) node [right] {$\\coordvar$} ;\n" 
    str +=  fmt % data
    # vertical
    str += "\\draw [->] (0,%(cmin)f-1) -- (0,%(cmax)f) " % { 'cmin' : CMIN, 'cmax' : CMAX }
    str += "node [above] {$\\cost(\\numarcs)$} ;\n"
    # horizontal ticks
    for tick in np.arange(ZMIN,ZMAX+1,1) :
        data = { 'tick' : tick, 'horz' : horz_level }
        str += "\\draw (%(tick)d,%(horz)f-.1) -- (%(tick)f,%(horz)f+.1) " % data
        str += "node [below=5] {$%d$} ;\n" % tick
        
    # draw dashed C curve
    def drawpieces( ZZ, style ) :
        str = ''
        for z1,z2 in zip( ZZ[:-1], ZZ[1:] ) :
            c1 = Cf(z1)
            c2 = Cf(z2)
            
            data = { 'z1' : z1, 'z2' : z2, 'c1' : c1, 'c2' : c2, 'style' : style }
            str += "\\draw [%(style)s] (%(z1)f,%(c1)f) -- (%(z2)f,%(c2)f) ;\n" % data
            
        return str
            
    # draw dashed C, all of it
    str += drawpieces( range(ZMIN,ZMAX+1), 'dashed' )
    # draw C so far
    lastz = int( np.floor( zplus ) )
    str += drawpieces( range(ZMIN,lastz+1) + [ zplus ], 'thick' )
    
    # add local vertical line
    str += texvline( zplus, horz_level-1, max( CMAX, Cf(zplus) ) + 1, 'dashed' )
    # add optimal
    str += "\\draw node at (%(z)f,%(C)f) {$\\star$} ;\n" % { 'z' : ZOPT, 'C' : COPT }
    
    str += "\\end{tikzpicture}\n"
    return str
    








def position( address, roadmap, pos, length_attr='length' ) :
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
    X = [ x for x,y in points ]
    Y = [ y for x,y in points ]
    return X, Y




def SHOWMATCH( match, S, T, roadmap, pos, length_attr='length', ax=None ) :
    # crucial labels
    VERTEX = 'v'
    POINT_IN_S = 'S'
    POINT_IN_T = 'T'
    
    # plot the roadmap skeleton (undirected)
    if ax is None : ax = plt.gca()
    skeleton = nx.convert.convert_to_undirected( roadmap )
    nx.draw_networkx_nodes( skeleton, pos, ax=ax )
    nx.draw_networkx_edges( skeleton, pos=pos, ax=ax, edge_color='g', alpha=.15 )
    
    road_labels = { (u,v) : road
                   for u,v,road in roadmap.edges_iter( keys=True ) }
    nx.draw_networkx_edge_labels( skeleton, pos=pos, ax=ax, 
                                  edge_labels=road_labels )
    #nx.draw( skeleton, pos=pos, ax=ax, edge_color='g', alpha=.0 )
    ax.set_aspect('equal')
    
    """ Now, want to get edges with proper thickness """
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
            if score <= 0 : continue
            
            posu = other_pos[uu]
            posv = other_pos[vv]
            
            xu, yu = posu
            xv, yv = posv
            ax.plot( [xu,xv], [yu,yv], color='k', alpha=.6, linewidth=score )
        
    # plot the points on top, so visible
    # show S points in red
    positions = [ position(addr, roadmap, pos) for addr in S ]
    X, Y = pointsToXY( positions )
    ax.scatter( X, Y, color='r' )
    # show T points in blue
    positions = [ position(addr, roadmap, pos).tolist() for addr in T ]
    X, Y = pointsToXY( positions )
    ax.scatter( X, Y, color='b' )





if __name__ == '__main__' :
    plt.close('all')
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument( '--number', type=int, default=10 )
    #parser.add_argument( '--Hout', type=str, default='Hout.tex' )
    #parser.add_argument( '--Cout', type=str, default='Cout.tex' )
    args = parser.parse_args()
    
    N = args.number
    #N = 10
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
    
    # plot the network
    #nx.draw(graph, pos=pos)
    
    
    # requires bintrees, too.
    import setiptah.roadgeometry.roadmap_basic as ROAD
    import setiptah.roadbm.bm as roadbm
    import setiptah.roadgeometry.probability as roadprob
    
    uniform = roadprob.UniformDist( roadmap )
    
    
    M = 500
    addresses = [ uniform.sample() for i in xrange(M) ]
    positions = [ position(addr, roadmap, pos) for addr in addresses ]
    
    # <codecell>
    if False :
        X = [ x for x,y in positions ]
        Y = [ y for x,y in positions ]
        
        nx.draw(graph, pos=pos)
        ax = plt.gca()
        ax.scatter(X,Y)
        ax.set_aspect('equal')
    
    unpack = lambda addr : ( addr.road, addr.coord )
    SS = [ unpack( uniform.sample() ) for i in xrange(M) ]
    TT = [ unpack( uniform.sample() ) for i in xrange(M) ]
    
    match = roadbm.ROADSBIPARTITEMATCH( SS, TT, roadmap )
    #match = [ (i,i) for i in xrange(M) ]
    SHOWMATCH( match, SS, TT, roadmap, pos=pos )
