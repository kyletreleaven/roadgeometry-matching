

import itertools

import numpy as np
import bintrees

import networkx as nx

""" my dependencies """
import setiptah.roadgeometry.roadmap_basic as ROAD
import setiptah.roadgeometry.astar_basic as ASTAR

import setiptah.bpmatch.roadmaps as roadbm




def drawHeightFunction( S, T, length ) :
    roadnet = nx.MultiDiGraph()
    roadnet.add_edge(0,1,'line', length=length )
    
    SS = ( ('line',s) for s in S )
    TT = ( ('line',t) for t in T )
    
    segments = roadbm.SEGMENTS( SS, TT, roadnet )
    
    
    return segments




if __name__ == '__main__' :
    YMIN = -4.
    YMAX = 4.
    WIDTH = YMAX - YMIN
    if False :
        zr = 0
        X = [ YMIN + WIDTH * np.random.rand() for i in xrange(5) ]
        O = [ YMIN + WIDTH * np.random.rand() for i in xrange(3) ]
    else :
        zr = 1
        #X = [ -3.5, 2, 2.75 ]
        #O = [ -3., -2., 1., 2.25, 3. ]
        X = [ -2.5, 2, ]
        O = [ -2., -1., 1., 3. ]
        
    segment = drawHeightFunction( X, O, YMAX )['line']
    intervals = roadbm.INTERVALS( segment )
    
    """ to tikz """
    str = ''
    str += "\\begin{tikzpicture}\n"
    
    # draw the axis
    # horizontal
    str += "\\draw [->] (%f,0) -- (%f,0) node [right] {$\\coordvar$} ;\n" % ( YMIN, YMAX+.5 )
    str += "\\draw [thick] (%(ymax)f,-.15) -- (%(ymax)f,.15) " % { 'ymax' : YMAX }
    str += "node [above] {$\\roadlen_\\roadvar$} ;\n"
    # vertical
    fmin = np.floor( min( intervals ) )
    fmax = np.ceil( max( intervals ) )
    str += "\\draw [->] (%(ymin)f,%(fmin)f) -- (%(ymin)f,%(fmax)f) " % { 'ymin' : YMIN, 'fmin' : fmin-.25, 'fmax' : fmax+.25 }
    str += "node [above] {$\\postcumarcs(y;\\roadvar) = \\cumarcs(y;\\roadvar) + \\numarcs_\\roadvar$} ;\n"
    # vertical ticks
    for tick in np.arange(fmin,fmax+1,1) :
        str += "\\draw (%(tR)f,%(tick)d) -- (%(tL)f,%(tick)d) " % { 'tR' : YMIN+.1, 'tL' : YMIN-.1, 'tick' : tick }
        str += "node [left=5] {$%d$} ;\n" % tick
        
    # place the X's and O's
    for y, item in segment.items() :
        phases = [ (item.P, '$\\times$' ), (item.Q, '$\\circ$') ]
        for Y, mark in phases :
            for i in Y : str += "\\draw (%f,0) node {%s} ;\n" % ( y, mark )
            
    
    
    # place the levels
    for f, ranges in intervals.items() :
        for a,b in ranges :
            label = ''
            if a == '-' :
                a = YMIN
                label = "node [%s] {$\\numarcs_\\roadvar$}"
                if True :       # zr >= 0?
                    label = label % "above"
                else :
                    label = label % "below"
            if b == '+' :
                b = YMAX
                str += "\\draw (%f,%f) node [right] {$\\numarcs_\\roadvar + \\surplus_\\roadvar$} ;\n" % (b,f)
                            
            levelstr = "\\draw [thick] (%(yl)f,%(z)f) -- %(extra)s (%(yr)f,%(z)f) ;\n"
            DATA = { 'yl' : a, 'yr' : b, 'z' : f, 'extra' : label }
            str += levelstr % DATA
            # place the shades
            str += "\\path [fill=black,opacity=.2] (%(yl)f,0) rectangle (%(yr)f,%(z)f) ;\n" % DATA
            
    



    str += "\\end{tikzpicture}\n"
    
    f = open( 'Hout.tex', 'w' )
    f.write( str )
    f.close()




