'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Nov 20, 2019

@author: Niels Lubbes
'''

import os

from surface_equivalence.class_se_tools import SETools

from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing

from linear_series.class_poly_ring import PolyRing
from linear_series.class_base_points import BasePointTree
from linear_series.class_linear_series import LinearSeries

from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_matrix
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_gcd


def usecase_B5():
    '''
    
    '''

    # f = 4e0-2e1-e2-e3 (class of hyperplane sections)
    # p = e0-e1         (class of conical family)
    # k = -3e0+e1+e2+e3 (canonical class)

    # basepoints in chart x0!=0;
    p1 = ( 0, 0 )
    p2 = ( 0, 1 )
    p3 = ( 1, 0 )


    # 0f+p = e0-e1
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p1, 1 )
    ls0_p1 = LinearSeries.get( [1], bp_tree )
    SETools.p( 'ls0_p1 =', ls0_p1 )

    # 1f-3p = e0+e1-e2-e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p2, 1 )
    bp = bp_tree.add( 'z', p3, 1 )
    ls1_m3 = LinearSeries.get( [1], bp_tree )
    SETools.p( 'ls1_m3 =', ls1_m3 )

    # 1f-2p = 2e0-e2-e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p2, 1 )
    bp = bp_tree.add( 'z', p3, 1 )
    ls1_m2 = LinearSeries.get( [2], bp_tree )
    SETools.p( 'ls1_m2 =', ls1_m2 )

    # 1f-1p = 3e0-e1-e2-e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p1, 1 )
    bp = bp_tree.add( 'z', p2, 1 )
    bp = bp_tree.add( 'z', p3, 1 )
    ls1_m1 = LinearSeries.get( [3], bp_tree )
    SETools.p( 'ls1_m1 =', ls1_m1 )

    # 1f-0p = 4e0-2e1-e2-e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p1, 2 )
    bp = bp_tree.add( 'z', p2, 1 )
    bp = bp_tree.add( 'z', p3, 1 )
    ls1_m0 = LinearSeries.get( [4], bp_tree )
    SETools.p( 'ls1_m0 =', ls1_m0 )

    # 2f-4p = 4e0-2e2-2e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p2, 2 )
    bp = bp_tree.add( 'z', p3, 2 )
    ls2_m4 = LinearSeries.get( [4], bp_tree )
    SETools.p( 'ls2_m4 =', ls2_m4 )

    # generators of graded ring
    s = ring( 'x' )
    t = ring( 'y' )
    u = ring( 'x+y-z' )
    v = ring( 'x*y' )
    w = ring( '(x+y-z)^2' )

    # find linear relation
    g = ring( ls2_m4.pol_lst )
    matg = sage_matrix( sage_QQ, SERing.get_matrix_P2( g ) )
    SETools.p( '\n' + str( matg ) )

    h = h0, h1, h2, h3, h4, h5, h6, h7, h8, h9 = [ v * v, v * w, w * w, s * u * v, s * u * w, t * u * v, t * u * w, s * s * u * u, s * t * u * u, t * t * u * u ]
    math = sage_matrix( sage_QQ, SERing.get_matrix_P2( h ) )
    SETools.p( '\n' + str( math ) )

    kerh = math.transpose().right_kernel().matrix()
    SETools.p( kerh )
    assert kerh * sage_vector( h ) == sage_vector( [0] )
    assert h1 - h8 == 0

    # construct map f
    f = ring( ls1_m0.pol_lst )
    SETools.p( 'f =', f )

    # construct map g
    # sage_Permutations(10).random_element().to_matrix().rows()
    matp = [( 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ), ( 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ), ( 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ),
            ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ), ( 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ), ( 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ),
            ( 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ), ( 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ), ( 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ),
            ( 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 )]
    matp = sage_matrix( matp )
    x, y, z, v, w = ring( 'x,y,z,v,w' )  # P2(x:y:z) and P1xP1(x:y;v:w)
    g = [comp.subs( {x:x * v, y:y * v, z:y * w} ) for comp in f]
    g = [comp / sage_gcd( g ) for comp in g]
    g = list( matp * sage_vector( g ) )
    SETools.p( 'g =', g )

    # determine and set basepoint tree
    ls = LinearSeries( g, PolyRing( 'x,y,v,w' ) )
    bp_tree = ls.get_bp_tree()
    SETools.p( 'bp_tree =', bp_tree )
    tree_211 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    tree_211.add( 'xw', ( 0, 0 ), 2 ).add( 't', ( 1, 0 ), 1 )
    tree_211.add( 'yv', ( 0, 1 ), 1 )

    # 1g+0q = 4a+2b-2e1-e2-e3
    g1m0 = LinearSeries.get( [4, 2], tree_211 )
    SETools.p( 'g1m0 =', g1m0 )

    # 1g-3q = (a+b-e1-e2-e3) + (b-e1)
    g1m3 = LinearSeries.get( [1, 2], tree_211 )
    SETools.p( 'g1m3 =', g1m3 )

    # 1g-2q = 2a+2b-2e1-e2-e3
    g1m2 = LinearSeries.get( [2, 2], tree_211 )
    SETools.p( 'g1m2 =', g1m2 )

    # 1g-1q = 3a+2b-2e1-e2-e3
    g1m1 = LinearSeries.get( [3, 2], tree_211 )
    SETools.p( 'g1m1 =', g1m1 )

    # 2g-4q = 4a+4b-4e1-2e2-2e3
    tree_422 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    tree_422.add( 'xw', ( 0, 0 ), 4 ).add( 't', ( 1, 0 ), 2 )
    tree_422.add( 'yv', ( 0, 1 ), 2 )
    g2m4 = LinearSeries.get( [4, 4], tree_422 )
    SETools.p( 'g2m4 =', g2m4 )

    # find linear relation
    s = ring( 'x' )
    t = ring( 'y' )
    u = ring( 'x*v^2+y*v^2-y*v*w' )
    v = ring( 'x*y*v^2' )
    w = ring( 'x^2*v^2+y^2*v^2-y^2*w^2' )

    xx = ring( 'x' )
    yy = ring( 'y' )
    ss = ring( 's' )
    tt = ring( 't' )
    uu = ring( 'u' )
    vv = ring( 'v' )
    ww = ring( 'w' )

    g = ring( g2m4.pol_lst )
    matg = sage_matrix( sage_QQ, SERing.get_matrix_P1xP1( [ comp.subs( {xx:ss, yy:tt, vv:uu, ww:vv} ) for comp in g] ) )
    SETools.p( '\n' + str( matg ) )

    h = h0, h1, h2, h3, h4, h5, h6, h7, h8, h9 = [ v * v, v * w, w * w, s * u * v, s * u * w, t * u * v, t * u * w, s * s * u * u, s * t * u * u, t * t * u * u ]
    math = sage_matrix( sage_QQ, SERing.get_matrix_P1xP1( [ comp.subs( {xx:ss, yy:tt, vv:uu, ww:vv} ) for comp in h] ) )
    SETools.p( '\n' + str( math ) )

    kerh = math.transpose().right_kernel().matrix()
    SETools.p( kerh )
    assert kerh * sage_vector( h ) == sage_vector( [0] )
    assert 2 * h0 + h1 - 2 * h3 - 2 * h5 + h8 == 0



if __name__ == '__main__':

    #  Debug output settings
    #
    mod_lst = []
    mod_lst += ['__main__.py']
    SETools.filter( mod_lst )  # output only from specified modules
    SETools.filter( None )  # print all verbose output, comment to disable.

    if 'OUTPUT_PATH' not in os.environ:
        os.environ['OUTPUT_PATH'] = './'
    # os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/mathematica/link/bin'

    SETools.start_timer()

    #########################################
    #                                       #
    # (Un)comment one or more use cases     #
    #                                       #
    #########################################

    usecase_B5()

    #########################################
    #                                       #
    # End of list of use case methods.      #
    #                                       #
    #########################################

    SETools.end_timer()
    SETools.p( 'The End' )









