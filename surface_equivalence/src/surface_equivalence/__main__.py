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
from surface_equivalence.sage_interface import sage_PolynomialRing
from surface_equivalence.sage_interface import sage_FractionField
from surface_equivalence.sage_interface import sage_ideal
from surface_equivalence.sage_interface import sage__eval
from surface_equivalence.sage_interface import sage_Compositions
from surface_equivalence.sage_interface import sage_solve
from surface_equivalence.sage_interface import sage_SR
from surface_equivalence.sage_interface import sage_lcm
from surface_equivalence.sage_interface import sage_denominator
from surface_equivalence.sage_interface import sage_identity_matrix


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
    f0p1 = SERing.conv( LinearSeries.get( [1], bp_tree ).pol_lst )
    SETools.p( 'f0p1 =', len( f0p1 ), f0p1 )

    # 1f-3p = e0+e1-e2-e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p2, 1 )
    bp = bp_tree.add( 'z', p3, 1 )
    f1m3 = SERing.conv( LinearSeries.get( [1], bp_tree ).pol_lst )
    SETools.p( 'f1m3 =', len( f1m3 ), f1m3 )

    # 1f-2p = 2e0-e2-e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p2, 1 )
    bp = bp_tree.add( 'z', p3, 1 )
    f1m2 = SERing.conv( LinearSeries.get( [2], bp_tree ).pol_lst )
    SETools.p( 'f1m2 =', len( f1m2 ), f1m2 )

    # 1f-1p = 3e0-e1-e2-e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p1, 1 )
    bp = bp_tree.add( 'z', p2, 1 )
    bp = bp_tree.add( 'z', p3, 1 )
    f1m1 = SERing.conv( LinearSeries.get( [3], bp_tree ).pol_lst )
    SETools.p( 'f1m1 =', len( f1m1 ), f1m1 )

    # 1f-0p = 4e0-2e1-e2-e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p1, 2 )
    bp = bp_tree.add( 'z', p2, 1 )
    bp = bp_tree.add( 'z', p3, 1 )
    f1m0 = SERing.conv( LinearSeries.get( [4], bp_tree ).pol_lst )
    SETools.p( 'f1m0 =', len( f1m0 ), f1m0 )

    # 2f-4p = 4e0-2e2-2e3
    bp_tree = BasePointTree()
    bp = bp_tree.add( 'z', p2, 2 )
    bp = bp_tree.add( 'z', p3, 2 )
    f2m4 = SERing.conv( LinearSeries.get( [4], bp_tree ).pol_lst )
    SETools.p( 'f2m4 =', len( f2m4 ), f2m4 )

    # generators of graded ring of f
    U = U0, U1, U2, U3, U4 = ring( 'x1' ), ring( 'x2' ), ring( 'x1+x2-x0' ), ring( 'x1*x2' ), ring( '(x1+x2-x0)^2' )

    # compute bidegree (2,d) in order to find a relation between the generators
    u = u0, u1, u2, u3, u4 = ring( 'u0,u1,u2,u3,u4' )
    SETools.p( 'Compare number of monomials of given bi-weigth with dimension predicted by the Riemann-Roch formula...' )
    for d in reversed( [-i for i in range( 8 )] ):
        w_lst = [( 0, 1 ), ( 0, 1 ), ( 1, -3 ), ( 1, -2 ), ( 1, -2 )]
        SETools.p( '\tweigth=', ( 2, d ), ',\t#monomials=', len( SERing.get_wmon_lst( u, w_lst, 2, d ) ), ',\tRR=', 29 + 5 * d )

    # template for generators of coordinate ring for weight (2,-1) and (1,0)
    T2m4 = ring( '[u3^2,u3*u4,u4^2,u0*u2*u3,u0*u2*u4,u1*u2*u3,u1*u2*u4,u0^2*u2^2,u0*u1*u2^2,u1^2*u2^2]' )
    T1m0 = ring( '[u1^2*u4,u1^2*u3,u1^3*u2,u0*u1*u4,u0*u1*u3,u0*u1^2*u2,u0^2*u4,u0^2*u3,u0^2*u1*u2,u0^3*u2]' )
    SETools.p( 'T2m4 =', T2m4 )
    SETools.p( 'T1m0 =', T1m0 )

    # find linear relation for f2m4
    a = a0, a1, a2, a3, a4, a5, a6, a7, a8, a9 = [elt.subs( {u[i]:U[i] for i in range( 5 )} ) for elt in T2m4 ]
    mata = sage_matrix( sage_QQ, SERing.get_matrix_P2( a ) )
    kera = mata.transpose().right_kernel().matrix()
    SETools.p( 'kera =', kera )
    assert kera * sage_vector( a ) == sage_vector( [0] )
    assert a1 - a8 == 0

    # construct map gg from ff
    # sage_Permutations(10).random_element().to_matrix().rows()
    ff = f1m0
    matp = [( 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ), ( 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ), ( 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ),
            ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ), ( 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 ), ( 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 ),
            ( 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ), ( 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ), ( 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ),
            ( 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 )]
    matp = sage_matrix( matp )
    x0, x1, x2, y0, y1, y2, y3 = ring( 'x0,x1,x2,y0,y1,y2,y3' )  # P2(x0:x1:x2) and P1xP1(y0:y1;y2:y3)
    gg = [comp.subs( {x0:y1 * y3, x1:y0 * y2, x2 :y1 * y2} ) for comp in ff]
    gg = list( matp * sage_vector( gg ) )
    gcd_gg = sage_gcd( gg )
    gg = [comp / gcd_gg for comp in gg]
    SETools.p( 'gcd_gg =', gcd_gg )
    SETools.p( 'ff     =', len( ff ), ff )
    SETools.p( 'gg     =', len( gg ), gg )


    # determine and set basepoint tree
    ls = LinearSeries( SERing.conv( gg ), PolyRing( 'x,y,v,w' ) )
    bp_tree = ls.get_bp_tree()
    SETools.p( 'bp_tree(gg) =', bp_tree )
    tree_211 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    tree_211.add( 'xw', ( 0, 0 ), 2 ).add( 't', ( 1, 0 ), 1 )
    tree_211.add( 'yv', ( 0, 1 ), 1 )

    # 1g+0q = 4a+2b-2e1-e2-e3
    g1m0 = SERing.conv( LinearSeries.get( [4, 2], tree_211 ).pol_lst )
    SETools.p( 'g1m0 =', len( g1m0 ), g1m0 )

    # 1g-3q = (a+b-e1-e2-e3) + (b-e1)
    g1m3 = SERing.conv( LinearSeries.get( [1, 2], tree_211 ).pol_lst )
    SETools.p( 'g1m3 =', len( g1m3 ), g1m3 )

    # 1g-2q = 2a+2b-2e1-e2-e3
    g1m2 = SERing.conv( LinearSeries.get( [2, 2], tree_211 ).pol_lst )
    SETools.p( 'g1m2 =', len( g1m2 ), g1m2 )

    # 1g-1q = 3a+2b-2e1-e2-e3
    g1m1 = SERing.conv( LinearSeries.get( [3, 2], tree_211 ).pol_lst )
    SETools.p( 'g1m1 =', len( g1m1 ), g1m1 )

    # 2g-4q = 4a+4b-4e1-2e2-2e3
    tree_422 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    tree_422.add( 'xw', ( 0, 0 ), 4 ).add( 't', ( 1, 0 ), 2 )
    tree_422.add( 'yv', ( 0, 1 ), 2 )
    g2m4 = SERing.conv( LinearSeries.get( [4, 4], tree_422 ).pol_lst )
    SETools.p( 'g2m4 =', len( g2m4 ), g2m4 )

    # set generators of graded coordinate ring of g
    V = V0, V1, V2, V3, V4 = ring( 'y0' ), ring( 'y1' ), ring( 'y0*y2^2+y1*y2^2-y1*y2*y3' ), ring( 'y0*y1*y2^2' ), ring( 'y0^2*y2^2+y1^2*y2^2-y1^2*y3^2' )

    # find linear relation for g2m4
    b = b0, b1, b2, b3, b4, b5, b6, b7, b8, b9 = [elt.subs( {u[i]:V[i] for i in range( 5 )} ) for elt in T2m4 ]
    matb = sage_matrix( sage_QQ, SERing.get_matrix_P1xP1( b ) )
    kerb = matb.transpose().right_kernel().matrix()
    SETools.p( 'kerb =', kerb )
    assert kerb * sage_vector( b ) == sage_vector( [0] )
    assert 2 * b0 + b1 - 2 * b3 - 2 * b5 + b8 == 0

    # compute inverse of G
    G = [ elt.subs( {u[i]:V[i] for i in range( 5 )} ) for elt in T1m0]
    z = ring( 'z0,z1,z2,z3,z4,z5,z6,z7,z8,z9' )
    t = ring( 't' )
    id = [ G[i] * z[0] - z[i] * G[0] for i in range( 10 ) ] + [t * G[0] - 1]
    I01 = sage_ideal( id ).elimination_ideal( [t, y2, y3 ] ).gens()
    I23 = sage_ideal( id ).elimination_ideal( [t, y0, y1 ] ).gens()
    I01 = [ elt for elt in I01 if elt.degree( y0 ) == 1 and elt.degree( y1 ) == 1 ][0]
    I23 = [ elt for elt in I23 if elt.degree( y2 ) == 1 and elt.degree( y3 ) == 1 ][0]
    Q0 = I01.coefficient( y1 )
    Q1 = -I01.coefficient( y0 )
    Q2 = I23.coefficient( y3 )
    Q3 = -I23.coefficient( y2 )
    Q = [Q0, Q1, Q2, Q3]
    SETools.p( 'Q =', Q )
    # [-j, -i, -i, -g - 2*h + i + j]

    QoG = [q.subs( { z[i]:G[i] for i in range( 10 ) } ) for q in Q ]
    gcd01 = sage_gcd( QoG[0], QoG[1] )
    gcd23 = sage_gcd( QoG[2], QoG[3] )
    QoG = [QoG[0] / gcd01, QoG[1] / gcd01, QoG[2] / gcd23, QoG[3] / gcd23]
    SETools.p( 'QoG =', QoG )
    assert QoG == [y0, y1, y2, y3]

    # compose F with projective isomorphism P
    c = c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12 = [ring( 'c' + str( i ) ) for i in range( 13 )]
    dctZ = { u[i]:z[i] for i in range( 5 )}
    dctP = { z[0]:c0 * u0 + c1 * u1,
             z[1]:c2 * u0 + c3 * u1,
             z[2]:c4 * u2,
             z[3]:c5 * u3 + c6 * u4 + c7 * u0 * u2 + c8 * u1 * u2,
             z[4]:c9 * u3 + c10 * u4 + c11 * u0 * u2 + c12 * u1 * u2 }
    PoF = [ comp.subs( dctZ ).subs( dctP ) for comp in T1m0]
    PoF = [ comp.subs( {u[i]:U[i] for i in range( 5 )} ) for comp in PoF]
    PoF = [ comp / sage_gcd( PoF ) for comp in PoF]
    SETools.p( 'PoF =', len( PoF ), PoF )

    # compose PoF with Q
    QoPoF = [ comp.subs( {z[i]:PoF[i] for i in range( 10 )} ) for comp in Q]
    gcd01 = sage_gcd( QoPoF[0], QoPoF[1] )
    gcd23 = sage_gcd( QoPoF[2], QoPoF[3] )
    QoPoF = [QoPoF[0] / gcd01, QoPoF[1] / gcd01, QoPoF[2] / gcd23, QoPoF[3] / gcd23]
    SETools.p( 'QoPoF =', len( QoPoF ), QoPoF )

    # create list of equations for the ci
    b = T2m4
    rel_g4m2 = 2 * b[0] + b[1] - 2 * b[3] - 2 * b[5] + b[8]
    SETools.p( 'rel_g4m2 =', rel_g4m2 )
    rel_g4m2 = rel_g4m2.subs( dctZ ).subs( dctP ).subs( {u[i]:U[i] for i in range( 5 )} )
    SETools.p( 'rel_g4m2 =', rel_g4m2 )
    rel_lst = []
    x = ring( '[x0,x1,x2]' )
    for exp in sage_Compositions( 4 + 3, length = 3 ):
        rel_lst += [rel_g4m2.coefficient( {x[i]:exp[i] - 1 for i in range( 3 )} )]
    SETools.p( 'rel_lst =', len( rel_lst ), rel_lst )
    t = ring( 't' )
    rel_lst += [ ( c0 * c3 - c1 * c2 ) * c4 * ( c5 * c10 - c9 * c6 ) * t - 1 ]

    # solve for ci
    prime_lst = sage_ideal( rel_lst ).elimination_ideal( t ).primary_decomposition()
    SETools.p( 'prime_lst =', len( prime_lst ) )
    for prime in prime_lst:
        SETools.p( '\t', prime.gens() )

    # SETools.p( '>', sage_ideal( rel_lst ).elimination_ideal( t ).triangular_decomposition() )


    # put solutions in dictionary form
    for gen_lst in [prime.gens() for prime in prime_lst]:
        sol_dct = sage_solve( [sage_SR( gen ) for gen in gen_lst], [sage_SR( elt ) for elt in c], solution_dict = True )
        SETools.p( '\t sol_dct =', sol_dct )
        assert len( sol_dct ) == 1
    prime_lst2 = []
    prime_lst2 += [prime_lst[0].gens() + [c0 - 1, c4 - 1]]
    prime_lst2 += [prime_lst[1].gens() + [c1 - 1, c4 - 1]]
    prime_lst2 += [prime_lst[2].gens() + [c1 - 1, c4 - 1]]
    prime_lst2 += [prime_lst[3].gens() + [c0 - 1, c4 - 1]]
    SETools.p( 'Added equations to prime_lst to simplify solutions:' )
    for prime in prime_lst2:
        SETools.p( '\t', prime )
    SETools.p( 'Simplified solutions:' )
    for gen_lst in prime_lst2:
        sol_dct = sage_solve( [sage_SR( gen ) for gen in gen_lst], [sage_SR( elt ) for elt in c], solution_dict = True )
        SETools.p( '\t sol_dct =', sol_dct )
        assert len( sol_dct ) == 1

    r0, r1 = ring( 'r0,r1' )
    sol0 = {c0:1, c1:0, c2:0, c3:-r0 * r1, c4:1, c5:0, c6:r0, c7:0, c8:0, c9:r1, c10:-2 * r0, c11:2, c12:-2 * r0 * r1}
    sol1 = {c0:0, c1:1, c2:-r0 * r1, c3:0, c4:1, c5:0, c6:r0, c7:0, c8:0, c9:r1, c10:-2 * r0, c11:-2 * r0 * r1, c12:2}
    sol2 = {c0:0, c1:1, c2:-r0 * r1, c3:0, c4:1, c5:r0, c6:0, c7:0, c8:0, c9:-2 * r0, c10:r1, c11:-2 * r0 * r1, c12:2}
    sol3 = {c0:1, c1:0, c2:0, c3:-r0 * r1, c4:1, c5:r0, c6:0, c7:0, c8:0, c9:-2 * r0, c10:r1, c11:2, c12:-2 * r0 * r1}
    sol_lst = [sol0, sol1, sol2, sol3]

    SETools.p( 'Simplified solutions by hand:' )
    for sol in sol_lst:
        SETools.p( '\t', sol )


    y = ring( '[y0,y1,y2,y3]' )
    gr_lst = []
    SETools.p( 'Computing (gg o r) for each sol in sol_lst...' )
    for sol in sol_lst:
        gr = [ comp.subs( {y[i]:QoPoF[i] for i in range( 4 )} ).subs( sol ) for comp in gg ]
        SETools.p( '\t gr =', gr )
        gcd_gr = sage_gcd( gr )
        SETools.p( '\t\tgcd_gr =', gcd_gr )
        gr_lst += [[ comp / gcd_gr for comp in gr ]]
    SETools.p( 'gr_lst =', len( gr_lst ) )
    for gr in gr_lst:
        SETools.p( '\t gr =', gr )

    # load("/home/niels/Desktop/n/src/git/surface_equivalence/surface_equivalence/src/surface_equivalence/se_tools")['gr_lst']


    # get coefficient matrix ff
    mff = SERing.get_matrix_P2( ff )
    kff = mff.right_kernel_matrix().T
    SETools.p( 'mff =', mff.dimensions(), list( mff ) )
    SETools.p( 'kff =', kff.dimensions(), list( kff ) )
    assert ( mff * kff ).is_zero()

    # get implicit equations for image of gg
    z = ring( 'z0,z1,z2,z3,z4,z5,z6,z7,z8,z9' )
    y = ring( 'y0,y1,y2,y3' )
    igg = SERing.R.ideal( [ z[i] - gg[i] for i in range( 10 )  ] ).elimination_ideal( [y[i] for i in range( 4 )] )
    SETools.p( 'igg =', igg )

    # Compute isomorphisms for each gr
    SETools.p( 'Compute projective isomorphism for each gr in gr_lst:' )
    for gr in gr_lst:

        mgr = SERing.get_matrix_P2( gr )
        mgk = mgr * kff
        assert mgk.is_zero()  # because the surfaces in P^9 are linearly normal

        Ef = sage_matrix( sage_QQ, mff.rows() + kff.T.rows() )
        Egr = sage_matrix( mgr.rows() + kff.T.rows() )
        UpI = Egr * Ef.inverse()
        assert ( UpI.submatrix( 10, 10 ) - sage_identity_matrix( 5 ) ).is_zero()

        U = UpI.submatrix( 0, 0, 10, 10 )
        SETools.p( '\tU =', U.dimensions(), list( U ) )

        # check if the answer is correct
        Uff = list( U * sage_vector( ff ) )
        iggs = igg.subs( {z[i]:Uff[i] for i in range( 10 )} )
        assert iggs.is_zero()



def usecase_invert_map():
    '''
    Inversion of maps using Groebner basis.
    '''

#     RP = sage_PolynomialRing( sage_QQ, 'y0,y1,y2,y3', order = 'deglex' )
#     y0, y1, y2, y3 = RP.gens()
#     RH = sage_PolynomialRing( sage_FractionField( RP ), 'x0,x1,x2,x3,x4', order = 'deglex' )
#     x0, x1, x2, x3, x4 = RH.gens()
#
#     X = [x1 ** 2 + x2 ** 2 + x3 ** 2 + x4 ** 2 - x0 ** 2]
#     f = [y0 - ( x0 - x4 ), y1 - x1, y2 - x2, y3 - x3]
#     G1 = sage_ideal( X + f ).groebner_basis()
#     SETools.p( 'G1 =', G1 )
#
#
#     RP = sage_PolynomialRing( sage_QQ, 'y1,y2,y3,y4', order = 'deglex' )
#     y1, y2, y3, y4 = RP.gens()
#     S = RP.quotient_ring( [ y1 ** 2 + y2 ** 2 + y3 ** 2 + y4 ** 2 - 1] )
#     RH = sage_PolynomialRing( sage_FractionField( S ), 'x1,x2,x3,x4', order = 'deglex' )
#     x1, x2, x3, x4 = RH.gens()
#
#     d = x1 ** 2 + x2 ** 2 + x3 ** 2
#     g = [ ( d + 1 ) * x4 - 1,
#           2 * x1 * x4 - y1,
#           2 * x2 * x4 - y2,
#           2 * x3 * x4 - y3,
#           ( d - 1 ) * x4 - y4]
#     G2 = sage_ideal( g ).groebner_basis()
#     SETools.p( 'G2 =', str( G2 ).replace( 'bar', '' ) )


    R = sage_PolynomialRing( sage_QQ, 'x0,x1,x2,x3,y0,y1,y2,y3,y4,t', order = 'deglex' )
    x0, x1, x2, x3, y0, y1, y2, y3, y4, t = R.gens()

    d = x1 ** 2 + x2 ** 2 + x3 ** 2
    f0 = d + x0 ** 2
    f1 = 2 * x0 * x1
    f2 = 2 * x0 * x2
    f3 = 2 * x0 * x3
    f4 = d - x0 ** 2
    eq = -y0 ** 2 + y1 ** 2 + y2 ** 2 + y3 ** 2 + y4 ** 2

    g = [ f1 * y0 - y1 * f0,
          f2 * y0 - y2 * f0,
          f3 * y0 - y3 * f0,
          f4 * y0 - y4 * f0,
          t * f0 - 1 ]

    SETools.p( sage_ideal( g ).elimination_ideal( [ t ] ).gens() )
    Ix1 = sage_ideal( g ).elimination_ideal( [t, x2, x3] ).gens()
    Ix2 = sage_ideal( g ).elimination_ideal( [t, x1, x3] ).gens()
    Ix3 = sage_ideal( g ).elimination_ideal( [t, x1, x2] ).gens()

    SETools.p( 'Ix1 =', Ix1 )
    SETools.p( 'Ix2 =', Ix2 )
    SETools.p( 'Ix3 =', Ix3 )

    Ix1 = [ elt for elt in Ix1 if elt.degree( x0 ) == 1 and elt.degree( x1 ) == 1 ][0]
    Ix2 = [ elt for elt in Ix2 if elt.degree( x0 ) == 1 and elt.degree( x2 ) == 1 ][0]
    Ix3 = [ elt for elt in Ix3 if elt.degree( x0 ) == 1 and elt.degree( x3 ) == 1 ][0]

    SETools.p( 'Ix1 =', Ix1 )
    SETools.p( 'Ix2 =', Ix2 )
    SETools.p( 'Ix3 =', Ix3 )

    g0 = Ix1.coefficient( x1 )
    g1 = -Ix1.coefficient( x0 )
    g2 = -Ix2.coefficient( x0 )
    g3 = -Ix3.coefficient( x0 )

    SETools.p( 'g = f^(-1) =', [g0, g1, g2, g3] )




if __name__ == '__main__':

    #  Debug output settings
    #
    mod_lst = []
    mod_lst += ['__main__.py']
    SETools.filter( mod_lst )  # output only from specified modules
    # SETools.filter( None )  # print all verbose output, comment to disable.

    if 'OUTPUT_PATH' not in os.environ:
        os.environ['OUTPUT_PATH'] = './'
    os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/mathematica/link/bin'

    SETools.start_timer()

    #########################################
    #                                       #
    # (Un)comment one or more use cases     #
    #                                       #
    #########################################

    usecase_B5()
    # usecase_invert_map()

    #########################################
    #                                       #
    # End of list of use case methods.      #
    #                                       #
    #########################################

    SETools.end_timer()
    SETools.p( 'The End' )









