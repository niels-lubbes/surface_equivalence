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
from linear_series.class_ls_tools import LSTools

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
from surface_equivalence.sage_interface import sage_diff


def usecase_B1():
    '''
    We compute the automorphisms of a Roman surface from a set of 
    compatible reparametrizations.        
    If the domain the projective plane, then a parametrization of a 
    Roman surface is basepoint free. 
    The compatible reparametrizations are in this case linear 
    automorphisms of the projective plane.
    '''
    y = ring( 'y0,y1,y2,y3' )
    x = ring( 'x0,x1,x2' )
    c = ring( 'c0,c1,c2,c3,c4,c5,c6,c7,c8' )

    # parametrizations f and g of Roman surface
    f = ring( '[x0^2+x1^2+x2^2,x0*x1,x0*x2,x1*x2]' )
    g = f

    # compatible reparametrizations are linear and indexed by c
    r = {x[0]:c[0] * y[0] + c[1] * y[1] + c[2] * y[2],
         x[1]:c[3] * y[0] + c[4] * y[1] + c[5] * y[2],
         x[2]:c[6] * y[0] + c[7] * y[1] + c[8] * y[2]}

    # compute kernel and coefficient matrix of f
    Mf = SERing.get_matrix_P2( f )
    Kf = Mf.right_kernel_matrix().T
    assert ( Mf * Kf ).is_zero()

    # compute the coefficient matrix of g composed with r
    gr = [ comp.subs( r ) for comp in g ]
    assert sage_gcd( gr ) == 1
    assert SERing.get_degree( gr, 'y0,y1,y2' ) == SERing.get_degree( f )
    Mgr = SERing.get_matrix_P2( gr, 'y0,y1,y2' )

    # output to screen
    SETools.p( 'f       =', f )
    SETools.p( 'g       =', g )
    SETools.p( 'r       =', r )
    SETools.p( 'Mf      =', Mf.dimensions(), list( Mf ) )
    SETools.p( 'Mgr     =', Mgr.dimensions(), list( Mgr ), '\n' + str( Mgr ) )

    # compute c such that Mgr*Kf==0
    ec_lst = ( Mgr * Kf ).list() + [  sage_matrix( SERing.R, 3, 3, c ).det() * ring( 't' ) - 1 ]
    pc_lst = sage_ideal( ec_lst ).elimination_ideal( ring( 't' ) ).primary_decomposition()
    SETools.p( 'sol_lst = ' )
    sol_lst = []
    for pc in pc_lst:
        s_lst = list( reversed( sorted( pc.gens() ) ) )
        s_dct = {}
        for s in s_lst:
            if s in c:
                s_dct.update( {s:0} )
            for ca in c:
                for cb in c:
                    if s == ca - cb and ca > cb: s_dct.update( {ca:cb} )
                    if s == ca - cb and ca < cb: s_dct.update( {cb:ca} )
                    if s == ca + cb and ca > cb: s_dct.update( {ca:-cb} )
                    if s == ca + cb and ca < cb: s_dct.update( {cb:-ca} )
        sol_lst += [s_dct]
        SETools.p( '\t\t', s_lst, '-->', s_dct )

    # compute implicit equation for image Y of g
    eqg = sage_ideal( [y[i] - g[i] for i in range( 4 )] ).elimination_ideal( x ).gens()
    SETools.p( 'eqg =', eqg )
    SETools.p( '    =', str( eqg.subs( {y[0]:1} ) ).replace( 'y1', 'x' ).replace( 'y2', 'y' ).replace( 'y3', 'z' ) )

    # computing U and test each sol in sol_lst
    SETools.p( 'Testing each sol in sol_lst...' )
    for sol in sol_lst:
        # compute the projective isomorphism in terms of parametrized matrix U
        Ef = sage_matrix( sage_QQ, list( Mf ) + list( Kf.T ) )
        Egr = sage_matrix( list( Mgr.subs( sol ) ) + list( Kf.T ) )
        UpI = Egr * ~Ef
        assert ( UpI.submatrix( 4, 4 ) - sage_identity_matrix( 2 ) ).is_zero()
        U = UpI.submatrix( 0, 0, 4, 4 )
        U = U / sage_gcd( U.list() )
        assert U.dimensions() == ( 4, 4 )

        # verify whether U*f is a parametrization for Y for all (c0,...,c7)
        Uf = list( U * sage_vector( f ) )
        eqg_sub = [ eq.subs( {y[i]:Uf[i] for i in range( 4 )} ) for eq in eqg ]
        assert eqg_sub == [0]

        # output U with corresponding solution
        SETools.p( '\t U =', list( U ), ', sol =', sol )


def usecase_B2_helper_bp( gr ):
    '''
    This is a helper method for usecase_B2().
    
    We return equations in c that give conditions so that 
    gr has the same basepoints as f.

    The basepoints and infinitely near basepoints of f are as follows:
    
        { 4, <<x^2*z^6, x^5*y^2*z, x^3*y^5, x*y^2*z^5 + 2*y^3*z^5>>, QQ[x, y, z] }
        chart=z, depth=0, mult=2, sol=(0, 0), { 4, <<x^2, x^5*y^2, x^3*y^5, x*y^2 + 2*y^3>>, QQ[x, y] }
            chart=t, depth=1, mult=1, sol=(0, 0), { 4, <<x^2, x^5*y^5, x^3*y^6, x*y + 2*y>>, QQ[x, y] }
                chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x, x^9*y^5, x^8*y^6, x*y + 2*y>>, QQ[x, y] }
        chart=x, depth=0, mult=3, sol=(0, 0), { 4, <<z^6, y^2*z, y^5, 2*y^3*z^5 + y^2*z^5>>, QQ[y, z] }
            chart=t, depth=1, mult=2, sol=(0, 0), { 4, <<z^3, y^2, y^5*z^2, 2*y^3*z^5 + y^2*z^4>>, QQ[y, z] }
                chart=t, depth=2, mult=1, sol=(0, 0), { 4, <<z, y^2, y^5*z^5, 2*y^3*z^6 + y^2*z^4>>, QQ[y, z] }
                    chart=s, depth=3, mult=1, sol=(0, 0), { 4, <<z, y, y^9*z^5, 2*y^8*z^6 + y^5*z^4>>, QQ[y, z] }
            chart=s, depth=1, mult=1, sol=(0, 0), { 4, <<y^3*z^6, z, y^2, 2*y^5*z^5 + y^4*z^5>>, QQ[y, z] }
                chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<y^8*z^6, z, y, 2*y^9*z^5 + y^8*z^5>>, QQ[y, z] }
        chart=y, depth=0, mult=3, sol=(0, 0), { 4, <<x^2*z^6, x^5*z, x^3, x*z^5 + 2*z^5>>, QQ[x, z] }
            chart=t, depth=1, mult=2, sol=(0, 0), { 4, <<x^2*z^5, x^5*z^3, x^3, x*z^3 + 2*z^2>>, QQ[x, z] }
                chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x^5*z^5, x^6*z^3, x, x^2*z^3 + 2*z^2>>, QQ[x, z] }
                    chart=t, depth=3, mult=1, sol=(0, 0), { 4, <<x^5*z^9, x^6*z^8, x, x^2*z^4 + 2*z>>, QQ[x, z] }
    
    '''
    x = [ring( 'x' + str( i ) ) for i in range( 3 )]

    # initialize the list to be returned
    eqn_lst = []

    #
    #     chart=z, depth=0, mult=2, sol=(0, 0), { 4, <<x^2, x^5*y^2, x^3*y^5, x*y^2 + 2*y^3>>, QQ[x, y] }
    #         chart=t, depth=1, mult=1, sol=(0, 0), { 4, <<x^2, x^5*y^5, x^3*y^6, x*y + 2*y>>, QQ[x, y] }
    #             chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x, x^9*y^5, x^8*y^6, x*y + 2*y>>, QQ[x, y] }
    #
    gr0 = [ comp.subs( {x[0]:1} ) for comp in gr ]
    gr0t = [ comp.subs( {x[1]:x[1] * x[2]} ).quo_rem( x[2] ** 2 )[0] for comp in gr0 ]
    gr0ts = [ comp.subs( {x[2]:x[1] * x[2]} ).quo_rem( x[1] )[0] for comp in gr0t ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 )]:
        eqn_lst += [ sage_diff( comp, x[1], a, x[2], b ).subs( {x[1]:0, x[2]:0} ) for comp in gr0 ]
    eqn_lst += [ comp.subs( {x[1]:0, x[2]:0} ) for comp in gr0t ]
    eqn_lst += [ comp.subs( {x[1]:0, x[2]:0} ) for comp in gr0ts ]

    #
    #     chart=x, depth=0, mult=3, sol=(0, 0), { 4, <<z^6, y^2*z, y^5, 2*y^3*z^5 + y^2*z^5>>, QQ[y, z] }
    #         chart=t, depth=1, mult=2, sol=(0, 0), { 4, <<z^3, y^2, y^5*z^2, 2*y^3*z^5 + y^2*z^4>>, QQ[y, z] }
    #             chart=t, depth=2, mult=1, sol=(0, 0), { 4, <<z, y^2, y^5*z^5, 2*y^3*z^6 + y^2*z^4>>, QQ[y, z] }
    #                 chart=s, depth=3, mult=1, sol=(0, 0), { 4, <<z, y, y^9*z^5, 2*y^8*z^6 + y^5*z^4>>, QQ[y, z] }
    #         chart=s, depth=1, mult=1, sol=(0, 0), { 4, <<y^3*z^6, z, y^2, 2*y^5*z^5 + y^4*z^5>>, QQ[y, z] }
    #             chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<y^8*z^6, z, y, 2*y^9*z^5 + y^8*z^5>>, QQ[y, z] }
    #
    gr1 = [ comp.subs( {x[1]:1} ) for comp in gr ]
    gr1t = [ comp.subs( {x[0]:x[0] * x[2]} ).quo_rem( x[2] ** 3 )[0] for comp in gr1 ]
    gr1tt = [ comp.subs( {x[0]:x[0] * x[2]} ).quo_rem( x[2] ** 2 )[0] for comp in gr1t ]
    gr1tts = [ comp.subs( {x[2]:x[0] * x[2]} ).quo_rem( x[0] ** 1 )[0] for comp in gr1tt ]
    gr1s = [ comp.subs( {x[2]:x[0] * x[2]} ).quo_rem( x[0] ** 3 )[0] for comp in gr1 ]
    gr1ss = [ comp.subs( {x[2]:x[0] * x[2]} ).quo_rem( x[0] ** 1 )[0] for comp in gr1s ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 ), ( 2, 0 ), ( 1, 1 ), ( 0, 2 )]:
        eqn_lst += [ sage_diff( comp, x[0], a, x[2], b ).subs( {x[0]:0, x[2]:0} ) for comp in gr1 ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 )]:
        eqn_lst += [ sage_diff( comp, x[0], a, x[2], b ).subs( {x[0]:0, x[2]:0} ) for comp in gr1t ]
    eqn_lst += [ comp.subs( {x[0]:0, x[2]:0} ) for comp in gr1tt ]
    eqn_lst += [ comp.subs( {x[0]:0, x[2]:0} ) for comp in gr1tts ]
    eqn_lst += [ comp.subs( {x[0]:0, x[2]:0} ) for comp in gr1s ]
    eqn_lst += [ comp.subs( {x[0]:0, x[2]:0} ) for comp in gr1ss ]

    #
    #     chart=y, depth=0, mult=3, sol=(0, 0), { 4, <<x^2*z^6, x^5*z, x^3, x*z^5 + 2*z^5>>, QQ[x, z] }
    #         chart=t, depth=1, mult=2, sol=(0, 0), { 4, <<x^2*z^5, x^5*z^3, x^3, x*z^3 + 2*z^2>>, QQ[x, z] }
    #             chart=s, depth=2, mult=1, sol=(0, 0), { 4, <<x^5*z^5, x^6*z^3, x, x^2*z^3 + 2*z^2>>, QQ[x, z] }
    #                 chart=t, depth=3, mult=1, sol=(0, 0), { 4, <<x^5*z^9, x^6*z^8, x, x^2*z^4 + 2*z>>, QQ[x, z] }
    #
    gr2 = [ comp.subs( {x[2]:1} ) for comp in gr ]
    gr2t = [ comp.subs( {x[0]:x[0] * x[1]} ).quo_rem( x[1] ** 3 )[0] for comp in gr2 ]
    gr2ts = [ comp.subs( {x[1]:x[0] * x[1]} ).quo_rem( x[0] ** 2 )[0] for comp in gr2t ]
    gr2tst = [ comp.subs( {x[0]:x[0] * x[1]} ).quo_rem( x[1] ** 1 )[0] for comp in gr2ts ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 ), ( 2, 0 ), ( 1, 1 ), ( 0, 2 )]:
        eqn_lst += [ sage_diff( comp, x[0], a, x[1], b ).subs( {x[0]:0, x[1]:0} ) for comp in gr2 ]
    for a, b in [( 0, 0 ), ( 0, 1 ), ( 1, 0 )]:
        eqn_lst += [ sage_diff( comp, x[0], a, x[1], b ).subs( {x[0]:0, x[1]:0} ) for comp in gr2t ]
    eqn_lst += [ comp.subs( {x[0]:0, x[1]:0} ) for comp in gr2ts ]
    eqn_lst += [ comp.subs( {x[0]:0, x[1]:0} ) for comp in gr2tst ]

    eqn_lst = sorted( list( set( eqn_lst ) ) )

    assert 'x0' not in str( eqn_lst )
    assert 'x1' not in str( eqn_lst )
    assert 'x2' not in str( eqn_lst )

    return eqn_lst


def usecase_B2():
    '''
    We compute the projective isomorphisms between
    the images of birational maps:
    
    f:P2--->X and g:P1xP1--->Y

    Further explanation of this example can be found 
    in the accompanying arxiv article on projective
    isomorphisms between rational surfaces.
    '''
    x = [ring( 'x' + str( i ) ) for i in range( 3 )]
    y = [ring( 'y' + str( i ) ) for i in range( 4 )]
    z = [ring( 'z' + str( i ) ) for i in range( 4 )]
    c = [ring( 'c' + str( i ) ) for i in range( 8 )]

    # the maps f and g are parametrizations of (linear projections of) toric surfaces
    f = ring( '[x0^6*x1^2,x0*x1^5*x2^2,x1^3*x2^5,x0^5*x2^3+x0^5*x2^3+x0^5*x1*x2^2]' )
    g = ring( '[y0^3*y1^2*y2^5,y1^5*y2^3*y3^2,y0^2*y1^3*y3^5,y0^5*y2^2*y3^3+y0^4*y1*y2^3*y3^2]' )
    g = [g[0], g[1] + g[0], g[2], g[3] + g[2]]
    SETools.p( 'f =', len( f ), f )
    SETools.p( 'g =', len( g ), g )
    assert sage_gcd( f ) == 1
    assert sage_gcd( g ) == 1

    # we compute the implicit equations of the images of the maps f and g
    eqf = sage_ideal( [z[i] - f[i] for i in range( 4 )] ).elimination_ideal( x ).gens()
    SETools.p( 'eqf =', eqf )
    assert len( eqf ) == 1
    assert eqf[0].degree() == 26
    eqg = sage_ideal( [z[i] - g[i] for i in range( 4 )] ).elimination_ideal( y ).gens()
    SETools.p( 'eqg =', eqg )
    assert len( eqg ) == 1
    assert eqg[0].degree() == 26

    # We compute the coefficient matrix Mf and its kernel Kf
    Mf = SERing.get_matrix_P2( f )
    Kf = Mf.right_kernel_matrix().T
    SETools.p( 'Mf  =', Mf.dimensions(), list( Mf ) )
    SETools.p( 'Kf  =', Kf.dimensions(), list( Kf ) )
    assert ( Mf * Kf ).is_zero()

    # we do a basepoint analysis for f and g
    bf = LinearSeries( SERing.conv( f ), PolyRing( 'x,y,z', True ) ).get_bp_tree()
    SETools.p( 'bf =', bf )
    bg = LinearSeries( SERing.conv( g ), PolyRing( 'x,y,v,w', True ) ).get_bp_tree()
    SETools.p( 'bg =', bg )

    # we compute maps to P1xP1 from two pencils
    PolyRing.reset_base_field()
    bpt = BasePointTree()
    bpt.add( 'y', ( 0, 0 ) , 1 )
    pen1 = SERing.conv( LinearSeries.get( [1], bpt ).pol_lst )
    SETools.p( 'pen1 =', pen1 )
    assert set( [x[0], x[1]] ) == set( pen1 )
    # thus the first pencil defines a map pen1: (x0:x1:x2) |--> (x0:x1)
    bpt = BasePointTree()
    bpt.add( 'x', ( 0, 0 ) , 1 )
    pen2 = SERing.conv( LinearSeries.get( [1], bpt ).pol_lst )
    SETools.p( 'pen2 =', pen2 )
    assert set( [x[0], x[2]] ) == set( pen2 )
    # thus the second pencil defines a map pen2: (x0:x1:x2) |--> (x0:x2)
    # We find that
    #     pen1 x pen2: P2-->P1xP1, (x0:x1:x2) |--> (x0:x1;x0:x2) and
    #     pen2 x pen1: P2-->P1xP1, (x0:x1:x2) |--> (x0:x2;x0:x1)

    # we find the following compatible reparametrizations
    # by composing the maps pen1 x pen2 and pen2 x pen1
    # with a parametrized map in the identity component of Aut(P1xP1).
    r0 = {y[0]:c[0] * x[0] + c[1] * x[1],
          y[1]:c[2] * x[0] + c[3] * x[1],
          y[2]:c[4] * x[0] + c[5] * x[2],
          y[3]:c[6] * x[0] + c[7] * x[2]}
    r1 = {y[0]:c[0] * x[0] + c[1] * x[2],
          y[1]:c[2] * x[0] + c[3] * x[2],
          y[2]:c[4] * x[0] + c[5] * x[1],
          y[3]:c[6] * x[0] + c[7] * x[1]}
    # Remark: all substitutions with .subs(...) are performed at the same time.

    ###################################################
    # reparametrization r0                            #
    ###################################################

    # compose g with reparametrization r0
    gcd0 = sage_gcd( [ comp.subs( r0 ) for comp in g ] )
    assert gcd0 == 1
    gr0 = [ comp.subs( r0 ) / gcd0 for comp in g ]
    SETools.p( 'gr0 =', len( gr0 ), gcd0, gr0 )
    assert SERing.get_degree( gr0 ) == 10
    assert SERing.get_degree( f ) == 8

    # find conditions on c so that gr0 has the same basepoints as f
    eqn0_lst = usecase_B2_helper_bp( gr0 )
    eqn0_lst += [ ring( '(c0*c3-c1*c2)*(c4*c7-c5*c6)*t-1' ) ]
    prime0_lst = sage_ideal( eqn0_lst ).elimination_ideal( ring( 't' ) ).primary_decomposition()
    SETools.p( 'eqn0_lst =', len( eqn0_lst ), eqn0_lst )
    for prime0 in prime0_lst:
        SETools.p( '\t', prime0.gens() )
    sol00 = {c[1]:0, c[2]:0, c[5]:0, c[6]:0, c[0]:1, c[4]:1}  # notice that wlog c0=c4=1
    sol01 = {c[0]:0, c[3]:0, c[4]:0, c[7]:0, c[1]:1, c[5]:1}  # notice that wlog c1=c5=1
    assert len( prime0_lst ) == 2
    assert set( [gen.subs( sol00 ) for gen in prime0_lst[0].gens()] ) == set( [0] )
    assert set( [gen.subs( sol01 ) for gen in prime0_lst[1].gens()] ) == set( [0] )

    # sol00: notice that c3!=0 and c7!=0
    gcd00 = sage_gcd( [ comp.subs( sol00 ) for comp in gr0] )
    assert gcd00 == x[0] * x[0]
    gr00 = [ comp.subs( sol00 ) / gcd00 for comp in gr0]
    SETools.p( 'gr00 =', len( gr00 ), gcd00, gr00 )
    assert SERing.get_degree( gr00 ) == 8
    Mgr00 = SERing.get_matrix_P2( gr00 )
    assert Mgr00.dimensions() == ( 4, 45 )
    # find conditions for c so that Mgr00 has the same kernel as the matrix of f
    p00_lst = sage_ideal( ( Mgr00 * Kf ).list() + [ring( 'c3*c7*t-1' )] ).elimination_ideal( ring( 't' ) ).primary_decomposition()
    assert [p00.gens() for p00 in p00_lst] == [[2 * c[3] - c[7]]]
    Mgr00 = Mgr00.subs( {c[7]:2 * c[3]} )
    SETools.p( 'Mgr00 =', Mgr00.dimensions(), list( Mgr00 ) )
    # found a solution: Mgr00

    # sol01: notice that c2!=0 and c6!=0
    gcd01 = sage_gcd( [ comp.subs( sol01 ) for comp in gr0] )
    assert gcd01 == x[0] * x[0]
    gr01 = [ comp.subs( sol01 ) / gcd01 for comp in gr0]
    SETools.p( 'gr01 =', len( gr01 ), gcd01, gr01 )
    assert SERing.get_degree( gr01 ) == 8
    assert [] == sage_ideal( ( SERing.get_matrix_P2( gr01 ) * Kf ).list() + [ring( 'c2*c6*t-1' )] ).elimination_ideal( ring( 't' ) ).primary_decomposition()
    # --> no solution


    ###################################################
    # reparametrization r1                            #
    ###################################################

    # compose g with reparametrization r1
    gcd1 = sage_gcd( [ comp.subs( r1 ) for comp in g ] )
    assert gcd1 == 1
    gr1 = [ comp.subs( r1 ) / gcd1 for comp in g ]
    SETools.p( 'gr1 =', gcd1, gr1 )
    assert SERing.get_degree( gr1 ) == 10
    assert SERing.get_degree( f ) == 8

    # find conditions on c so that gr1 has the same basepoints as f
    eqn1_lst = usecase_B2_helper_bp( gr1 )
    eqn1_lst += [ ring( '(c0*c3-c1*c2)*(c4*c7-c5*c6)*t-1' ) ]
    SETools.p( 'eqn1_lst =', len( eqn1_lst ), eqn1_lst )
    prime1_lst = sage_ideal( eqn1_lst ).elimination_ideal( ring( 't' ) ).primary_decomposition()
    for prime1 in prime1_lst:
        SETools.p( '\t', prime1.gens() )
    sol10 = {c[0]:0, c[3]:0, c[5]:0, c[6]:0, c[1]:1, c[4]:1}  # notice that wlog c1=c4=1
    sol11 = {c[1]:0, c[2]:0, c[4]:0, c[7]:0, c[0]:1, c[5]:1}  # notice that wlog c0=c5=1
    assert len( prime1_lst ) == 2
    assert set( [gen.subs( sol10 ) for gen in prime1_lst[0].gens()] ) == set( [0] )
    assert set( [gen.subs( sol11 ) for gen in prime1_lst[1].gens()] ) == set( [0] )

    # sol10: notice that c2!=0 and c7!=0
    gcd10 = sage_gcd( [ comp.subs( sol10 ) for comp in gr1] )
    assert gcd10 == x[0] * x[0]
    gr10 = [ comp.subs( sol10 ) / gcd10 for comp in gr1]
    SETools.p( 'gr10 =', len( gr10 ), gcd10, gr10 )
    assert SERing.get_degree( gr10 ) == 8
    assert [] == sage_ideal( ( SERing.get_matrix_P2( gr10 ) * Kf ).list() + [ring( 'c2*c7*t-1' )] ).elimination_ideal( ring( 't' ) ).primary_decomposition()
    # --> no solution

    # sol11: notice that c3!=0 and c6!=0
    gcd11 = sage_gcd( [ comp.subs( sol11 ) for comp in gr1] )
    assert gcd11 == x[0] * x[0]
    gr11 = [ comp.subs( sol11 ) / gcd11 for comp in gr1]
    SETools.p( 'gr11 =', len( gr11 ), gcd11, gr11 )
    assert SERing.get_degree( gr11 ) == 8
    assert [] == sage_ideal( ( SERing.get_matrix_P2( gr11 ) * Kf ).list() + [ring( 'c3*c6*t-1' )] ).elimination_ideal( ring( 't' ) ).primary_decomposition()
    # --> no solution

    ###################################################
    # compute extended matrices                       #
    ###################################################

    # Mgr00 is the only case we have to consider as other cases had no solution
    Mgr = Mgr00

    # compute the projective isomorphism between the images of f and g
    # in terms of parametrized matrix U
    Ef = sage_matrix( sage_QQ, list( Mf ) + list( Kf.T ) )
    Egr = sage_matrix( list( Mgr ) + list( Kf.T ) )
    UpI = Egr * ~Ef
    assert ( UpI.submatrix( 4, 4 ) - sage_identity_matrix( 41 ) ).is_zero()
    U = UpI.submatrix( 0, 0, 4, 4 )
    U = U / sage_gcd( U.list() )
    SETools.p( 'U =', U.dimensions(), list( U ), '\n' + str( U ) )

    # verify whether U*f is a parametrization for X for all (c0,...,c7)
    Uf = list( U * sage_vector( f ) )
    eqg_sub = [ eq.subs( {z[i]:Uf[i] for i in range( 4 )} ) for eq in eqg ]
    if eqg_sub != [0]:
        SETools.p( 'Uf  =', len( Uf ), Uf )
        SETools.p( 'eqg_sub=', len( eqg_sub ), eqg_sub )
    assert eqg_sub == [0]


def usecase_B4():
    '''
    We compute the projective automorphism of the 
    rational normal scrolls that is parametrized 
    the birational map f: P2 ---> X.
     
    Further explanation of this example can be found 
    in the accompanying arxiv article on projective
    isomorphisms between rational surfaces.
    '''

    # e0-e1
    p1 = ( 0, 0 );p2 = ( 1, 0 );p3 = ( 0, 1 )
    PolyRing.reset_base_field()
    bpt = BasePointTree()
    bpt.add( 'z', p1, 1 )
    f0p1 = SERing.conv( LinearSeries.get( [1], bpt ).pol_lst )
    SETools.p( 'f0p1 =', len( f0p1 ), f0p1 )

    # e0-e2-e3
    bpt = BasePointTree()
    bpt.add( 'z', p2, 1 )
    bpt.add( 'z', p3, 1 )
    f1m2 = SERing.conv( LinearSeries.get( [1], bpt ).pol_lst )
    SETools.p( 'f1m2 =', len( f1m2 ), f1m2 )

    # 2e0-e1-e2-e3
    bpt = BasePointTree()
    bpt.add( 'z', p1, 1 )
    bpt.add( 'z', p2, 1 )
    bpt.add( 'z', p3, 1 )
    f1m1 = SERing.conv( LinearSeries.get( [2], bpt ).pol_lst )
    SETools.p( 'f1m1 =', len( f1m1 ), f1m1 )

    # 3e0-2e1-e2-e3
    bpt = BasePointTree()
    bpt.add( 'z', p1, 2 )
    bpt.add( 'z', p2, 1 )
    bpt.add( 'z', p3, 1 )
    f1m0 = SERing.conv( LinearSeries.get( [3], bpt ).pol_lst )
    SETools.p( 'f1m0 =', len( f1m0 ), f1m0 )

    # set generators for the graded coordinate ring of f
    U = [ring( 'x1' ), ring( 'x2' ), ring( 'x1+x2-x0' ), ring( 'x1*x2' )]
    u = ring( 'u0,u1,u2,u3' )

    # obtain monomials of weight (1,0)
    w_lst = [( 0, 1 ), ( 0, 1 ), ( 1, -2 ), ( 1, -1 )]
    M1m0 = SERing.get_wmon_lst( u, w_lst, 1, 0 )
    SETools.p( 'M1m0 =', len( M1m0 ), M1m0 )

    # compose F with projective isomorphism P
    z = [ring( 'z' + str( i ) ) for i in range( 6 )]
    c = [ring( 'c' + str( i ) ) for i in range( 8 )]
    dctZ = { u[i]:z[i] for i in range( 4 )}
    dctP = { z[0]:c[0] * u[0] + c[1] * u[1],
             z[1]:c[2] * u[0] + c[3] * u[1],
             z[2]:c[4] * u[2],
             z[3]:c[5] * u[3] + c[6] * u[0] * u[2] + c[7] * u[1] * u[2]}
    PoF = [ comp.subs( dctZ ).subs( dctP ) for comp in M1m0]
    SETools.p( 'PoF wrt u =', len( PoF ), PoF )
    PoF = [ comp.subs( {u[i]:U[i] for i in range( 4 )} ) for comp in PoF]
    PoF = [ comp / sage_gcd( PoF ) for comp in PoF]
    SETools.p( 'PoF =', len( PoF ), PoF )

    # recover matrix for automorphism P
    F = [ comp.subs( {u[i]:U[i] for i in range( 4 )} ) for comp in M1m0 ]
    M = []
    for pol in [ comp.subs( dctZ ).subs( dctP ) for comp in M1m0]:
        row = []
        for mon in M1m0:
            row += [ pol.coefficient( mon ) ]
        M += [row]
    M = sage_matrix( M )
    SETools.p( 'M =', M.dimensions(), '\n' + str( M ) )
    MF = list( M * sage_vector( F ) )
    assert MF == PoF

    # compute the inverse Q of F
    t = ring( 't' )
    x = ring( 'x0,x1,x2' )
    id = [ F[i] * z[0] - z[i] * F[0] for i in range( 5 ) ] + [t * F[0] - 1]
    I1 = sage_ideal( id ).elimination_ideal( [t, x[2] ] ).gens()
    I2 = sage_ideal( id ).elimination_ideal( [t, x[1] ] ).gens()
    I1 = [ elt for elt in I1 if elt.degree( x[0] ) == 1 and elt.degree( x[1] ) == 1 ][0]
    I2 = [ elt for elt in I2 if elt.degree( x[0] ) == 1 and elt.degree( x[2] ) == 1 ][0]
    Q0 = I1.coefficient( x[1] )
    Q1 = -I1.coefficient( x[0] )
    Q2 = I2.coefficient( x[0] )
    Q = [Q0, Q1, Q2]
    SETools.p( 'Q =', Q )
    QoF = [comp.subs( {z[i]:F[i] for i in range( 5 )} ) for comp in Q]
    QoF = [comp / sage_gcd( QoF ) for comp in QoF]
    assert QoF == [x[0], x[1], x[2]]

    # compute the composition
    QoPoF = [ comp.subs( {z[i]:PoF[i] for i in range( 5 )} ) for comp in Q]
    QoPoF = [ comp / sage_gcd( QoPoF ) for comp in QoPoF ]
    SETools.p( 'QoPoF =', len( QoPoF ), QoPoF )

    # from the compatible reparametrizations QoPoF we compute the
    # projective automorphisms U of X.
    f = f1m0
    gr = [ comp.subs( {x[i]:QoPoF[i] for i in range( 3 )} ) for comp in f]
    gcd_gr = sage_gcd( gr )

    gr = [ comp / gcd_gr for comp in gr ]
    Mf = SERing.get_matrix_P2( f )
    Mgr = SERing.get_matrix_P2( gr )
    Kf = Mf.right_kernel_matrix().T
    SETools.p( 'f   =', len( f ), f )
    SETools.p( 'gr  =', len( gr ), gr )
    SETools.p( '\t gcd_gr  =', gcd_gr )
    SETools.p( 'Mf  =', Mf.dimensions(), list( Mf ) )
    SETools.p( 'Mgr =', Mgr.dimensions(), list( Mgr ) )
    SETools.p( 'Kf  =', Kf.dimensions(), list( Kf ) )

    assert ( Mf * Kf ).is_zero()
    assert ( Mgr * Kf ).is_zero()

    Ef = sage_matrix( sage_QQ, list( Mf ) + list( Kf.T ) )
    Egr = sage_matrix( list( Mgr ) + list( Kf.T ) )
    UpI = Egr * ~Ef
    assert ( UpI.submatrix( 5, 5 ) - sage_identity_matrix( 5 ) ).is_zero()
    U = UpI.submatrix( 0, 0, 5, 5 )
    SETools.p( 'UpI =', UpI.dimensions(), list( UpI ) )
    SETools.p( 'U   =', U.dimensions(), list( U ) )

    # verify whether U*f is a parametrization for X for all (c0,...,c7)
    Uf = list( U * sage_vector( f ) )
    SETools.p( 'Uf  =', len( Uf ), Uf )
    eqX = sage_ideal( [ z[i] - f[i] for i in range( 5 )] ).elimination_ideal( [x[0], x[1], x[2]] ).gens()
    eqXs = [ eq.subs( {z[i]:Uf[i] for i in range( 5 )} ) for eq in eqX ]
    SETools.p( 'eqX =', len( eqX ), eqX )
    SETools.p( 'eqXs=', len( eqXs ), eqXs )
    assert eqXs == [0, 0, 0]


def usecase_B5():
    '''
    We compute the projective isomorphism between two 
    conic bundles that are parametrized by the 
    birational maps 
    
    ff: P2 ---> X     and    gg: P1xP1 ---> Y
     
    Further explanation of this example can be found 
    in the accompanying arxiv article on projective
    isomorphisms between rational surfaces.
    '''

    # we construct linear series associated to ff in order to determine
    # the generators of the graded coordinate ring of conic bundle X

    # basepoints in chart x0!=0;
    p1 = ( 0, 0 );p2 = ( 0, 1 );p3 = ( 1, 0 )

    # 0f+p = e0-e1
    PolyRing.reset_base_field()
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

    # by inspection we recover the generators of graded ring of ff
    U = U0, U1, U2, U3, U4 = ring( 'x1' ), ring( 'x2' ), ring( 'x1+x2-x0' ), ring( 'x1*x2' ), ring( '(x1+x2-x0)^2' )

    # compute bidegree (2,d) in order to find a relation between the generators
    u = u0, u1, u2, u3, u4 = ring( 'u0,u1,u2,u3,u4' )
    SETools.p( 'Compare number of monomials of given bi-weight with dimension predicted by the Riemann-Roch formula...' )
    for d in reversed( [-i for i in range( 8 )] ):
        w_lst = [( 0, 1 ), ( 0, 1 ), ( 1, -3 ), ( 1, -2 ), ( 1, -2 )]
        SETools.p( '\tweight=', ( 2, d ), ',\t#monomials=', len( SERing.get_wmon_lst( u, w_lst, 2, d ) ), ',\tRR=', 29 + 5 * d )

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

    # we construct linear series associated to gg in order to determine
    # the generators of the graded coordinate ring of conic bundle Y

    # determine and set basepoint tree
    ls = LinearSeries( SERing.conv( gg ), PolyRing( 'x,y,v,w' ) )
    bp_tree = ls.get_bp_tree()
    SETools.p( 'bp_tree(gg) =', bp_tree )
    tree_211 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    tree_211.add( 'xw', ( 0, 0 ), 2 ).add( 't', ( 1, 0 ), 1 )
    tree_211.add( 'yv', ( 0, 1 ), 1 )

    # 1g+0q = 4l0+2l1-2e1-e2-e3
    g1m0 = SERing.conv( LinearSeries.get( [4, 2], tree_211 ).pol_lst )
    SETools.p( 'g1m0 =', len( g1m0 ), g1m0 )

    # 1g-3q = (l0+l1-e1-e2-e3) + (b-e1)
    g1m3 = SERing.conv( LinearSeries.get( [1, 2], tree_211 ).pol_lst )
    SETools.p( 'g1m3 =', len( g1m3 ), g1m3 )

    # 1g-2q = 2l0+2l1-2e1-e2-e3
    g1m2 = SERing.conv( LinearSeries.get( [2, 2], tree_211 ).pol_lst )
    SETools.p( 'g1m2 =', len( g1m2 ), g1m2 )

    # 1g-1q = 3l0+2l1-2e1-e2-e3
    g1m1 = SERing.conv( LinearSeries.get( [3, 2], tree_211 ).pol_lst )
    SETools.p( 'g1m1 =', len( g1m1 ), g1m1 )

    # 2g-4q = 4l0+4l1-4e1-2e2-2e3
    tree_422 = BasePointTree( ['xv', 'xw', 'yv', 'yw'] )
    tree_422.add( 'xw', ( 0, 0 ), 4 ).add( 't', ( 1, 0 ), 2 )
    tree_422.add( 'yv', ( 0, 1 ), 2 )
    g2m4 = SERing.conv( LinearSeries.get( [4, 4], tree_422 ).pol_lst )
    SETools.p( 'g2m4 =', len( g2m4 ), g2m4 )

    # by inspection we recover the generators of graded ring of gg
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
    SETools.p( 'Q =', Q )  # [-z9, -z8, -z8, -z6 - 2*z7 + z8 + z9]

    # check the inverse
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

    # create a list of equations for the ci
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

    # solve for ci and put the solutions in dictionary form
    prime_lst = sage_ideal( rel_lst ).elimination_ideal( t ).primary_decomposition()
    SETools.p( 'prime_lst =', len( prime_lst ) )
    for prime in prime_lst:
        SETools.p( '\t', prime.gens() )
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

    #  compose compatible reparametrizations with gg
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

    # get coefficient matrix of ff and its kernel
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

        # check if the answer is correct by substituting into the equations of Y
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


    R = sage_PolynomialRing( sage_QQ, 'x0,x1,x2,x3,z0,z1,z2,z3,z4,t', order = 'deglex' )
    x0, x1, x2, x3, z0, z1, z2, z3, z4, t = R.gens()

    d = x1 ** 2 + x2 ** 2 + x3 ** 2
    f0 = d + x0 ** 2
    f1 = 2 * x0 * x1
    f2 = 2 * x0 * x2
    f3 = 2 * x0 * x3
    f4 = d - x0 ** 2
    eq = -z0 ** 2 + z1 ** 2 + z2 ** 2 + z3 ** 2 + z4 ** 2

    g = [ f1 * z0 - z1 * f0,
          f2 * z0 - z2 * f0,
          f3 * z0 - z3 * f0,
          f4 * z0 - z4 * f0,
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

    usecase_B1()
    # usecase_B2()
    # usecase_B4()
    # usecase_B5()
    # usecase_invert_map()

    #########################################
    #                                       #
    # End of list of use case methods.      #
    #                                       #
    #########################################

    SETools.end_timer()
    SETools.p( 'The End' )









