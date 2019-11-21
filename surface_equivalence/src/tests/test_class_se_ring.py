'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Dec 31, 2018
@author: Niels Lubbes
'''

import os

from tests.class_test_tools import TestTools

from surface_equivalence.class_se_tools import SETools

from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing
from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_SR
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_matrix


class TestSERing( TestTools ):


    def test__random_ZZ( self ):
        out = SERing.random_ZZ( 10 )
        assert out >= -10 and out <= 10


    def test__random_QQ( self ):
        out = SERing.random_QQ()
        assert out in sage_QQ


    def test__random_elt( self ):
        assert SERing.random_elt( range( 50 ) ) in range( 50 )


    def test__random_matrix_QQ( self ):
        m, n = 4, 9
        mat = SERing.random_matrix_QQ( m, n )
        print( mat )
        print( mat.rows() )

        assert mat.nrows() == m
        assert mat.ncols() == n
        for row in mat:
            for col in row:
                assert col in sage_QQ


    def test__random_inv_matrix_QQ( self ):
        m = 4
        mat = SERing.random_inv_matrix_QQ( m )
        print( mat.rows() )
        assert mat.is_invertible()


    def test__get_mon_P1xP1__22( self ):
        out = SERing.get_mon_P1xP1( 2, 2 )
        chk = '[s^2*u^2, s^2*u*v, s^2*v^2, s*t*u^2, s*t*u*v, s*t*v^2, t^2*u^2, t^2*u*v, t^2*v^2]'
        print( out )
        s, t, u, v = ring( 's,t,u,v' )
        for p in out:
            assert p.degree( s ) + p.degree( t ) == 2
            assert p.degree( u ) + p.degree( v ) == 2
        assert str( out ) == chk


    def test__get_mon_P1xP1__23( self ):
        out = SERing.get_mon_P1xP1( 2, 3 )
        chk = '[s^2*u^3, s^2*u^2*v, s^2*u*v^2, s^2*v^3, s*t*u^3, s*t*u^2*v, s*t*u*v^2, s*t*v^3, t^2*u^3, t^2*u^2*v, t^2*u*v^2, t^2*v^3]'
        print( out )
        s, t, u, v = ring( 's,t,u,v' )
        for p in out:
            assert p.degree( s ) + p.degree( t ) == 2
            assert p.degree( u ) + p.degree( v ) == 3
        assert str( out ) == chk


    def test__get_mon_P2( self ):
        out = SERing.get_mon_P2( 2 )
        chk = '[x^2, x*y, x*z, y^2, y*z, z^2]'
        print( out )
        x, y, z = ring( 'x,y,z' )
        for p in out:
            assert p.degree() == 2
        assert str( out ) == chk


    def test__get_bidegree( self ):
        m, n = 3, 5
        mon_lst = SERing.get_mon_P1xP1( m, n )
        assert SERing.get_bidegree( mon_lst ) == ( 3, 5 )

        mat = SERing.random_matrix_QQ( 4, len( mon_lst ) )
        pol_lst = list( mat * sage_vector( mon_lst ) )
        assert SERing.get_bidegree( pol_lst ) == ( 3, 5 )


    def test__get_matrix_P1P1( self ):
        mat = '[(13, 1, -2, 1/2, 1, 3, 3, 0, -1), (2/5, -2/5, -2, -1, -2, -1, -5/49, 19, 0), (1, -1, 0, -3/2, 1/2, -1/22, 2, 5, -8/9), (7/23, -1/2, 0, -1, 1, 32, 1/2, 0, 1)]'
        mat = sage_matrix( sage_QQ, ring( mat ) )
        mon_lst = SERing.get_mon_P1xP1( 2, 2 )
        pol_lst = list( mat * sage_vector( mon_lst ) )
        out = SERing.get_matrix_P1xP1( pol_lst )
        print( mat )
        print( mon_lst )
        print( pol_lst )
        print( out )
        assert mat == out


    def test__get_matrix_P2( self ):
        mat = '[(8, -1, 0, 1/17, 2, 0, 1, 4/15, 2, 0), (1, 1/6, -1, 1, -5/7, -3/2, -5, -78, -1, 0), (-1, 0, -3/2, -9/2, -1, -1, 1, -1, 0, -35), (-11/4, 0, 0, 0, -1, 1, 0, 1/2, -5/4, -1/2), (1/7, -7/11, -2/13, 1, -3, 1/2, -1, 1/2, 0, 1/2), (4, -1, -1, 4, 0, 7/4, 3/5, 0, 1, -35), (-1/5, 0, -2, -55/2, -1/3, 2, -2, -2/11, -1, -1/2), (-1/6, -1/2, 2, 0, 2, 2, 0, -1, -1, 8), (-1, -2/3, 5/33, 2, 0, -1/48, -1, -1/15, 1/2, 1), (0, 3/2, 0, -2, 1/5, -4/7, 1/2, 1/4, 1, 0)]'
        mat = sage_matrix( sage_QQ, ring( mat ) )
        mon_lst = SERing.get_mon_P2( 3 )
        pol_lst = list( mat * sage_vector( mon_lst ) )
        out = SERing.get_matrix_P2( pol_lst )
        print( mat )
        print( mon_lst )
        print( pol_lst )
        print( out )
        assert mat == out


    def test__compose_aut_P1P1__1( self ):

        matL = '[(3,2),(7,5)]'
        matL = sage_matrix( sage_QQ, ring( matL ) )
        assert matL.is_invertible()

        matR = '[(1,2),(3,7)]'
        matR = sage_matrix( sage_QQ, ring( matR ) )
        assert matR.is_invertible()

        pol_lst = SERing.get_mon_P1xP1( 1, 1 )
        out = SERing.compose_aut_P1P1( pol_lst, matL, matR, False )

        print( pol_lst )
        print( out )

        check = '[3*s*u + 2*t*u + 6*s*v + 4*t*v, 9*s*u + 6*t*u + 21*s*v + 14*t*v, 7*s*u + 5*t*u + 14*s*v + 10*t*v, 21*s*u + 15*t*u + 49*s*v + 35*t*v]'
        assert out == ring( check )


    def test__compose_aut_P1P1__2( self ):

        matL = '[(3,2),(7,5)]'
        matL = sage_matrix( sage_QQ, ring( matL ) )
        assert matL.is_invertible()

        matR = '[(1,2),(3,7)]'
        matR = sage_matrix( sage_QQ, ring( matR ) )
        assert matR.is_invertible()

        pol_lst = SERing.get_mon_P1xP1( 1, 1 )
        out = SERing.compose_aut_P1P1( pol_lst, matL, matR, True )

        print( pol_lst )
        print( out )

        #        [3*s*u + 2*t*u + 6*s*v + 4*t*v, 9*s*u + 6*t*u + 21*s*v + 14*t*v, 7*s*u + 5*t*u + 14*s*v + 10*t*v, 21*s*u + 15*t*u + 49*s*v + 35*t*v]
        check = '[3*s*u + 6*t*u + 2*s*v + 4*t*v, 9*s*u + 21*t*u + 6*s*v + 14*t*v, 7*s*u + 14*t*u + 5*s*v + 10*t*v, 21*s*u + 49*t*u + 15*s*v + 35*t*v]'
        assert out == ring( check )


    def test__maple_is_installed( self ):

        assert self.maple_is_installed()


    def test__solve( self ):

        if not self.maple_is_installed():
            return

        # s, t, u, v = ring( 's,t,u,v' )
        # a, b, c, d = ring( 'a,b,c,d' )
        # e, f, g, h = ring( 'e,f,g,h' )

        pol_lst = [ ring( 'a^2-b^2' ), ring( 'b+1' ), ring( '(c^2+1)*(c-1)' )]
        var_lst = ring( '[a,b,c]' )
        sol_lst = SERing.solve( pol_lst, var_lst )

        a, b, c = ring( 'a,b,c' )
        for sol in sol_lst:
            assert sol[a] ** 2 - sol[b] ** 2 == 0
            assert sol[b] + 1 == 0
            assert ( sol[c] ** 2 + 1 ) * ( sol[c] - 1 ) == 0


if __name__ == '__main__':

    SETools.filter( None )
    SETools.start_timer()

#     TestSERing().test__get_mon_P1xP1__22()
#     TestSERing().test__get_mon_P1xP1__23()
#     TestSERing().test__get_mon_P2()
#     TestSERing().test__random_matrix_QQ()
#     TestSERing().test__get_bidegree()
#     TestSERing().test__get_matrix_P1xP1()
    TestSERing().test__get_matrix_P2()
#     TestSERing().test__compose_aut_P1P1__1()
#     TestSERing().test__compose_aut_P1P1__2()
#     TestSERing().test__random_ZZ()
#     TestSERing().test__random_elt()
#     TestSERing().test__random_inv_matrix_QQ()
#     TestSERing().test__solve()

    SETools.end_timer()
