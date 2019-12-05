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


    def test__conv__xyvw( self ):

        inp = '[x,y,v,w,x*y*v*w,x+y+v+w,a0*a1*a2*a3*x^2*y^2+v^2*w^2]'
        out = SERing.conv( inp )
        chk = '[y0, y1, y2, y3, y0*y1*y2*y3, y0 + y1 + y2 + y3, a0*a1*a2*a3*y0^2*y1^2 + y2^2*y3^2]'

        print( out )
        assert str( out ) == chk


    def test__conv__xyz( self ):

        inp = '[x,y,z,x*y*z,x+y+z,a0*a1*a2*a3*x^2*y^2+z^2]'
        out = SERing.conv( inp )
        chk = '[x1, x2, x0, x0*x1*x2, x0 + x1 + x2, a0*a1*a2*a3*x1^2*x2^2 + x0^2]'

        print( out )
        assert str( out ) == chk


    def test__conv__xi( self ):

        inp = '[x1, x2, x0, x0*x1*x2, x0 + x1 + x2, a0*a1*a2*a3*x1^2*x2^2 + x0^2]'
        out = SERing.conv( inp )
        chk = '[x, y, z, x*y*z, x + y + z, a0*a1*a2*a3*x^2*y^2 + z^2]'


        print( out )
        assert str( out ) == chk


    def test__conv__yi( self ):

        inp = '[y0, y1, y2, y3, y0*y1*y2*y3, y0 + y1 + y2 + y3, a0*a1*a2*a3*y0^2*y1^2 + y2^2*y3^2]'
        out = SERing.conv( inp )
        chk = '[x, y, v, w, x*y*v*w, x + y + v + w, a0*a1*a2*a3*x^2*y^2 + v^2*w^2]'

        print( out )
        assert str( out ) == chk


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
        out = SERing.get_mon_P1xP1( 2, 2, 'y0,y1,y2,y3' )
        chk = '[y0^2*y2^2, y0^2*y2*y3, y0^2*y3^2, y0*y1*y2^2, y0*y1*y2*y3, y0*y1*y3^2, y1^2*y2^2, y1^2*y2*y3, y1^2*y3^2]'
        print( out )
        y0, y1, y2, y3 = ring( 'y0,y1,y2,y3' )
        for p in out:
            assert p.degree( y0 ) + p.degree( y1 ) == 2
            assert p.degree( y2 ) + p.degree( y3 ) == 2
        assert str( out ) == chk


    def test__get_mon_P1xP1__23( self ):
        out = SERing.get_mon_P1xP1( 2, 3, 'y0,y1,y2,y3' )
        chk = '[y0^2*y2^3, y0^2*y2^2*y3, y0^2*y2*y3^2, y0^2*y3^3, y0*y1*y2^3, y0*y1*y2^2*y3, y0*y1*y2*y3^2, y0*y1*y3^3, y1^2*y2^3, y1^2*y2^2*y3, y1^2*y2*y3^2, y1^2*y3^3]'
        print( out )
        y0, y1, y2, y3 = ring( 'y0,y1,y2,y3' )
        for p in out:
            assert p.degree( y0 ) + p.degree( y1 ) == 2
            assert p.degree( y2 ) + p.degree( y3 ) == 3
        assert str( out ) == chk


    def test__get_mon_P2( self ):
        out = SERing.get_mon_P2( 2, 'x0,x1,x2' )
        chk = '[x0^2, x0*x1, x0*x2, x1^2, x1*x2, x2^2]'
        print( out )
        x0, x1, x2 = ring( 'x0,x1,x2' )
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


    def test__get_degree( self ):
        d = 4
        mon_lst = SERing.get_mon_P2( d )
        print( mon_lst )
        out = SERing.get_degree( mon_lst )
        print( out )
        assert out == d

        mat = SERing.random_matrix_QQ( 4, len( mon_lst ) )
        pol_lst = list( mat * sage_vector( mon_lst ) )
        print( pol_lst )
        out = SERing.get_degree( pol_lst )
        print( out )
        assert out == d


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


    def test__maple_is_installed( self ):

        assert self.maple_is_installed()



if __name__ == '__main__':

    SETools.filter( None )
    SETools.start_timer()

#     TestSERing().test__conv__xyvw()
#     TestSERing().test__conv__xyz()
#     TestSERing().test__conv__xi()
#     TestSERing().test__conv__yi()
#     TestSERing().test__get_mon_P1xP1__22()
#     TestSERing().test__get_mon_P1xP1__23()
#     TestSERing().test__get_mon_P2()
#     TestSERing().test__random_matrix_QQ()
#     TestSERing().test__get_bidegree()
    TestSERing().test__get_degree()
#     TestSERing().test__get_matrix_P1xP1()
#     TestSERing().test__get_matrix_P2()
#     TestSERing().test__random_ZZ()
#     TestSERing().test__random_elt()
#     TestSERing().test__random_inv_matrix_QQ()


    SETools.end_timer()
