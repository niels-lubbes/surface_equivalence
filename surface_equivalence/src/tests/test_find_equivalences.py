'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Dec 31, 2018
@author: Niels Lubbes
'''

from surface_equivalence.find_equivalences import iso_P1xP1

from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing

from surface_equivalence.class_se_tools import SETools

from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_matrix

from tests.class_test_tools import TestTools


class TestFindEquivalences( TestTools ):

    def test__iso_P1xP1__1( self ):

        if not self.maple_is_installed():
            return


        mat1 = '[(13, 1, -2, 1/2, 1, 3, 3, 0, -1), (2/5, -2/5, -2, -1, -2, -1, -5/49, 19, 0), (1, -1, 0, -3/2, 1/2, -1/22, 2, 5, -8/9), (7/23, -1/2, 0, -1, 1, 32, 1/2, 0, 1)]'
        mat1 = sage_matrix( sage_QQ, ring( mat1 ) )

        mati = '[(-12/17, 2, 1/18, 5/7), (-12/7, 0, -1/235, 1/2), (1/7, 0, 1/13, 1/16), (-1, 6, 19, 9)]'
        mati = sage_matrix( sage_QQ, ring( mati ) )

        assert mati.is_invertible()

        p1_lst = list( mat1 * sage_vector( SERing.get_mon_P1xP1( 2, 2 ) ) )
        p2_lst = list( mati * sage_vector( p1_lst ) )

        print( mat1 )
        print( mati )

        assert list( mati * mat1 * sage_vector( SERing.get_mon_P1xP1( 2, 2 ) ) ) == p2_lst

        mat_out = iso_P1xP1( p1_lst, p2_lst )

        assert mat_out == mati


    def test__iso_P1xP1__2( self ):

        if not self.maple_is_installed():
            return


        mat1 = '[(13, 1, -2, 1/2, 1, 3, 3, 0, -1), (2/5, -2/5, -2, -1, -2, -1, -5/49, 19, 0), (1, -1, 0, -3/2, 1/2, -1/22, 2, 5, -8/9), (7/23, -1/2, 0, -1, 1, 32, 1/2, 0, 1)]'
        mat1 = sage_matrix( sage_QQ, ring( mat1 ) )

        mati = '[(-12/17, 2, 1/18, 5/7), (-12/7, 0, -1/235, 1/2), (1/7, 0, 1/13, 1/16), (-1, 6, 19, 9)]'
        mati = sage_matrix( sage_QQ, ring( mati ) )
        assert mati.is_invertible()

        matL = '[(3,2),(7,5)]'
        matL = sage_matrix( sage_QQ, ring( matL ) )
        assert matL.is_invertible()

        matR = '[(1,2),(3,7)]'
        matR = sage_matrix( sage_QQ, ring( matR ) )
        assert matR.is_invertible()

        p1_lst = list( mat1 * sage_vector( SERing.get_mon_P1xP1( 2, 2 ) ) )

        p2_lst = SERing.compose_aut_P1P1( p1_lst, matL, matR, False )
        p2_lst = list( mati * sage_vector( p2_lst ) )

        print( mat1 )
        print( '---' )
        print( mati )
        print( '---' )
        print( matL )
        print( '---' )
        print( matR )
        print( '---' )

        mat_out = iso_P1xP1( p1_lst, p2_lst )

        assert ( 20825 / 17 ) * mat_out == mati


    def test__iso_P1xP1__3( self ):

        if not self.maple_is_installed():
            return

        mat1 = '[(13, 1, -2, 1/2, 1, 3, 3, 0, -1), (2/5, -2/5, -2, -1, -2, -1, -5/49, 19, 0), (1, -1, 0, -3/2, 1/2, -1/22, 2, 5, -8/9), (7/23, -1/2, 0, -1, 1, 32, 1/2, 0, 1)]'
        mat1 = sage_matrix( sage_QQ, ring( mat1 ) )

        mati = '[(-12/17, 2, 1/18, 5/7), (-12/7, 0, -1/235, 1/2), (1/7, 0, 1/13, 1/16), (-1, 6, 19, 9)]'
        mati = sage_matrix( sage_QQ, ring( mati ) )
        assert mati.is_invertible()

        matL = '[(3,2),(7,5)]'
        matL = sage_matrix( sage_QQ, ring( matL ) )
        assert matL.is_invertible()

        matR = '[(1,2),(3,7)]'
        matR = sage_matrix( sage_QQ, ring( matR ) )
        assert matR.is_invertible()

        p1_lst = list( mat1 * sage_vector( SERing.get_mon_P1xP1( 2, 2 ) ) )

        p2_lst = SERing.compose_aut_P1P1( p1_lst, matL, matR, True )
        p2_lst = list( mati * sage_vector( p2_lst ) )

        print( mat1 )
        print( '---' )
        print( mati )
        print( '---' )
        print( matL )
        print( '---' )
        print( matR )
        print( '---' )

        mat_out = iso_P1xP1( p1_lst, p2_lst )

        assert ( 20825 / 17 ) * mat_out == mati


    def test__iso_P1xP1__random( self ):

        if not self.maple_is_installed():
            return

        d_lst = [2, 3]
        d1, d2 = SERing.random_elt( d_lst ), SERing.random_elt( d_lst )

        mon_lst = SERing.get_mon_P1xP1( d1, d2 )
        n = len( mon_lst )
        m = 4  # SERing.random_elt( range( 4, n + 1 ) )

        flip = SERing.random_elt( [True, False] )

        mat1 = SERing.random_matrix_QQ( m, n )
        while True:
            mati = SERing.random_inv_matrix_QQ( m )
            matL = SERing.random_inv_matrix_QQ( 2 )
            matR = SERing.random_inv_matrix_QQ( 2 )
            if mati[0, 0] != 0 and matL[0, 0] != 0 and matR[0, 0] != 0:
                break

        print( ' --- random parameters ---' )
        print( 'd1, d2 =', d1, d2 )
        print( 'm, n   =', m, n )
        print( 'flip   =', flip )
        print( '--- mat1 ---' )
        print( mat1 )
        print( '--- mati ---' )
        print( mati )
        print( '--- matL ---' )
        print( matL )
        print( '--- matR ---' )
        print( matR )
        print( '---' )


        p1_lst = list( mat1 * sage_vector( mon_lst ) )
        p2_lst = SERing.compose_aut_P1P1( p1_lst, matL, matR, flip )
        p2_lst = list( mati * sage_vector( p2_lst ) )

        mato = iso_P1xP1( p1_lst, p2_lst )

        print( mato )

        a00 = mati[0, 0]
        b00 = mato[0, 0]

        assert ( a00 / b00 ) * mato == mati



if __name__ == '__main__':

    # SETools.filter( ['find_equivalences.py'] )
    SETools.filter( None )
    SETools.start_timer()

    # TestFindEquivalences().test__iso_P1xP1__1()
    # TestFindEquivalences().test__iso_P1xP1__2()
    # TestFindEquivalences().test__iso_P1xP1__3()
    TestFindEquivalences().test__iso_P1xP1__random()

    SETools.end_timer()
