'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Dec 31, 2018
@author: Niels Lubbes
'''

from surface_equivalence.class_se_tools import SETools

from surface_equivalence.class_se_ring import ring
from surface_equivalence.class_se_ring import SERing

from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_matrix
from surface_equivalence.sage_interface import sage_vector


def iso_P1xP1( p1_lst, p2_lst, aff_dct = '{a:1, e:1}' ):
    '''
    Takes as input two birational morphisms 
    F: P1xP1--->X and G: P1xP1--->Y, 
    where X, Y are surfaces in Pn and P1xP1 denotes 
    the fibre product the projective line P1 with itself.
    Notice that the F and G are morphisms and thus basepoint free.  
        
    Parameters
    ----------
    p1_lst: list<sage_POLY>
        List of polynomials in QQ[s,t,u,v] representing the map F.
        F should not be of bidegree (1,1)        
         
    p2_lst: list<sage_POLY>
        List of polynomials in QQ[s,t,u,v] representing the map G.
          
    aff_dct: string
        A string of a dictionary that selects an affine chart 
        of an automorphism of P1xP1 that is defined by 
        two 2x2 matrices:         
           ( [ a  b ]   [ e  f ] )
           ( [ c  d ] , [ g  h ] )
        The default value is {a:1, e:1}, which means that we
        assume that the automorphism is of the form:
           ( [ 1  b ]   [ 1  f ] )
           ( [ c  d ] , [ g  h ] )        
          
    Returns
    -------
    list<sage_Matrix>        
        A matrix M which represent a projective isomorphism 
        M: Pn--->Pn so that M(X)=Y. 
        If X and Y are not projective isomorphic, then None
        is returned.  
        
    '''
    #
    # obtain bi-degrees and verify whether both are the same
    #
    d1, d2 = SERing.get_bidegree( p1_lst )
    if set( SERing.get_bidegree( p2_lst ) ) != set( [ d1, d2 ] ):
        SETools.p( 'bidegree(p2_lst) != bidegree(p1_lst) =', ( d1, d2 ) )
        return None

    #
    # Create dictionaries for precomposing with a paramelement of Aut(P1xP1).
    #
    s, t, u, v = ring( 's,t,u,v' )
    a, b, c, d = ring( 'a,b,c,d' )
    e, f, g, h = ring( 'e,f,g,h' )
    dct1 = {}  # automorphism is in identity component
    dct1.update( {s:a * s + b * t, t:c * s + d * t} )
    dct1.update( {u:e * u + f * v, v:g * u + h * v} )
    dct2 = {}  # automorphism flips two factors of P1xP1
    dct2.update( {s:e * u + f * v, t:g * u + h * v} )
    dct2.update( {u:a * s + b * t, v:c * s + d * t} )

    #
    # Try for both dictionaries
    #
    for dct in [dct1, dct2]:
        SETools.p( 'dct  =', dct )

        #
        # Precompose p2_lst with parametrized element of Aut(P1xP1).
        #
        p2s_lst = [ p2.subs( dct ) for p2 in p2_lst ]
        SETools.p( 'p1_lst  =', len( p1_lst ), p1_lst )
        SETools.p( 'p2_lst  =', len( p2_lst ), p2_lst )
        SETools.p( 'p2s_lst =', len( p2s_lst ), p2s_lst )

        #
        # Obtain matrices of maps.
        #
        mat1 = SERing.get_matrix_P1xP1( p1_lst )
        mat2 = SERing.get_matrix_P1xP1( p2_lst )
        mat2s = SERing.get_matrix_P1xP1( p2s_lst )
        SETools.p( mat1.dimensions(), ', mat1  =', mat1.rows() )
        SETools.p( mat2.dimensions(), ', mat2  =', mat2.rows() )
        SETools.p( mat2s.dimensions(), ', mat2s =', mat2s.rows() )

        #
        # Construct equations for parametrized element of Aut(P1xP1).
        #
        # Suppose that the linear normalizations XN, YN of X and Y live
        # in projective r-space Pr and that X and Y are projective isomorphic.
        # In this case X and Y are projections of XN and YN from the same
        # kernel "ker1". Thus we require that "mat2s * ker1.T" is the
        # zero matrix.
        #
        ker1 = mat1.right_kernel_matrix()
        assert ( mat1 * ker1.T ).is_zero()
        zmat = mat2s * ker1.T
        SETools.p( zmat.dimensions(), ', zmat  =', zmat.rows() )

        #
        # take affine chart the equations
        #
        aff_dct = ring( aff_dct )
        var_lst = [ var for var in [a, b, c, d, e, f, g, h] if var not in aff_dct.keys() ]
        pol_lst = [ p.subs( aff_dct ) for p in zmat.list() ]

        #
        # solve the equations
        #
        sol_lst = SERing.solve( pol_lst, var_lst )
        if sol_lst == []:
            # either the surfaces are not projectively equivalent
            # or the automorphism of P1xP1 flips the two factors
            continue

        #
        # update the parametrized matrix mat2s with the 1st solution
        #
        for sol in sol_lst:
            sol.update( aff_dct )
        mat3 = mat2s.subs( sol_lst[0] )
        SETools.p( mat3.dimensions(), ', mat3  =', mat3.rows() )

        #
        # extend the matrices mat1 and mat3 with their common kernel
        # and find a matrix mat4 so that mat4*mat1_ext == mat3_ext
        #
        mat1_ext = sage_matrix( sage_QQ, mat1.rows() + ker1.rows() )
        mat3_ext = sage_matrix( sage_QQ, mat3.rows() + ker1.rows() )
        mat4 = mat3_ext * mat1_ext.inverse()
        SETools.p( mat4.dimensions(), ', mat4  =', mat4.rows() )

        #
        # obtain the submatrix of mat4 that corresponds to the required
        # projective isomorphism
        #
        mat5 = mat4.submatrix( 0, 0, mat3.nrows(), mat3.nrows() )
        SETools.p( mat5.dimensions(), ', mat5  =', mat5.rows() )

        return mat5

    return None


