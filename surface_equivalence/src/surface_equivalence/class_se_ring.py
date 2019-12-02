'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Dec 31, 2018
@author: Niels Lubbes
'''

from surface_equivalence.class_se_tools import SETools

from surface_equivalence.sage_interface import sage_PolynomialRing
from surface_equivalence.sage_interface import sage_QQ
from surface_equivalence.sage_interface import sage_ZZ
from surface_equivalence.sage_interface import sage__eval
from surface_equivalence.sage_interface import sage_solve
from surface_equivalence.sage_interface import sage_Compositions
from surface_equivalence.sage_interface import sage_matrix
from surface_equivalence.sage_interface import sage_vector
from surface_equivalence.sage_interface import sage_solve
from surface_equivalence.sage_interface import sage_SR
from surface_equivalence.sage_interface import sage_maple

def ring( s ):
    return SERing.ring( s )


class SERing:
    '''
    This class represents a polynomial ring 
    
    Attributes
    ----------
    R : sage_PolynomialRing
        Polynomial ring QQ[a0,...,a4,c0,...,c19,u0,...,u4,x0,...,x2,y0,...,y3,z0,...,z19].
    '''
    x_lst = ['x' + str( i ) for i in range( 3 )]
    y_lst = ['y' + str( i ) for i in range( 4 )]
    z_lst = ['z' + str( i ) for i in range( 20 )]
    a_lst = ['a' + str( i ) for i in range( 5 )]
    c_lst = ['c' + str( i ) for i in range( 20 )]
    u_lst = ['u' + str( i ) for i in range( 5 )]

    R = sage_PolynomialRing( sage_QQ, a_lst + c_lst + u_lst + x_lst + y_lst + z_lst )

    @staticmethod
    def ring( expr ):
        return sage__eval( str( expr ), SERing.R.gens_dict() )


    @staticmethod
    def conv( pol_lst ):
        '''
        Renames variables as follows:
        x0,x1,x2     <---> z,x,y
        or
        y0,y1,y2,y3  <---> x,y,v,w        
        '''
        pol_lst = str( pol_lst )

        NR = sage_PolynomialRing( sage_QQ, 'a0,a1,a2,a3,a4,x0,x1,x2,y0,y1,y2,y3,x,y,z,v,w' )
        a0, a1, a2, a3, a4, x0, x1, x2, y0, y1, y2, y3, x, y, z, v, w = NR.gens()

        inring = True
        if 'x0' in pol_lst or 'x1' in pol_lst or 'x2' in pol_lst:
            dct = {x0:z, x1:x, x2:y}
            inring = False
        elif 'y0' in pol_lst or 'y1' in pol_lst or 'y2' in pol_lst or 'y3' in pol_lst:
            dct = {y0:x, y1:y, y2:v, y3:w}
            inring = False
        elif 'v' in pol_lst or 'w' in pol_lst:
            dct = {x:y0, y:y1, v:y2, w:y3}
        else:
            dct = {x:x1, y:x2, z:x0}

        pol_lst = sage__eval( pol_lst, NR.gens_dict() )
        pol_lst = [ pol.subs( dct ) for pol in pol_lst ]

        if inring:
            return SERing.ring( pol_lst )
        else:
            return pol_lst


    @staticmethod
    def get_mon_P1xP1( m, n, vars = 'y0,y1,y2,y3' ):
        '''        
        Parameters
        ----------
        m: int
        
        n: int
        
        vars: string
            Names of the variables.            
    
        Returns
        -------
        list<sage_POLY>
            All monomials of the form y0^%*y1^%*y2^%*y3^% of bi-degree (m,n) 
            in (y0,y1) and (y2,y3), respectively.
        '''

        y0, y1, y2, y3 = ring( vars )
        mon_lst = []
        for a, b, c, d in sage_Compositions( m + n + 4, length = 4 ):
                if a + b == m + 2 and c + d == n + 2:
                    mon_lst += [ y0 ** ( a - 1 ) * y1 ** ( b - 1 ) * y2 ** ( c - 1 ) * y3 ** ( d - 1 ) ]
        return mon_lst


    @staticmethod
    def get_mon_P2( d, vars = 'x0,x1,x2' ):
        '''        
        Parameters
        ----------
        d: int
        
        vars: string
            Names of the variables.
            
        Returns
        -------
        list<sage_POLY>
            All monomials of the form x0^%*x1^%*x2^% of degree d.         
        '''
        x0, x1, x2 = ring( vars )
        mon_lst = []
        for a, b, c in sage_Compositions( d + 3, length = 3 ):
            mon_lst += [ x0 ** ( a - 1 ) * x1 ** ( b - 1 ) * x2 ** ( c - 1 ) ]
        return mon_lst


    @staticmethod
    def get_bidegree( pol_lst, vars = 'y0,y1,y2,y3' ):
        '''
        Parameters
        ----------        
        pol_lst : list<SERing.R> 
            List of homogeneous polynomials of equal bi-degree (m,n) in 
            (y0:y1;y2:y3) where the variables are specified by vars.
        
        vars: string
            Names of the variables.
                        
        Returns
        -------
        (int,int)
            A pair of integers (m,n) defining the bi-degree of input polynomials.
        '''
        pol_lst = ring( pol_lst )
        y0, y1, y2, y3 = ring( 'y0,y1,y2,y3' )
        m = max( [ p.degree( y0 ) for p in pol_lst] )
        n = max( [ p.degree( y2 ) for p in pol_lst] )
        return m, n


    @staticmethod
    def get_matrix_P1xP1( pol_lst, vars = 'y0,y1,y2,y3' ):
        '''
        Obtains the matrix M so that wrt monomial basis vector v defined 
        by SERing.get_mon_P1xP1() we recover pol_lst as M*v.
        
                
        Parameters
        ----------        
        pol_lst : list<SERing.R> 
            List of homogeneous polynomials of equal bi-degree in 
            (y0:y1;y2:y3) where the variables are specified by vars. 
            
        vars: string
            Names of the variables.           
            
        Returns
        -------  
        sage_Matrix<SERing.R> 
            A matrix with polynomials in QQ[a,b,c,d,e,f,g,h].
        '''
        pol_lst = ring( pol_lst )
        m, n = SERing.get_bidegree( pol_lst, vars )
        mon_lst = SERing.get_mon_P1xP1( m, n, vars )

        SETools.p( 'm = ' + str( m ) + ', n = ' + str( n ) )

        mat = []
        for pol in pol_lst:
            row = []
            for mon in mon_lst:
                row += [ pol.coefficient( mon ) ]
            mat += [row]

        return sage_matrix( mat )


    @staticmethod
    def get_matrix_P2( pol_lst, vars = 'x0,x1,x2' ):
        '''
        Obtains the matrix M so that wrt monomial basis vector v defined 
        by SERing.get_mon_P1xP1() we recover pol_lst as M*v.
        
                
        Parameters
        ----------        
        pol_lst : list<SERing.R> 
            List of homogeneous polynomials of equal degree in (x0:x1:x2), 
            where the variables are specified by vars. 
            
        vars: string
            Names of the variables.    
            
        Returns
        -------  
        sage_Matrix<SERing.R> 
            The coefficient matrix of pol_lst whose entries are polynomials.
        '''
        pol_lst = ring( pol_lst )
        d = pol_lst[0].degree()
        mon_lst = SERing.get_mon_P2( d, vars )

        SETools.p( 'd = ' + str( d ) )

        mat = []
        for pol in pol_lst:
            row = []
            for mon in mon_lst:
                row += [ pol.coefficient( mon ) ]
            mat += [row]

        return sage_matrix( mat )


    @staticmethod
    def compose_aut_P1P1( pol_lst, matL, matR, flip = False ):
        '''
        Precomposes map P1xP1--->X with an automorphism P1xP1.
        We represent Aut(P1xP1) with two 2x2 matrices matL and matR, that
        act on the left and right factor of P1xP1 respectively. 
        If flip is True,  then the factors are flipped as well.

        Notes
        -----
        Internally the variables a,b,c,d are used as temporary place holders.
                        
        Parameters
        ----------        
        pol_lst : list<SERing.R> 
            List of homogeneous polynomials of equal bi-degree in QQ[s,t,u,v].
            This list defines the map P1xP1--->X.
        
        matL : sage_matrix<QQ>
            A 2x2 matrix representing an automorphism of projective line P1.
        
        matR : sage_matrix<QQ>
            A 2x2 matrix representing an automorphism of projective line P1.

        flip : boolean
            Flips the left and right factor of P1xP1.
            
        Returns
        -------  
        list<SERing.R>
            A list of polynomials in s,t,u,v which are defined by pol_lst
            composed with an automorphism of P1P1
        '''
        s, t, u, v = ring( 's,t,u,v' )
        a, b, c, d = ring( 'a,b,c,d' )

        if not flip:
            sn, tn = matL * sage_vector( [a, b] )
            un, vn = matR * sage_vector( [c, d] )
        else:
            sn, tn = matL * sage_vector( [c, d] )
            un, vn = matR * sage_vector( [a, b] )

        dct1 = {s:sn, t:tn, u:un, v:vn}
        dct2 = {a:s, b:t, c:u, d:v}

        out_lst = [pol.subs( dct1 ).subs( dct2 ) for pol in pol_lst ]
        SETools.p( 'out_lst =', out_lst )

        return out_lst


    @staticmethod
    def random_ZZ( val ):
        '''
        Parameters
        ----------
        val : int  
            An integer.
        
        Returns
        -------
        int
            A random element in the interval [-val,val]
        '''
        return int( sage_ZZ.random_element( -val, val + 1 ) )


    @staticmethod
    def random_QQ():
        '''        
        Returns
        -------
        sage_QQ
            A random rational number
        '''
        return sage_QQ.random_element()


    @staticmethod
    def random_matrix_QQ( m, n ):
        '''
        Parameters
        ----------
        m : int
        n : int
        
        Returns
        -------
        sage_matrix<sage_QQ>        
            A random mxn-matrix with rational number as entries. 
        '''
        mat = []
        for ri in range( m ):
            row = []
            for ci in range( n ):
                row += [ SERing.random_QQ() ]
            mat += [row]
        mat = sage_matrix( sage_QQ, mat )
        return mat


    @staticmethod
    def random_inv_matrix_QQ( m ):
        '''
        Parameters
        ----------
        m : int
        
        Returns
        -------
        sage_matrix<sage_QQ>        
            A random invertible mxm-matrix with rational number as entries. 
        '''
        while True:
            mat = SERing.random_matrix_QQ( m, m )
            if mat.is_invertible():
                return mat


    @staticmethod
    def random_elt( lst ):
        '''
        Parameters
        ----------
        lst : list
            A list.
        
        Returns
        -------
        object        
            A random element in "lst".
        '''
        idx = int( sage_ZZ.random_element( 0, len( lst ) ) )
        return lst[idx]


    @staticmethod
    def solve( pol_lst, var_lst, sol_QQ = True ):
        '''
        Parameters
        ----------
        pol_lst : list<SERing.R> 
            List of polynomials.
    
        var_lst : list<SERing.R>
            List of variables corresponding to genetors "SERing.R". 
        
        sol_QQ : boolean
            If True, then only return solutions over QQ.
        
        Returns
        -------
        list<dict>
            A list of solutions for the zeroset of the ideal generated by the polynomials in "pol_lst".  
            The solutions may be parametrized in terms of variables r%, where % is an integer.
            Each solution is represented by a dictionary where 
            * keys are elements of "var_lst" (coerced to sage_SR if sol_QQ==False)             
            * values are solutions (over QQ if sol_QQ==True and over sage_SR otherwise)
            If sol_QQ==True, then everything is coerced to SERing.R.
        '''
        # pol_lst = ring( pol_lst )
        # var_lst = ring( var_lst )
        SETools.p( 'pol_lst =', pol_lst )
        SETools.p( 'var_lst =', var_lst )

        #
        # compute groebner basis with Maple
        #
        sage_maple.eval( 'with(Groebner);' )
        sage_maple.eval( 'gb := Basis( ' + str( pol_lst ) + ', plex(' + str( var_lst )[1:-1] + ') );' )
        pol_lst = ring( sage_maple.eval( 'lprint(gb);' ) )
        SETools.p( 'pol_lst (Groebner basis) =', pol_lst )

        #
        # cast to Sage symbolic ring SR and call solve.
        #
        spol_lst = [ sage_SR( str( pol ) ) for pol in pol_lst ]
        svar_lst = [ sage_SR( str( var ) ) for var in var_lst ]
        sol_lst = sage_solve( spol_lst, svar_lst, solution_dict = True )
        SETools.p( 'sol_lst (CC) =', sol_QQ, sol_lst )

        if not sol_QQ:
            return sol_lst

        out_lst = []
        for sol in sol_lst:
            sol_in_QQ = True
            for key in sol:
                sol_in_QQ &= sol[key] in sage_QQ
            if sol_in_QQ:
                out_lst += [sol]
        sol_lst = ring( out_lst )
        SETools.p( 'sol_lst (QQ) =', sol_QQ, sol_lst )

        return ring( out_lst )



