'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Jan 1, 2019
@author: Niels Lubbes
'''

import os

from surface_equivalence.sage_interface import sage_maple


class TestTools( object ):

    def __clean__( self, str0, s_lst ):
        for s in s_lst:
            while s in str0:
                str0 = str0.replace( s, '' )
        return str0


    def equal_output_strings( self, str1, str2, s_lst = ' \n\t' ):
        '''
        This method can be used to test outputs of methods, 
        where the output is considered as a String and where 
        characters such as spaces, newlines or tabs do not 
        matter for testing the correctness. 
                
        Parameters
        ----------
        str1: string
        str2: string
        s_lst: list<string>    
            A list of strings. A String is considered as a list 
            of strings of length one (ie. characters).
        Returns
        -------
        boolean
            Returns True if "str1" and "str2" are equal after we replace 
            all occurrences of strings in "s_lst" with the empty string ''.  
            This method returns False otherwise.
        '''

        return self.__clean__( str1, s_lst ) == self.__clean__( str2, s_lst )


    def maple_is_installed( self ):
        '''
        Returns
        -------
        boolean
            True if maple is installed and false otherwise.
        '''

        os.environ['PATH'] += os.pathsep + '/home/niels/Desktop/n/app/maple/link/bin'

        try:
            sage_maple.eval( '1 + 1' )
            return True
        except:
            print( 'WARNING: Maple is not installed.' )
            return False





if __name__ == '__main__':

    # print( TestTools().equal_output_strings( 'a', 'b' ) )
    pass

