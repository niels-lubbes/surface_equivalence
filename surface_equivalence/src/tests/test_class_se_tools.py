'''
Use of this source code is governed by a MIT-style license that can be found in the LICENSE file.
Created on Dec 31, 2018
@author: Niels Lubbes
'''
from surface_equivalence.class_se_tools import SETools

class TestClassSETools:


    def test__p( self ):

        SETools.filter( None )
        assert SETools.p( 'Hello world!' ) != None

        SETools.filter( ['another_class.py'] )
        assert SETools.p( 'No output since called from another class.' ) == None

        SETools.filter_unset()
        assert SETools.p( 'Filter is disabled so output this string.' ) != None

        SETools.filter_reset()
        assert SETools.p( 'Filter is enabled again so do not output.' ) == None

        SETools.filter( ['test_class_se_tools.py'] )
        assert SETools.p( 'Only output if called from this class' ) != None


    def test__tool_dct( self ):

        tool1 = SETools()
        tool2 = SETools()

        # watch out to not use the default file name
        # otherwise it might take long to load the data
        test_fname = 'test_tools'
        key = 'test__tool_dct'

        dct = tool1.get_tool_dct( fname = test_fname )
        dct[key] = True
        tool1.save_tool_dct( fname = test_fname )

        assert key in tool1.get_tool_dct( fname = test_fname )
        assert key in tool2.get_tool_dct( fname = test_fname )

        tool1.set_enable_tool_dct( False )
        assert key not in tool1.get_tool_dct( fname = test_fname )
        assert key not in tool2.get_tool_dct( fname = test_fname )

        tool1.set_enable_tool_dct( True )
        assert key in tool1.get_tool_dct( fname = test_fname )
        assert key in tool2.get_tool_dct( fname = test_fname )


