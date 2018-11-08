"""This module contains a variety of tests for the csp.py module"""
import unittest

from cone_search_plus import csp
import astropy.units as q

 
class TestStuff(unittest.TestCase):
    """Tests for SourceList class"""
    def setUp(self):
        # Trappist-1
        self.ra = 346.6223683553692
        self.dec = -05.0413976917903

        # Look for targets within 2 arcminutes
        self.radius = 2*q.arcmin

    def test_SourceList(self):
        """Test if Filter object is created properly"""
        # Perform the search
        self.sourcelist = csp.SourceList([self.ra, self.dec], self.radius)

        # Make sure the target is Trappist-1
        target = self.sourcelist.target['MAIN_ID'][0].decode('UTF-8')
        self.assertTrue(target == 'TRAPPIST-1')

        # Make sure we find 6 sources
        self.assertEqual(len(self.sourcelist.sources), 6)