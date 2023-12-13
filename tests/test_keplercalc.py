from unittest import TestCase
from keplercalc import Keplercalc

class TestKeplercalc(TestCase):
    def test_hello1(self):
        output = Keplercalc.hello('Spencer')
        self.assertEqual(output, 'Hello Spencer!')
