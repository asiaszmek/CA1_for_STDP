import unittest
from CA1 import CA1_PC
import utility_functions as uf
from neuron import h


class TestAddingSpines(unittest.TestCase):
    add_ER = False
    where_ca = ["soma", "apical"]
    @classmethod
    def setUpClass(cls):
        cls.cell = CA1_PC(add_ER=cls.add_ER, where_ca=cls.where_ca,
                          where_spines=["rad_t2"],
                          spine_pos={"rad_t2": [0.5, 0.51, 0.52, 0.53]},
                          receptor_list=[])
        
    def test_add_in_specified_positions(self):
        spine_no = len(self.cell.cell_filter("head"))
        self.assertEqual(4, spine_no)
        

if __name__ == "__main__":
    unittest.main()

