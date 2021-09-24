import unittest
from CA1 import CA1_PC
import utility_functions as uf
from neuron import h


class TestAddingSpecies(unittest.TestCase):
    add_ER = False
    where_ca = ["soma", "apical"]
    standard_list = ["Calmodulin", "Calbindin", "Fixed"]
    where_spines = []

    @classmethod
    def setUpClass(cls):
        buffer_list = cls.standard_list + ["OGB1", "BF2", "Fluo3",
                                           "Mg Green"]
        cls.cell = CA1_PC(add_ER=cls.add_ER, where_ca=cls.where_ca,
                      where_spines=cls.where_spines,
                      buffer_list=buffer_list, recompile=False)
        
    def test_add_in_OGB1_buffer(self):
        self.assertEqual(2, len(self.cell.buffers["OGB1"]))

    def test_add_in_BF2_buffer(self):
        self.assertEqual(2, len(self.cell.buffers["BF2"]))

    def test_add_in_Mg_Green_buffer(self):
        self.assertEqual(2, len(self.cell.buffers["Mg Green"]))

    def test_add_in_Fluo3_buffer(self):
        self.assertEqual(2, len(self.cell.buffers["Fluo3"]))

    def test_add_in_Calbindin_buffer(self):
        self.assertEqual(2, len(self.cell.buffers["Calbindin"]))

    def test_add_in_Calmodulin_buffer(self):
        self.assertEqual(3, len(self.cell.buffers["Calmodulin"]))

    def test_add_in_fixed_buffer(self):
        self.assertEqual(2, len(self.cell.buffers["Fixed"]))
        
    def test_add_all_buffer(self):
        self.assertEqual(7, len(self.cell.buffers.keys()))


if __name__ == "__main__":
   
    
    unittest.main()

