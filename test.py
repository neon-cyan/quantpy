import unittest
from quatics_lexers import QuanticsParsers
from logall import ParseLogAll
from logall import I_ImportLogalls

class Test_quantpy(unittest.TestCase):
    def test_allene_inp_parse(self):
        with open('testfiles/allene_nopulse.inp') as f:
            data = f.read()
        all_in = QuanticsParsers.parse_input(data)

        self.assertEqual(all_in['name'], 'allene_nopulse')
        self.assertEqual(all_in['ngwp'], 16)
        self.assertEqual(all_in['data'], 'dd_data_nm')
        self.assertEqual(all_in['freqf'], 'allene_simple_32_freg.log')
        self.assertEqual(all_in['tfinal'], 40)
        self.assertEqual(all_in['tpsi'], 0.1)

    def test_allene_out_parse(self):
        with open('testfiles/output_allene') as f:
            data = f.read()
        all_out = QuanticsParsers.parse_output(data)

        self.assertEqual(all_out[0]['time'], 0.0)
        self.assertEqual(all_out[-1]['time'], 40.0)

        self.assertEqual(all_out[0]['DiagD'][0], 5395.1270)
        self.assertEqual(all_out[0]['DiagD'][-1], 48.1523)

        self.assertEqual(all_out[0]['GGP'][0], 15.7242)
        self.assertEqual(all_out[0]['GGP'][-1], -0.3816)

    def test_gly_logall_single_parse(self):
        with open('testfiles/out_testdir/glycine/dd_data_nm/gwp1_V1/gwp1_V1_dd_data_nm.logall') as f:
            data = f.read()
        all_out, _ = ParseLogAll.parse(data)
        with open('testfiles/out_testdir/glycine/dd_data_nm/gwp2_V1/gwp2_V1_dd_data_nm.logall') as f:
            data = f.read()
        all_out, _ = ParseLogAll.parse(data)

    def test_allene_logall_single_parse(self):
        with open('testfiles/out_testdir/allene/dd_data_nm/gwp1_V1/gwp1_V1_dd_data_nm.logall') as f:
            data = f.read()
        all_out, _ = ParseLogAll.parse(data)
        with open('testfiles/out_testdir/allene/dd_data_nm/gwp2_V1/gwp2_V1_dd_data_nm.logall') as f:
            data = f.read()
        all_out, _ = ParseLogAll.parse(data)

    def test_pullall_gly_partial(self):
        data = I_ImportLogalls('testfiles/out_testdir/glycine/dd_data_nm', 15, step_lim=10)

    def test_pullall_gly_full(self):
        data = I_ImportLogalls('testfiles/out_testdir/glycine/dd_data_nm', 15)

if __name__ == '__main__':
    unittest.main()