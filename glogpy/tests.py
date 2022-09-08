import unittest
from glogpy.freqency_job import frequency_job as fj
from glogpy.dynamics_job import dynamics_job as dj
from glogpy.l118_traj_job import l118_job
from glogpy.job import gaussian_job
from glogpy.jobio import gaussian_jobio

class TestJobs(unittest.TestCase):
    # PART I : Test on a few different generic odd jobs

    def test_aim_irc(self):
        with open('glogpy/tests/ircaim.log', 'r') as f:
            data = f.read()
        gaussian_job(data)

    def test_aim_irc_io(self):
        gaussian_jobio('glogpy/tests/ircaim.log')

    def test_isop_dft(self):
        with open('glogpy/tests/isop.log', 'r') as f:
            data = f.read()
        gaussian_job(data)

    def test_isop_dft_io(self):
        gaussian_jobio('glogpy/tests/isop.log')

    # Test a partial IO job

    def test_benz_partial_io(self):
        gaussian_jobio('glogpy/tests/benzene_freq_partial.log',allow_partial=True)

    # PART II : Test frequency parser

    def test_freq_job_o2_nosymm(self):
        with open('glogpy/tests/o2_freq.log', 'r') as f:
            data = f.read()
        gj = fj(data)
        ans = gj.parse()
        self.assertEqual(ans['vibsyms'][0],'A')
        self.assertDictEqual(ans['atomnos'],{1: 8, 2: 8})
        self.assertEqual(ans['vibfreqs'][0],1643.2361)

    def test_freq_job_gly_nosymm(self):
        with open('glogpy/tests/gly_freq_nosymm.log', 'r') as f:
            data = f.read()
        gj = fj(data)
        ans = gj.parse()
        self.assertEqual(ans['vibsyms'][0],'A')
        self.assertEqual(ans['vibsyms'][-1],'A')

        self.assertDictEqual(ans['atomnos'],{1: 6, 2: 8, 3: 8, 4: 6, 5: 7, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1})

        self.assertEqual(ans['vibfreqs'][0],63.5772)
        self.assertEqual(ans['vibfreqs'][-1],3678.9799)

    def test_freq_job_benzene_symm(self):
        with open('glogpy/tests/benzene_freq.log', 'r') as f:
            data = f.read()
        gj = fj(data)
        ans = gj.parse()

        self.assertEqual(ans['vibsyms'][0],'E2U')
        self.assertEqual(ans['vibsyms'][-1],'A1G')

        self.assertDictEqual(ans['atomnos'],{1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 6, 8: 6, 9: 6, 10: 6, 11: 6, 12: 6})

        self.assertEqual(ans['vibfreqs'][0],415.5995)
        self.assertEqual(ans['vibfreqs'][-1],3231.3848)

    # TEST III : Test a few QuEh dynamics jobs

    def test_gly_dynamics(self):
        with open('glogpy/tests/qdyn_gly.log', 'r') as f:
            data = f.read()
        gj = dj(data)
        ans = gj.parse()

        self.assertDictEqual(ans['atomnos'],{1: 6, 2: 8, 3: 8, 4: 6, 5: 7, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1})
        #(ans)

    def test_allene_dynamics(self):
        with open('glogpy/tests/qdyn_allene.log', 'r') as f:
            data = f.read()
        gj = dj(data)
        ans = gj.parse()
        self.assertDictEqual(ans['atomnos'],{1: 6, 2: 6, 3: 6, 4: 1, 5: 1, 6: 1, 7: 1})

    def test_allene_dynamics_withCIMat(self):
        with open('glogpy/tests/qdyn_allene.log', 'r') as f:
            data = f.read()
        gj = dj(data)
        ans = gj.parse(do_CI_States=True)
        # print(ans)
        self.assertDictEqual(ans['atomnos'],{1: 6, 2: 6, 3: 6, 4: 1, 5: 1, 6: 1, 7: 1})

    # TEST IV : L118 tests
    def test_formaldehyde_l118(self):
        gj = l118_job('glogpy/tests/methanal_s1_ts.log')
        ans = gj.parse(print_info=False)

    def test_formaldehyde_l118_partial(self):
        gj = l118_job('glogpy/tests/methanal_s1_partial.log', allow_partial=True)
        ans = gj.parse(print_info=False, allow_truncate=True)

    def test_isop_l118_spindens(self):
        gj = l118_job('glogpy/tests/isoprop_sd.log')
        ans = gj.parse(print_info=False)

if __name__ == '__main__':
    unittest.main()