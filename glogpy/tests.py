import unittest
from glogpy.freqency_job import frequency_job as fj
from glogpy.dynamics_job import dynamics_job as dj

class TestJobs(unittest.TestCase):

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

    def test_gly_dynamics(self):
        with open('glogpy/tests/qdyn_gly.log', 'r') as f:
            data = f.read()
        gj = dj(data)
        ans = gj.parse()

        self.assertDictEqual(ans['atomnos'],{1: 6, 2: 8, 3: 8, 4: 6, 5: 7, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1})
        print(ans)
    def test_allene_dynamics(self):
        with open('glogpy/tests/qdyn_allene.log', 'r') as f:
            data = f.read()
        gj = dj(data)
        ans = gj.parse()
        self.assertDictEqual(ans['atomnos'],{1: 6, 2: 6, 3: 6, 4: 1, 5: 1, 6: 1, 7: 1})

if __name__ == '__main__':
    unittest.main()