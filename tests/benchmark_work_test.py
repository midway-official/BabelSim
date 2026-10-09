"""The timed outer passes and accepted physical steps have different meanings."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from benchmark_solver_optimization import piso_work


class WorkTest(unittest.TestCase):
    def test_standard_piso(self):
        log = "\n".join(f"PISO 1 time={t} accepted=true linear=ok converged=true"
                        for t in (0.01, 0.02, 0.03))
        self.assertEqual(piso_work(log), (3, 3))

    def test_nonlinear_and_failed_steps(self):
        log = """PISO 1 time=0.01 accepted=true converged=false
PISO 72 time=0.01 accepted=true converged=true
PISO 1 time=0.02 accepted=true converged=false
PISO 26 time=0.02 accepted=true converged=true
PISO 1 time=0.03 accepted=false converged=false
PISO corrector=1 coupling=0.01
"""
        self.assertEqual(piso_work(log), (2, 98))


if __name__ == "__main__":
    unittest.main()
