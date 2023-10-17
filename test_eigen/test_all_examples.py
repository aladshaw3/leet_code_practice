## Pytest for all Eigen Examples Object in C++

import pytest
import sys
from subprocess import run

# -----------------------------------------------------------------------------
# Start test for Ex1
class Test_Ex1:
    @pytest.fixture(scope="class")
    def case(self):
        run_cmd = ['./run_ex1']

        return run_cmd

    def test_run_exe(self, case):
        p = run(case)
        assert(p.returncode==0)

# Start test for Ex2
class Test_Ex2:
    @pytest.fixture(scope="class")
    def case(self):
        run_cmd = ['./run_ex2']

        return run_cmd

    def test_run_exe(self, case):
        p = run(case)
        assert(p.returncode==0)


# Start test for Ex3
class Test_Ex3:
    @pytest.fixture(scope="class")
    def case(self):
        run_cmd = ['./run_ex3']

        return run_cmd

    def test_run_exe(self, case):
        p = run(case)
        assert(p.returncode==0)

# Start test for Ex4
class Test_Ex4:
    @pytest.fixture(scope="class")
    def case(self):
        run_cmd = ['./run_ex4']

        return run_cmd

    def test_run_exe(self, case):
        p = run(case)
        assert(p.returncode==0)

# Start test for Ex5
class Test_Ex5:
    @pytest.fixture(scope="class")
    def case(self):
        run_cmd = ['./run_ex5']

        return run_cmd

    def test_run_exe(self, case):
        p = run(case)
        assert(p.returncode==0)


# Start test for Ex6
class Test_Ex6:
    @pytest.fixture(scope="class")
    def case(self):
        run_cmd = ['./run_ex6']

        return run_cmd

    def test_run_exe(self, case):
        p = run(case)
        assert(p.returncode==0)


# Start test for Ex7
class Test_Ex7:
    @pytest.fixture(scope="class")
    def case(self):
        run_cmd = ['./run_ex7']

        return run_cmd

    def test_run_exe(self, case):
        p = run(case)
        assert(p.returncode==0)

# Start test for Ex8
class Test_Ex8:
    @pytest.fixture(scope="class")
    def case(self):
        run_cmd = ['./run_ex8']

        return run_cmd

    def test_run_exe(self, case):
        p = run(case)
        assert(p.returncode==0)
