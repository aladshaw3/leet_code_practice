[![Matlab](https://github.com/aladshaw3/leet_code_practice/actions/workflows/matlab-script.yml/badge.svg)](https://github.com/aladshaw3/leet_code_practice/actions/workflows/matlab-script.yml)
[![Eigen](https://github.com/aladshaw3/leet_code_practice/actions/workflows/eigen-examples.yml/badge.svg)](https://github.com/aladshaw3/leet_code_practice/actions/workflows/eigen-examples.yml)

# leet_code_practice

This is a repo full of practice problems and other stuff I want to save

**NOTE on Eigen Tests**

In order to better maintain the tests, I have removed Eigen as a submodule and instead build the tests against the latest `master` branch of the Eigen repo in CI.

To run tests locally, you need to have `eigen` in a folder parallel with this repo.

```
	~/projects
	  |
	   -- eigen
	  | 
	   -- leet_code_practice
```

Both `eigen` and `leet_code_practice` should have same parent directory (here named `projects`). Install `eigen` by cloning from GitLab.

```
git clone https://gitlab.com/libeigen/eigen.git
```
