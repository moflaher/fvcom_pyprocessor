from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()

setup(name='fvcom_pyprocessor',
      version='0.1',
      description='A Python package for the creation of FVCOM input files and related tasks',
      long_description=readme(),
      url='https://github.com/moflaher/fvcom_pyprocessor',
      author='Mitchell O\'Flaherty-Sproul',
      author_email='073208o@acadiau.ca',
      packages=find_packages(),
      zip_safe=False)
