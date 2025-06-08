
from setuptools import setup

scripts = ['pymdpbsa', 'pytleap']

if __name__ == '__main__':

    setup(name='AmberLite',
          version='22.0',
          description='Various modules needed for AmberTools Python programs',
          author='Romain Wolf, Pawel Janowski, and Jason M. Swails',
          author_email='jason.swails -at- gmail.com',
          url='http://ambermd.org',
          license='GPL v2 or later',
          packages=[], py_modules=[],
          scripts=scripts)
