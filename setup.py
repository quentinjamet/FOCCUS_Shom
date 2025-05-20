from setuptools import setup

setup(name='FOCCUS-Shom',
      version='0.0.1',
      description='Routines used in FOCCUS project at Shom',
      url='https://github.com/quentinjamet/FOCCUS-Shom',
      author='qjamet',
      author_email='quentin.jamet@shom.fr',
      license='MIT',
      packages=['FOCCUS-Shom'],
      install_requires=['numpy', 'scipy', 'xarray', 'datetime'],
      zip_safe=False)
