from setuptools import setup

setup(name='FOCCUS_Shom',
      version='0.0.1',
      description='Routines used in FOCCUS project at Shom',
      url='https://github.com/quentinjamet/FOCCUS_Shom',
      author='qjamet',
      author_email='quentin.jamet@shom.fr',
      license='MIT',
      packages=['FOCCUS_Shom'],
      install_requires=['numpy', 'scipy', 'xarray', 'datetime'],
      zip_safe=False)
