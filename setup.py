from setuptools import setup

setup(name='massoc',
      version='0.3.0',
      packages=['massoc', 'massoc.GUI', 'massoc.scripts'],
      package_dir={'massoc': 'massoc',
                   'massoc.GUI': 'massoc/GUI',
                   'massoc.scripts': 'massoc/scripts'},
      description='Microbial association network platform',
      author='Lisa Röttjers',
      author_email='lisa.rottjers@kuleuven.be',
      url='https://github.com/ramellose/massoc',
      license='Apache-2.0',
      include_package_data=True,
      summary='Analysis platform for microbial association networks.',
      entry_points={
          'console_scripts': [
              'massoc = massoc.scripts.run_massoc:main'
          ]
      },
      )
