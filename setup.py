from setuptools import setup,find_packages
from FLED.__version__ import __version__

setup(name='FLED',
      version=__version__,
      description='full length ecDNA Detection tools',
      author='Lifuyu',
      url='https://github.com/FuyuLi/FLED',
      packages=find_packages(),
      install_requires=[
          'pysam==0.16','networkx==2.5','progressbar==2.5','biopython==1.76','numpy>=1.19.1',
          'pyspoa==0.0.5','scipy>=1.5.3'
      ],
      entry_points={
          'console_scripts': [
              'FLED = FLED.longreadscircDNAdetection:main'
          ],

      },
      classifiers=[
          'License :: .'
      ],
     )
