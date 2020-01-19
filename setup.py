import setuptools


setuptools.setup(
  name = 'GraphHierarchy',
  packages = ['GraphHierarchy'],   
  version = '0.1',      
  license='MIT',        
  description = 'A module calcualting quantities related to a network metric known as graph hierachy, see Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). Graph hierarchy and spread of infections. arXiv preprint arXiv:1908.04358.',   # Give a short description about your library
  author = 'Choudhry Shuaib',                   
  author_email = 'cshuaib@outlook.com',      
  url = 'https://github.com/shuaib7860/GraphHierarchy',   
  download_url = 'https://github.com/shuaib7860/GraphHierarchy/archive/v0.1.tar.gz',    
  keywords = ['Graph Hierarchy', 'Trophic Coherence', 'Hierarchical Coherence', 'Trophic Levels'],
  install_requires=[
          'numpy',
          'scipy',
          'networkx'
          'statistics'
      ],
  classifiers=[
    'Development Status :: 5 - Production/Stable',     
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Physics',
    'Topic :: Scientific/Engineering :: Mathematics',
    'Topic :: Scientific/Engineering :: Information Analysis',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
  ],
)