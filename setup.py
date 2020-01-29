import setuptools


LONG_DESCRIPTION = """
**GraphHierarchy** is a python package that calculates the hierarchical levels of nodes in a network as well as hierarchical coherence of a network structure. 
Hierarchical levels are the mathematical generalisation of the trophic analysis of networks. Trophic levels and hence trophic coherence can be defined only on networks with well defined sources and so trophic analysis of networks had been restricted to the ecological domain, until now. 
Graph Hierarchy is a python package that allows for analysis of all network structures via the trophic levels and coherence approach. 
Trophic coherence, a measure of a network’s hierarchical organisation, has been shown to be linked to a network’s structural and dynamical aspects. 
In GraphHierarchy we have developed the python code which implements the mathematical generalisation of the trophic coherence theory to all networks. See citation paper for more details.
.. _GitHub: https://github.com/shuaib7860/GraphHierarchy
"""

setuptools.setup(
  name = 'GraphHierarchy',
  packages = ['GraphHierarchy'],   
  version = '0.3',      
  license='MIT',        
  description = 'A module calculating quantities related to a network metric known as graph hierachy, see Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). Graph hierarchy and spread of infections. arXiv preprint arXiv:1908.04358.',   # Give a short description about your library
  long_description = LONG_DESCRIPTION,
  author = 'Choudhry Shuaib',                   
  author_email = 'cshuaib@outlook.com',      
  url = 'https://github.com/shuaib7860/GraphHierarchy',   
  download_url = 'https://github.com/shuaib7860/GraphHierarchy/archive/v0.3.tar.gz',
  keywords = ['Graph Hierarchy', 'Trophic Coherence', 'Hierarchical Coherence', 'Trophic Levels'],
  install_requires=[
          'numpy',
          'scipy',
          'networkx',
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