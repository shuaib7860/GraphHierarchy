import setuptools


LONG_DESCRIPTION = """
**GraphHierarchy** is a python package that calculates the hierarchical level of nodes in a network, the associated hierarchical differences for the edges and the hierarchical coherence of the network. 
Hierarchical levels are the mathematical generalisation of the trophic analysis of networks. Trophic levels and hence trophic coherence can be defined only on networks with well defined sources, known as basal nodes. 
Trophic coherence, a measure of a network’s hierarchical organisation, has been shown to be linked to a network’s structural and dynamical properties. Thus trophic analysis of networks had been restricted to the ecological domain, until now.  
Graph Hierarchy is a python package that implements this mathematical generalisation that allows for analysis of all network structures via the trophic approach. 
See Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). Graph hierarchy and spread of infections. arXiv preprint arXiv:1908.04358 for more details.
.. _GitHub: https://github.com/shuaib7860/GraphHierarchy
"""

setuptools.setup(
  name = 'GraphHierarchy',
  packages = ['GraphHierarchy'],   
  version = '1.2',      
  license='MIT',        
  description = 'A module calculating quantities related to a network metric known as trophic coherence but now generalised to all networks, see Moutsinas, G., Shuaib, C., Guo, W., & Jarvis, S. (2019). Graph hierarchy and spread of infections. arXiv preprint arXiv:1908.04358 for more details.',
  long_description = LONG_DESCRIPTION,
  author = 'Choudhry Shuaib',                   
  author_email = 'cshuaib@outlook.com',      
  url = 'https://github.com/shuaib7860/GraphHierarchy',   
  download_url = 'https://github.com/shuaib7860/GraphHierarchy/archive/v1.2.tar.gz',
  keywords = ['Graph Hierarchy', 'Trophic Coherence', 'Hierarchical Coherence', 'Trophic Levels', 'Hierarchical Levels', 'Influence Centrality', 'Democracy Coefficient'],
  install_requires=[
          'numpy',
          'scipy',
          'networkx'
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