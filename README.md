# CSI-Jupyter
Jupyter notebooks for CSI software

Some simple Jupyter notebooks for using the Classic Source Inversion (CSI) software, which has a version at https://github.com/jolivetr/csi

These examples are for a basic slip inversion on the 2019 Ridgecrest earthquake in California, using interferograms available online at the Harvard Dataverse (doi:10.7910/DVN/JL9YMS). See this reference for description of the interferograms:
Fielding, E. J., Z. Liu, O. L. Stephenson, M. Zhong, C. Liang, A. Moore, S.-H. Yun, and M. Simons (2020), Surface deformation related to the 2019 Mw 7.1 and 6.4 Ridgecrest Earthquakes in California from GPS, SAR interferometry, and SAR pixel offsets, Seismol. Res. Lett., 91(4), 2035-2046, doi:10.1785/0220190302.

The Jupyter notebooks should be run starting with the covarSAR-A2D166.ipynb, then the downsampleSAR-A2D166.ipynb, and finally the slipInversion-A2D166.ipynb 

The notebooks have also been exported to Python scripts that can be run on the command line.