# Stacking deep ALMA observations of robust sources above redshift 8

## Background

In an effort to reveal a composite dust detection above redshift 8, this page provides stacking code that can combine multiple observations in the image plane. As a community, this could enable the study of dust production limits in the early Universe from archival and sample-wide z > 8 studies. For convenience, the code is written in Python, with hardcoded source descriptions for convenience. The additional dependencies are numpy, matplotlib, astropy, and spectral_cube. 

## Stacking code
The data is stored in the 'cont'(inuum) folder. 
1. The first step is the python 'weightingFactors.py'. This code reprojects all data to the same cell size and shape, to ensure correct stacking. It also uses the redshift and stellar mass to convert the Jy / beam to M_dust / beam and M_dust / Mstar / beam. The code can be adapted by changing the lists 'sourceName', 'fitsfiles', 'stellarMass', 'redshiftValues', 'magnificationFactor', 'freq', 'radius', and 'betaUV'. Sources can be masked using the 'sourceMask' list. The files are written to the folders 'reprojectCont', 'dustMassFits' and 'dustToStellarMassFits'.
2. The python code stack/stack_things.py subsequently performs the stack preducing images and fits files, according to the classifications in the text files stored in stack. The individual text files describe the source files to stack (all.txt), and any of the sub-classifications. 

## Referencing the data
This project uses data from different groups to study the advent of dust in the Universe. It thus contains many separate efforts, including through the production of stacking techniques ([Lindroos et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.3502L/abstract), [Jolly et al. 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.3992J/abstract)) and the data itself. Please refer to the manuscript from Bakx et al. in prep., for more detailed citations.

## Getting in touch
Beyond referencing previous work, the code and data are free to use. If you wish to get in touch, I am more than happy to participate in discussions and further the studies of dust in the early Universe. Find me at tombakx (a) chalmers.se.