# dust
This code calculates the interstellar dust attenuation coefficient of ELG luminosities in [MultiDark-Galaxies](https://www.cosmosim.org/cms/documentation/projects/galaxies/). It is based on the work of Favole et al. 2019b (in prep.) for SDSS Halpha emitters at mean redshift z=0.1. A similar implementation for [OII] ELGs is documented in [Favole et al. 2019a](https://ui.adsabs.harvard.edu/abs/2019arXiv190805626F/abstract).

The MultiDark-Galaxies mocks have been run with three different semi-analytic models of galaxy formation and evolution: SAG, SAGE and Galacticus. In these SAMs most star formation occurs in the disc, that is where the gaseous component is located, thus we implement dust attenuation only for the galaxy disc component.

In order to run dust.py, you need to import **numpy**, **pytfits** and **astroML.cosmology**. 

