# dust
This code calculates the interstellar dust attenuation coefficient of ELG luminosities in [MultiDark-Galaxies](https://www.cosmosim.org/cms/documentation/projects/galaxies/). It is based on the work of Favole et al. 2019b (in prep.) for SDSS Halpha emitters at mean redshift z=0.1. A similar implementation for [OII] ELGs is documented in [Favole et al. 2019a](https://ui.adsabs.harvard.edu/abs/2019arXiv190805626F/abstract).

The MultiDark-Galaxies mocks have been run with three different semi-analytic models of galaxy formation and evolution: SAG, SAGE and Galacticus. In these SAMs most star formation occurs in the disc, that is where the gaseous component is located, thus we implement dust attenuation only for the galaxy disc component.

In order to run dust.py, you need to import **numpy**, **pytfits** and **astroML.cosmology**. 

An input catalogue *SAGinput_z=0.09.fits* is provided for testing. This is a random sampling from SAG model galaxies at z=0.09. The input columns are: 
1. Non-attenuated Halpha luminosity in log scale: *log10(LHa/erg s^-1)*
2. Cold gas mass in the disc in *[Msun]*
3. Half mass radius of the disc in *[Mpc]*
4. Metallicity of the disc

The output table provides the following columns:
1. Non-attenuated Halpha luminosity in log scale
2. Attenuated Halpha luminosity in log scale
3. Non-attenuated Halpha flux in log scale: *log10(FHa/erg s^-1 cm^-2)*
4. Attenuated Halpha flux in log scale
5. Attenuation coefficient of the disc *AA*, which is used as. *logLatt=logL-0.4AA*

