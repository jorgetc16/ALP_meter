# ALP_meter 

This is ALP_meter, a code designed to search for periodic signals in astrophysical polarimetry data and set constraints to Ultralight Dark Matter (ULDM) coupling to photons.

ALP_meter is free software, you can redistribute it and/or modify it under the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE

## Main features

ALP_meter has two main features:

- Calculate the Lomb-Scargle periodogram and find significative peaks (MainFindPeaks.py).

- Set cosntraints on periodic signals coming from the interaction between ULDM and photons (MainConstraints.py)
## Data

The data should be in a file of format .dat and organised in three columns: Date of measurement (days) - Polarisation measurement (deg) - 1 sigma uncertainty (deg).

In the folder InputData 10 data sources can be found as example of data files. These are not real observations but simulations based on the distributions followed by real observations.

## How to cite

- Andrés Castillo, Jorge Martin-Camalich, Jorge Terol-Calvo et al.: Searching for dark-matter waves with PPTA and QUIJOTE pulsar polarimetry arxiv.org/abs/2201.03422

## Contributors 

- Andrés Castillo

- Jorge Martin Camalich

- Jorge Terol Calvo
