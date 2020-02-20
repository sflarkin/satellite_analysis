# satellite_analysis (need to think of what to actaully call it.)

Developed by Sean Larkin

This is the (instert name), a collection of pipeline and analysis scripts for analyzing population of galaxies in the 
Generation 3 and Generation 6 VELA simulations, run by Daniel Ceverino. These scripts were designed to be used via terminal
on the NAS pleiades and NERSC Edison supercomputers. 

## Table of Contents

1. [Running Rockstar on the VELA Simulations](#Running-Rockstar)

2. [Accessing my Rockstar Runs on Pleiades](#Accessing-my-Rockstar-Runs-on-Pleiades)

3. Catalog Readers
- Mass Relations

4. Consistent Scripts






## Running Rockstar 

### Art to Ascii

[Rockstar](https://bitbucket.org/gfcstanford/rockstar/src/master/README.md), developed by Peter Behroozi, is a 6-D Friends of
Friends halo finder. While rockstar had been used to analyze the Bolshoi Planck cosmological simulations, another ART simulation,
the formatting of the two art formats were different. As a result, rockstar could not be directly run on the VELA simulations,
as their ART headers varied signifigantly. In order to get around this, I used the yt-project to reformat the VELA particle data
to the basic ASCII format that works with rockstar.

To do this, the script [arttoascii.py](VELAscripts/arttoascii.py) was used.






## Accessing my Rockstar Runs on Pleiades

The Rockstar runs I completed on pleiades will soon be available for download from the NAS data team site (list site here).
Available are the Generation 3 VELA 6 through 15 runs, and all of the Generation 6 runs will be available soon, as the runs
themselves finish in the coming weeks. Each catalog is around 2 GB per simulation.

(need to add a discussion of where and what everything is once I have it formatted correctly for the NASA people)






