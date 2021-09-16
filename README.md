# VELAHalos

Developed by Sean Larkin

This is the VELAHalos, a collection of pipeline and analysis scripts for studying galaxies in the Generation 3 and Generation 6 VELA simulations, run by Daniel Ceverino. These scripts were designed to be used via terminal on the NAS pleiades and NERSC Edison supercomputers. 

# TO DO
-Update the GEN3snaps.py file to be more customizable on call with terminal?

## Table of Contents

1. [Running Rockstar on the VELA Simulations](#Running-Rockstar)
    - Finding the VELA files to analyze
    - Art to Ascii
    

2. [Accessing my Rockstar Runs on Pleiades](#Accessing-my-Rockstar-Runs-on-Pleiades)

3. [Catalog Readers](#Catalog-Readers)
- Mass Relations

4. [Consistent Scripts](#Consistent-Scripts)

## Common Terms

VELA_dir: The location where the VELA simulation you want to analyze is located. This should be the folder with the 10MpcBox, PMcrda0, PMcrs0a0, and stars_a0 files are.

consistent_dir: The location where the consistent trees hlist output files are. Should be in the form path/to/hlists.

out_dir: The location where you want the outputs (images, rockstar catalogs, etc.) placed.

VELA_number: The number of the VELA simulation you are working with currently. For example, 07 for VELA07, or 14 for VELA14.

## Running Rockstar 

### Finding the VELA files to Analyze

For the Generation 3 VELA simulations, all of the timesteps were saved. For my analysis, I only focused on the timesteps closest to delta-a of 0.0025, starting at a of 0.1 and ending at a of 0.5. In order to sort through all of the files, and find the ~30% of total snapshots I was interested in on Pleiades, I used the script [GEN3snaps.py](VELAscripts/GEN3snaps.py).

This script takes the timing logs of all of the GEN3 simulations, and iterates through each of them, finding the snaps closest to the times you are interested in, and then writing all of the simulation numbers and scale factors to an .ascii file, where it also includes a shiftc command to move all of them to a directory at once. To change the timesteps it looks for, change the values of the list generator on line 7 to 

`
list_of_wanted_snaps = [i for i in np.arange(starting_a, ending_a + delta_a, delta_a)]
`

where starting_a is the starting scale factor you are interested in, ending_a is the last snapshot you are interested in, and delta_a is the time resolution you want. 

### Art to Ascii

[Rockstar](https://bitbucket.org/gfcstanford/rockstar/src/master/README.md), developed by Peter Behroozi, is a 6-D Friends of
Friends halo finder. While rockstar had been used to analyze the Bolshoi Planck cosmological simulations, another ART simulation, the formatting of the two art formats are incompatible. As a result, rockstar could not be directly run on the VELA simulations, as their ART headers varied significantly. In order to get around this, I used the yt-project to reformat the VELA particle data to the basic ASCII format that works with rockstar.

To do this, the script [arttoascii.py](VELAscripts/arttoascii.py) was used. To call this script, use the command

```
/path/to/satellite_analysis/VELAscripts/arttoascii.py VELA_dir out_dir
```

I often used the VELA_dir as the out_dir to keep the simulation files together, but if space is a concern, they can be kept seperate for easy management.

### Getting the Stellar Mass and Gas Data

After rockstar is run, the stellarmassrelation.py script can be run on the consistent trees data to extract the stellar mass and gas mass for each of the galaxies found. This script find the stellar mass within 3 multiples of the Virial Radius, .1, .15, .2 and 1.0 of the RVir. It also finds the gas within 1 Rvir. The process for finding the gas mass is very time intensive, due to yt needing to deconstruct the oc-tree of the VELA sims, so while other radii can be looked into, I have not yet done such analysis. 

This script is an early draft used for some prelimary investigations, but is not the most accurate for the VELA simulations. I reccomend using the post processing scripts included here.


## Catalog Readers

The most called programs I have made are my catalog readers. These programs are usually the first thing called in my analysis scripts. They are used to parse through the entire catalog of halos for a given VELA, and get all of the halos of the wanted parameters, usually for graphing. 

### Rockstar and Consistent Catalog Readers

These two script do the same process for the two different catalogs made by the rockstar/consistent trees halo finding. I recommend using the consistent trees catalogs, but if you have only run rockstar and are not planning on using consistent trees, you can use the rockstar reader with a few modifications. NOTE: if you do plan on using the rockstar reader for some of the scripts that only exist for consistent trees, the index for the two catalogs rows differ, so they will need to be changed.

To call the Consistent Catalog Reader use 

```
from satellite_analysis.catalogreaders import consistentcatalogreader as consistent

consistent.conistent_catalog_reader(consistent_dir, remove_subhalos='False', halo_mass=mass)
```

The remove_subhalos and halo_mass are optional arugments. If you want the reader to remove all of the subhalos from the halo data, set remove_subhalos='True', the default is 'False', which keeps the subhalos in. The halo_mass is the mass in MSun that the reader will remove halos smaller than that size. The default is set to 10e+08.

This script generates a list of lists, with length of the total number of snaps for a given simulation. Each list is contains a series of tuples, with each tuple including all of the data included in the consistent trees hlist catalogs, sorted by halo mass, with the largest halo found in the first tuple, and so on. To call the largest halo from the third snapshot, use

```
consistent.halo_data_sorted[3][0]
```

If you want a list of all of the scale factors for each snapshot, use 

```
consistent.consistent_file_index
```

### Tomer Catalogs

For several of my scripts, I do a comparison to others analysis of the VELA simulations. The data for these comparisons comes from Tomer Nessbaum's pandas catalogs. They contain data gathered by Nir Mandelkir, Sharon Lapiner, Tomer Nessbaum, and others. For more info on this catalog, contact Sharon Lapiner.

To call this catalog use

```
import pandas as pd
from satellite_analysis.catalogreaders import tomercatalogreader as tomer
```

For more info on what is in the catalog, and to see how to call the info you want from pandas data structure, see the using this catalog that is included in it.


## Rockstar Post Processing Scripts

The rockstar version I used was the single mass rockstar, which can only read particle data if they share the Mass. As the VELA simulations use a variety of particle masses for the stars and dark matter, rockstar was run only on the smallest mass dark matter particles. As a result of this, its calculations of several galaxy properties like the Virial Radius, Virial Mass, etc were smaller than if the stars and gas were included. Because of this, I used the galaxy outputs of rockstar, and ran a suite of post processing scripts to find more accurate halo measurements of several paramaters; Mvir/Rvir, SFR, Higher Mass Dark Matter Contamination, and Metallicity. Below I include a decription of these scripts and their dependancies.

### Proper Mvir Calc

As rockstar underestimates the size of the halos it finds as a result of only looking at the dark matter, to properly analyze the halo population, recalculating the Virial Mass and Radius are of upmost importance. This script takes the halos found in the rockstar catalog larger than a given mass, and recalculates the virial mass. To do this, it creates a profile of the gas, stars, and dark matter of the galaxy out to 2-2.5x times the virial radius found by rockstar. It then finds when the right delta_vir for the current epoch is based on the Virial function from Byron and Norman ???? within 1/5000 of the virial radius of rockstar, which is usually around 0.002 kpc, which is much lower than the error caused by the fitting to a not truly spherical profile or the gas-mesh resolution of 17.5 pc for the inner regions at late times.

After it finds this 'proper' Mvir and Rvir point, it saves the profiles it generated for the gas and stars, and calculates the gas, stars, and dark matter found within 10%, 15%, and 20% of the newly found rvir, and saves them as well.

This also generates images of the 3 2-d projections of the stars and dark matter within each halo, and saves them. The file names tell you about the halo in question, and are formatted in this way: aXXX_Y_ZZZZZZ.0.png where XXX is the scale factor of the current time, Y is the halo ordering based on size of those for this snapshot (0 being the largest, 1 the second largest, 2 the third largest, and so on), and ZZZZZZ being the halo_id as definded by consistent trees. 

Here is what this looks like for the Generation 6 VELA07 Largest Halo at scale factor 0.500.

![Propermvircalc image](READMEfigures/a500_0_id937383.0.png)

To run this script, you need the location of the consistent-trees hlists, and the loaction of the VELA simulation files, and run the following command.

```
/path/to/satellite_analysis/VELAscripts/propermvircalc.py VELA_dir out_dir
```

*NOTE ABOUT Skipping Halos: In my analysis, I found one halo where this script failed for finding its virial radius and mass. I could not find the reason for this, so this script is designed to skip over halos where this process fails. If you want to find any that have failed, you can compare the number of images made to the number of halos in the output file, as the images are created before this discrimination happens. The other scipts in this section will still work if any are skipped.

### SFR and Contamination

Another important property to understanding halo evolution is the Star Formation Rate (SFR). The script here takes the proper rvir values from the previous script, and calculates the SFR and SSFR for several stellar populations and radii. All populations are calculated for those populations withing 10%, 15%, and 20% of the virial radius. The stellar populations considered are as follows:
-The Stars Formed since the last timestep
-The Stars Formed withing the last 2, 3, 4, and 5 timesteps
-The stars formed within 60 and 100 Myr
-The stars formed 5, 60, and 100 Myr according to the formula from Tacchella et al. 2016

In addition to the SFRs, this script also sees if any of the higher mass dark matter particles from the less resolved outer regions of the ART code have contaminated any of the halos. This is used primarily as an exclusion paramater for those halos with high contamination (mass of contaminating particles reaching 10% of the total Virial Mass). 

### Gas Abundance and Metallicity

This script calculates the metallicity and metallicity profile for each halo. Need to clean this one up as some of the fields are wrong.

## SQL Database Creation

For several of my post processing scripts, I use a SQL database 

## Consistent Scripts 

If you do not have access or interest in creating a SQL database for this data, there are a number of scripts that work from the .ascii outputs of the consistent trees and post processing scripts. These are earlier versions of scripts found in the SQL section, so they may not have the full functionality of them, but they preform mostly the same processes.



## Accessing my Rockstar Runs on Pleiades

The Rockstar runs I completed on pleiades will soon be available for download from the NAS data team site (list site here). Available are the Generation 3 VELA 6 through 15 runs, and all of the Generation 6 runs will be available soon, as the runs themselves finish in the coming weeks. Each catalog is around 2 GB per simulation.

(need to add a discussion of where and what everything is once I have it formatted correctly for the NASA people)
