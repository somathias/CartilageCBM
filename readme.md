This repository contains code to simulate a center-based model of cartilage formation and growth. This code was used to generate results and figures for the manuscript ["Contributions of cell behavior to geometric order in embryonic cartilage"](https://www.biorxiv.org/content/10.1101/2022.06.27.497736). 


The code is organized as a [user project](https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides/UserProjects) for the simulation software [Chaste](https://www.cs.ox.ac.uk/chaste/ "Chaste"). In particular, it works with Chaste version 2019.1. 

Once Chaste is set up, this repository can be cloned and linked into the projects folder in the Chaste source, i.e. at `/path/to/Chaste/projects/`. (Re-)configuring Chaste using 

    ccmake path/to/Chaste
    
inside the Chaste build folder should then enable generation of the relevant executables using 

    cd chaste_build/
    make MesenchymalCondensationSimulation
    make CartilageSheetSimulation
    make MissingColumnSimulation

The jupyter notebook files found in the `experiments` folders show how to run the simulations.

The simulation data used to generate the result figures (~2GB) can be downloaded from figshare at [https://doi.org/10.6084/m9.figshare.21731804](https://doi.org/10.6084/m9.figshare.21731804). If unpacked and copied to `experiments/data`, the jupyter notebooks can be used to generate the figures found in the *Results* section of the article. The notebooks also list the exact commands used to generate the data including all parameter settings. 
