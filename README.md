This directory provides matlab files to simulate the thermal motion 
of many particles interacting via hydrodynamic interactions and 
hard-core repulsion. The motion is considered in two dimensions 
and the Rotne-Prager approximation is used for the hydrodynamic 
interactions. Based on the generated trajectories, the area enclosing 
rate is computed.

For details see:
S. Thapa, D. Zaretzky, R. Vatash, G. Gradziuk, C. Broedersz, Y. Shokef & Y. Roichman,
Observing non-equilibrium probability currents in optically driven colloidal suspensions
https://arxiv.org/pdf/2310.12718

Files:

sim_skel.m
-main file to run. Outputs the AER as a function of imaging rate.

rotne_get_frames.m
-main function that solves the Langevin equation with hydrodynamic interactions
and outputs positions of the particles.

rotne_prager.m
-function that computes the Diffusion matrix in case of Rotne-Prager interaction.

output_aer_2part_fps.m
-function that computes the AER from the positions and gives the AER as a function of 
the imaging rate.


