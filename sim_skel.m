function[out_fps]=sim_skel(dist,drivin_dist)
%
%dist - distance between particles in units of particle radius
%driving_dist - b_o/0.0364 um %[1.5 3.022] gives b_0=[55 110] nm
%out_fps - matrix with column 1- FPS, column 2 - AER
tic
rng('shuffle')
particle_count = 2; % number of particles
A=[10 10 4 4].*10^(-6);  %N/m stiffness matrix [k1x k1y k2x k2y] 
repos=36; %Hz, repositing rate 
Total_time = 600; % sec, is enough with this frame rate (total real time)
collision_time = 1e-5; % sec, time between collisions
camera_t = 1/10000; % sec, time between camera photos (1/FPS)
repos_time = 1/repos; %sec 

Kb = physconst('boltzmann'); % J/K
T = 300; %K, ~room temp
KbT = Kb*T; % J

R = 0.75*1e-6; % m, particle radius
eta = 0.89e-3; % Pa*sec, of water %eta=5e-3
gamma = 6*pi*eta*R; % N*sec/m
D = KbT/gamma; % m^2/sec
distan=dist*R; 
pos = [0.0 0.0 distan 0.0];
[frame]=rotne_get_frames_mex(particle_count,KbT,Total_time,collision_time,camera_t,repos_time,drivin_dist,pos,R,D,A);
frame = frame(any(frame,2),:);

 out_fps=output_aer_2part_fps(frame);

toc
 


end


