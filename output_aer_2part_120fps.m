function[out]=output_aer_2part_120fps(frame_full)

parttotal=2;
frameind=83;           
frame=frame_full(1:frameind:end,:);
FPS=10000;
FPS=FPS/frameind;      
tr = sim_to_tr(frame,parttotal);
st = AER_mult(tr,parttotal,FPS);
out=st.particle_12(end,1);


end