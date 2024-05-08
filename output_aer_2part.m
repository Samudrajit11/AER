function[out]=output_aer_2part(frame)

FPS=10000;
parttotal=2;
tr = sim_to_tr(frame,parttotal);
st = AER_mult(tr,parttotal,FPS);
out=st.particle_12(end,1);
%out=[out12 out13 out23];

end