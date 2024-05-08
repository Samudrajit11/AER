function[outmat12]=output_aer_2part_fps(frame_full)

parttotal=2;
outmat12=[];

for frameind=1:2000
    frame=[];
    frame=frame_full(1:frameind:end,:);
    FPS=10000;
    FPS=FPS/frameind;      
    tr = sim_to_tr(frame,parttotal);
    st = AER_mult(tr,parttotal,FPS);
    out12=st.particle_12(end,1);
    
    %out1=(drivemat(kk)*0.01);

    outmat12=[outmat12;FPS out12];
    
end 

end