function tr = sim_to_tr(frame,particle_num)
L = length(frame);
tr = zeros(L*particle_num,6);
frame2 = [frame(:,1:2:end), frame(:,2:2:end)];
tr(:,1:2) = reshape(frame2,[L*particle_num,2])*10^6/0.0582;
tr(:,3:4) = ones(L*particle_num,2);

tr(:,5) = repmat(cumsum(ones(L,1)),particle_num,1);
tr_6 = ones(L,particle_num);
for i = 1:particle_num
    tr_6(:,i) = tr_6(:,i)*i;
end
tr(:,6) = tr_6(:);
