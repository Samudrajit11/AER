function[frame]=rotne_get_frames(particle_count,KbT,Total_time,collision_time,camera_t,repos_time,drivin_dist,pos,R,D,A) %#codegen
%tic
%particle_count = 2;
time_b = round(camera_t/collision_time); % every time_b runs we take a frame
time_a = round(repos_time/collision_time); % every time_a trap reposition
frame = nan(round((Total_time/camera_t)),2*particle_count);
initiate_trap = pos; % For traps
trap=zeros(size(pos));
ax_num=1;
i=1;
while i<(Total_time/collision_time)
    D_new = rotne_prager(pos,R,D);
    A_new = chol(D_new,'lower');

    % Calculate the displacement given the random noise
    random_noise_displacement = (A_new*randn(2*particle_count,1))'.*sqrt(2*collision_time);

    % Calculate the trap force, and total displacement due to trap
    trap(1,1:2*particle_count) =  A .* (pos-initiate_trap);
    trap_displacement = (D_new*trap')'.*(1/KbT).*collision_time;

    % Calculate the new position after one dt
    pos = pos  + random_noise_displacement - trap_displacement;

    if ~mod(i,time_a)
        % reposition the trap with the relevant conversion between pixel and um using some chosen sigma (~size of displacement)
        initiate_trap(ax_num) = 0.0364*randn*drivin_dist*10^(-6);
%        repos_place(i/time_a) = initiate_trap(ax_num); % save the new positions for analysis later
    end

    % Save the position as a frame taken by the camera
    if ~mod(i,time_b)
        frame(i/time_b,1:2*particle_count) = pos;
    end

    i=i+1;

end
%toc
end