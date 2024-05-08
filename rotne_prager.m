function D_hydro = rotne_prager(particle_position,R,D_without_hydro) %#codegen
%D_hydro:
 %       going along the rows or columns: particle 1 x, particle 1 y, particle 2 x, particle 2 y...
%Defining the diffusion matrix using hydrodynamic interactions
%(Rotne-Prager) for a 2D matrix. The matrix is build as X1 Y1 X2 Y2 .... in the row and
%column. It is symmetrical so we will calculate only one triangle of the
%matrix and multiply it into a full square in the end. No interaction
%between different axis of the same particle.


% These are the coefficients taken from Rotne-Prager. We work with D_without_hydro, and Radius in um.
c1 = 3*R/4; %um
c2 = R^3/2; %um

% Create an initial I matrix for all particles. the diagonal is 1, because
% the interaction of the same particle on the same axis with itself is 1.
D_hydro = eye(length(particle_position));

% We will run on all rows up to the diagonal in columns to create the lower
% triangle
for row_particle=1:length(particle_position)
    for column_particle=1:row_particle
        
        %If we have the same particle our matrix won't change (1 for same
        %axis, 0 for different axis).
        same_particle_check = row_particle == column_particle || column_particle+1 == row_particle && mod(row_particle,2) == 0;
        if ~same_particle_check
            
            % this way we get the particle number so we can find both
            % its x axis and y axis
            first_particle = floor((row_particle / 2)+0.5);
            second_particle = floor((column_particle / 2)+0.5);
            
            % get the x, y axis of both particles
            dist_by_x_axis = abs(particle_position(first_particle*2-1) - particle_position(second_particle*2-1));
            dist_by_y_axis = abs(particle_position(first_particle*2) - particle_position(second_particle*2));
            
            % calculate the real distance
            dist_real = sqrt(sum(dist_by_x_axis^2+dist_by_y_axis^2));
            
            
            %for different axis
            if mod(row_particle,2) ~= mod(column_particle,2)
                D_hydro(row_particle,column_particle) = c1*(dist_by_x_axis*dist_by_y_axis/dist_real^3) - c2*(3*dist_by_x_axis*dist_by_y_axis/dist_real^5);
                
                %for y axis
            elseif mod(row_particle,2) == 0
                D_hydro(row_particle,column_particle) = c1*((1/dist_real)+(dist_by_y_axis^2/dist_real^3)) + c2*((1/dist_real^3)-(3*dist_by_y_axis^2/dist_real^5));
                
                %for x axis
            elseif mod(row_particle,2) == 1
                D_hydro(row_particle,column_particle) = c1*((1/dist_real)+(dist_by_x_axis^2/dist_real^3)) + c2*((1/dist_real^3)-(3*dist_by_x_axis^2/dist_real^5));
                
                %if we had a problem with all conditions
            else
                msg = 'No condition fits the chosen row, and column of the particle.';
                error(msg)
                
            end
            
        end
    end
end

% Multiply by the D from Einstein relation (the Rotne-Prager
% calculation is for D_hydro/D_without_hydro).
D_hydro = D_hydro.*D_without_hydro;

% From the triangle create the whole matrix, remove the diagonal
% because it gets doubled
D_hydro = D_hydro'+D_hydro-diag(diag(D_hydro));