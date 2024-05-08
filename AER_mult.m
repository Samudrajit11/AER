function st = AER_mult(tr,NOF,FPS)

%;------
%; Calculate AER for the given FPS, for all particle combinations.
%; output is the a struct, containing (for every pair),
%; a [nm^2/msec] vector for XX and for YY.
%; tr: is just the regular tr for all of the particles.
%; FPS: is the FPS
%; NOF: number of particles in the tr.
%;------


%Create a matrice of all possible combinations of 2 particles
C = nchoosek(1:NOF,2);

for i=1:size(C,1)
    
    %Create the relevant tr for the AER (X,Y for the relevant particles)
    m_a = find(tr(:,6)==C(i,1)); m_b = find(tr(:,6)==C(i,2));
    tr_a = tr(m_a,1:2); tr_b = tr(m_b,1:2);
    tr2 = [tr_a tr_b];
    
    %Run the AER
    [XX,YY,XY,YX] = calc_AER(tr2,FPS);
    T = cumsum(ones(length(XX),1)).*1/FPS;
    %Create a fieldname according to the particles
    A = char(strjoin(string(C(i,:))));
    fieldname = A(~isspace(A));
    %saves into a struct, which will be outputed
    st.('particle_'+string(fieldname)) = [XX YY XY YX T];   
end


end