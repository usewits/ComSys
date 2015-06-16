global potential;
global dpotentialdx;
global dpotentialdy;
global energy;
global grad_energy;


function result = distance(a,b)
    result = norm(a-b,2);
endfunction

function result = t_hat(phi, i)
    result = (phi{i+1}-phi{i})/norm(phi{i+1}-phi{i}); %Does this kind of normalisation make any sense?
endfunction

function result = perp_component_grad_phi(phi, i)
    global grad_energy;
    t = t_hat(phi,i);
    g = grad_energy(phi{i});
    result = g - (g*t')*t;
endfunction

function write_to_file(a, fileID)
    for i = 1:length(a)
        if i != 1
            fprintf(fileID,"\t");
        endif
        fprintf(fileID,"%f",a(i));
    endfor
    fprintf(fileID,"\n");
endfunction

function write_energies(phi, fileID)
    global energy
    for i = 1:length(phi)
        if i != 1
            fprintf(fileID,"\t");
        endif
        fprintf(fileID,"%f",energy(phi{i}));
    endfor
    fprintf(fileID,"\n");
endfunction

function result = phi_at_pos(phi, fractional_index)
    index = floor(fractional_index);
    alpha = fractional_index - index;
    if index >= length(phi)
        if alpha < 1e-10    %NOTE: epsilon!
            result = phi{index};
        else
            disp("Error: phi_at_pos called on invalid fractional_index!")
        endif
    else
        result = phi{index}*(1-alpha) + phi{index+1}*alpha;
    endif
endfunction

function result = get_fractional_index_at_dist(phi, fractional_index_a, dist)
    index_a = floor(fractional_index_a);%The integer index just before fractional_index_a
    index_b = index_a;%Goes wrong if last integer index is more than dist away..

    %Find index_b, such that the segment [phi{index_b}, phi{index_b+1}]
    %contains a point at distance dist
    pos_a = phi_at_pos(phi, fractional_index_a);
    while distance(pos_a, phi{index_b+1}) < dist
        index_b += 1;
        if index_b == length(phi)
            result = length(phi) + 1; %overflow!
            return
        endif
    endwhile

    %find point on [phi{index_b}, phi{index_b+1}] at distance dist from pos_a
    %(intersection of line-segment with circle)
    b = phi{index_b};
    bp = phi{index_b+1};
    a = pos_a;

    ba = b*a';
    bbp = b*bp';
    bpa = bp*a';
    a2 = a*a';
    b2 = b*b';
    bp2 = bp*bp';

    determinant = 4*(ba-bbp+bp2-bpa)^2 - 4*(b2-2*bbp+bp2)*(a2+bp2-2*bpa-dist^2);
    if determinant < 0
        result = -1;
        disp("Error! Could not find intersection! (negative determinant)");
        fractional_index_a
        dist
        disp("!")
        return
    endif
    alpha1 =  1 - (ba-bbp+bp2-bpa-0.5*sqrt(determinant)) / (b2-2*bbp+bp2);
    alpha2 =  1 - (ba-bbp+bp2-bpa+0.5*sqrt(determinant)) / (b2-2*bbp+bp2);

    epsilon = 1e-8;%NOTE: this should be bigger than the other epsilon used
    if alpha1 >= -epsilon && alpha1 < 0
        alpha1 = 0;
        disp("alpha1 = -0, correcting!");
    endif
    if alpha1 > 1 && alpha1 <= 1+epsilon
        alpha1 = 1;
        disp("alpha1 = +1, correcting!");
    endif
    if alpha2 >= -epsilon && alpha2 < 0
        alpha2 = 0;
        disp("alpha2 = -0, correcting!");
    endif
    if alpha2 > 1 && alpha2 <= 1+epsilon
        alpha2 = 1;
        disp("alpha2 = +1, correcting!");
    endif

    if alpha1 >= 0 && alpha1 <= 1
        result = index_b+alpha1;
    elseif alpha2 >= 0 && alpha2 <= 1
        result = index_b+alpha2;
    else
        result = -1;
        disp("Error! Could not find intersection! (no valid alpha)");
        fractional_index_a
        dist
        disp("!")
    endif

endfunction

function result = reparametrize(phi)
    d_max = -Inf;
    d_min =  Inf;
    d_mid =  0;
    for i = 1:(length(phi)-1)
        d = distance(phi{i+1}, phi{i});
        if d < d_min
            d_min = d;
        endif
        if d > d_max
            d_max = d;
        endif
    endfor

    epsilon = 1e-10;
    indices = 1:length(phi);

    while d_max - d_min >= epsilon
        d_mid = (d_min+d_max)/2;
        
        frac_index = 1;
        for i = 1:length(phi)-1%Get evenly spaced points along phi
            indices(i) = frac_index;
            frac_index = get_fractional_index_at_dist(phi, frac_index, d_mid);
            if frac_index > length(phi)
                break
            endif
        endfor

        if frac_index > length(phi) %d_mid is too big!
            d_max = d_mid;
        else %d_mid is too small (or exactly right)!
            d_min = d_mid;
        endif
    endwhile
    
    indices(1) = 1; %Force the endpoints to have correct indices
    indices(length(phi)) = length(phi);
    
    result = cell(0,length(phi));%Obtain positions along phi to return
    for i = 1:length(phi)
        result{i} = phi_at_pos(phi, indices(i));
    endfor
    
endfunction


%Settings of the simulation
N = 10;          %Number of positions in string
n_iters = 100;
stepsize = 0.001;
simulation = "peaks";
%simulation = "particles";



if strcmpi(simulation, "particles")
    %Begin and end state (a and b)
    h = sqrt(3)/2;   %Height of unit-triangle
    b = 1/2;
    state_a = [ -1, 0,  -b, h,   b, h,   1, 0,   b,-h,   -b,-h,   0, 0 ];
    state_b = [ -1, 0,   0, 0,  -b, h,   1, 0,   b,-h,   -b,-h,   b, h ]; 
                      %2 -> 7,  3 -> 2,                          7 -> 3

    %Potential and derivative
    function result = potentialParticles(x,y)
        r = sqrt(x*x + y*y);
        result = (1/r)^4 - 2*(1/r)^2;
    endfunction

    function result = dpotentialdxParticles(x,y) %Differentiated to x
        rsq = x*x + y*y;
        result = 4*x*(-1+rsq)/(rsq^3);
    endfunction

    function result = dpotentialdyParticles(x,y) %Differentiated to x
        rsq = x*x + y*y;
        result = 4*y*(-1+rsq)/(rsq^3);
    endfunction

    function result = energyParticles(a)        %Assuming 2D sum of local potentials
        global potential;
        result = 0;
        for ix = 1:2:length(a)
            iy = ix+1;
            for jx = (ix+2):2:length(a)
                jy = jx+1;
                
                dx = abs(a(ix)-a(jx));
                dy = abs(a(iy)-a(jy));
                result += potential(dx, dy);
            endfor
        endfor
    endfunction

    function result = grad_energyParticles(a)    %Assuming 2D sum of local potentials
        global dpotentialdx;
        global dpotentialdy;
        result = zeros(1,length(a));
        for ix = 1:2:length(a)
            iy = ix+1;
            for jx = 1:2:length(a)
                if ix == jx
                    continue;
                endif
                jy = jx+1;
                dx = abs(a(ix)-a(jx));%is this correct way to calc grad?
                dy = abs(a(iy)-a(jy));
                result(ix) += dpotentialdx(dx, dy);
                result(iy) += dpotentialdy(dx, dy);
            endfor
        endfor
    endfunction

    potential = @potentialParticles;
    dpotentialdx = @dpotentialdxParticles;
    dpotentialdy = @dpotentialdyParticles;
    energy = @energyParticles;
    grad_energy = @grad_energyParticles;
endif

if strcmpi(simulation, "peaks")
    %Begin and end state (a and b) (this are local minima, where found numerically using mathematica)
    state_a = [-1.3474,   0.204519];
    state_b = [ 0.228279,-1.62553 ];

    %Potential and derivative
    function result = potentialPeaks(x,y)
        result = peaks(x,y);
    endfunction

    function result = dpotentialdxPeaks(x,y) %Differentiated to x
        result = -6*exp(-x^2-(1+y)^2)*(1-x)-6*exp(-x^2-(1+y)^2)*(1-x)^2*x+2/3*exp(-(1+x)^2-y^2)*(1+x)-10*exp(-x^2-y^2)*(1/5-3*x^2)+20*exp(-x^2-y^2)*x*(x/5-x^3-y^5);
    endfunction

    function result = dpotentialdyPeaks(x,y) %Differentiated to x
        result = 2/3*exp(-(1+x)^2-y^2)*y+50*exp(-x^2-y^2)*y^4-6*exp(-x^2-(1+y)^2)*(1-x)^2*(1+y)+20*exp(-x^2-y^2)*y*(x/5-x^3-y^5);
    endfunction

    function result = energyPeaks(a)         %Assuming 2D global static potential
        global potential;
        result = 0;
        for ix = 1:2:length(a)
            iy = ix+1;
            result += potential(a(ix), a(iy));
        endfor
    endfunction

    function result = grad_energyPeaks(a)    %Assuming 2D global static potential
        global dpotentialdx;
        global dpotentialdy;
        result = zeros(1,length(a));
        for ix = 1:2:length(a)
            iy = ix+1;
            result(ix) = dpotentialdx(a(ix), a(iy));
            result(iy) = dpotentialdy(a(ix), a(iy));
        endfor
    endfunction

    potential = @potentialPeaks;
    dpotentialdx = @dpotentialdxPeaks;
    dpotentialdy = @dpotentialdyPeaks;
    energy = @energyPeaks;
    grad_energy = @grad_energyPeaks;
endif



%Find intermediate states using linear interpolation
phi = cell(0,N);        %phi{1} and phi{N} are fixed, phi{2}..phi{N-1} are not
for phii = 1:N
    alpha = (phii-1)/(N-1);
    phi{phii} = alpha*state_a + (1-alpha)*state_b;
endfor

%Open files to write results
fileID_pos = fopen('pos.txt','w');
fileID_energy = fopen('energy.txt','w');

phi_hist = cell(0, n_iters);      %Store all string positions (for debugging)

for iter = 0:n_iters
    iter                         %Display iteration number
    dphi = cell(0,N);
    t_hats = cell(0,N);

    for phii = 2:(N-1)             %Calculate direction of potential and string
        dphi{phii} = perp_component_grad_phi(phi, phii);
        t_hats{phii} = t_hat(phi, phii);
    endfor

    phi_hist{iter+1} = phi;      %Store all strings for review
    write_energies(phi, fileID_energy);

    for phii = 1:N
        write_to_file(phi{phii}, fileID_pos);
    endfor

    for phii = 2:(N-1)           %Update the string
        cur_t_hat = t_hats{phii};
        cur_dphi = dphi{phii};
        cur_phi = phi{phii};
        phi{phii} -= stepsize*cur_dphi;  %TODO: some GC method moving along this comp
    endfor

    phi = reparametrize(phi);   %Make sure states in string are equidistant
endfor

fclose(fileID_pos);
fclose(fileID_energy);
