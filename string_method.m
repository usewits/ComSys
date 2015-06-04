function result = distance(a,b)
    result = norm(a-b,2);
endfunction

function result = potential(r)
    result = (1/r)^4 - 2*(1/r)^2;
endfunction

function result = dpotentialdx(x,y) %Differentiated to x
    rsq = x*x + y*y;
    result = 4*x*(-1+rsq)/(rsq^3);
endfunction

function result = energy(a)
    result = 0;
    for ix = 1:2:length(a)
        iy = ix+1;
        for jx = (ix+2):2:length(a)
            jy = jx+1;
            result += potential(distance([a(ix),a(iy)],[a(jx),a(jy)]));
        endfor
    endfor
endfunction

function result = grad_energy(a)
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
            result(iy) += dpotentialdx(dy, dx);
        endfor
    endfor
endfunction

function result = t_hat(phi, i)
    result = (phi{i+1}-phi{i})/norm(phi{i+1}-phi{i}); %Does this kind of normalisation make any sense?
endfunction

function result = perp_component_grad_phi(phi, i)
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
endfunction %Tested, works!

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

endfunction %Function works (at least on initial/default phi)

function result = reparametrize(phi)
    d_max = -Inf;
    d_min =  Inf;
    d_mid =  0;    %TODO: test this function
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
    
    while d_max - d_min >= epsilon
        d_mid = (d_min+d_max)/2;
        
        frac_index = 1;
        for i = 1:length(phi)-1%Get evenly spaced points along phi
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
    
    indices = 1:length(phi);
    frac_index = 1;
    for i = 1:length(phi)-1%Get evenly spaced points along phi
        indices(i) = frac_index;
        frac_index = get_fractional_index_at_dist(phi, frac_index, d_mid);
    endfor
    indices(1) = 1; %Force the endpoints to have correct indices
    indices(length(phi)) = length(phi);
    indices
    
    result = cell(0,length(phi));%Obtain positions along phi to return
    for i = 1:length(phi)
        result{i} = phi_at_pos(phi, indices(i));
    endfor
    
endfunction



N = 100;          %Number of positions in string
h = sqrt(3)/2;   %Height of unit-triangle
b = 1/2;
n_iters = 4;

state_a = [ -1, 0,  -b, h,   b, h,   1, 0,   b,-h,   -b,-h,   0, 0 ];
state_b = [ -1, 0,   0, 0,  -b, h,   1, 0,   b,-h,   -b,-h,   b, h ]; 
                  %2 -> 7,  3 -> 2,                          7 -> 3

phi_hist = cell(0,n_iters);
phi = cell(0,N);
for phii = 0:N
    alpha = phii/N;
    phi{phii+1} = alpha*state_a + (1-alpha)*state_b;
endfor

fileID_pos = fopen('pos.txt','w');
fileID_energy = fopen('energy.txt','w');

for iter = 0:n_iters
    dphi = cell(0,N);
    t_hats = cell(0,N);

    for phii = 1:(N-1)
        dphi{phii+1} = perp_component_grad_phi(phi, phii+1);
        t_hats{phii+1} = t_hat(phi, phii+1);
    endfor

    phi_hist{iter+1} = phi;
    write_energies(phi, fileID_energy);

    for phii = 1:(N-1)
        cur_t_hat = t_hats{phii+1}; %This is what is should be!
        cur_dphi = dphi{phii+1};
        cur_phi = phi{phii+1};
        write_to_file(cur_phi, fileID_pos);
        phi{phii+1} -= 0.001*cur_dphi;  %TODO: some GC method moving along this comp
    endfor
    %phi = reparametrize(phi);
endfor

fclose(fileID_pos);
fclose(fileID_energy);
