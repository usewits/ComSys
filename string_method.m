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

N = 10;          %Number of positions in string
h = sqrt(3)/2;   %Height of unit-triangle
b = 1/2;
n_iters = 100;

state_a = [ -1, 0,  -b, h,   b, h,   1, 0,   b,-h,   -b,-h,   0, 0 ];
state_b = [ -1, 0,   0, 0,  -b, h,   1, 0,   b,-h,   -b,-h,   b, h ]; 
                  %2 -> 7,  3 -> 2,                          7 -> 3

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

    write_energies(phi, fileID_energy);

    for phii = 1:(N-1)
        cur_t_hat = t_hats{phii+1}; %This is what is should be!
        cur_dphi = dphi{phii+1};
        cur_phi = phi{phii+1};
        write_to_file(cur_phi, fileID_pos);
        phi{phii+1} -= 0.01*cur_dphi;  %TODO: some GC method moving along this comp
    endfor
endfor

fclose(fileID_pos);
fclose(fileID_energy);
