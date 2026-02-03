% Demo: compare Stokeslet traction MEX against MATLAB reference

rng(1);
nt = 10000;
ns = 1000;

xt = rand(nt,3);
xs = rand(ns,3);
nnt = rand(nt,3);
% normalize normals
nnt = nnt ./ vecnorm(nnt,2,2);

qs_3D = rand(ns,3); % source strengths
qs = qs_3D';

% Evaluate Stokeslet traction at targets via direct sum
t_mex = tic;
U = SE0P_Stokestraction_direct_full_ext(xs, qs_3D, nnt, struct('eval_ext_x', xt));
tmex = toc(t_mex);
U = U';
res = 1/(8*pi)*U(:);

% Reference traction matrix
t_ref = tic;
T = stokes_traction(xs, xt, nnt); % source at xs, target at xt
tref = toc(t_ref);
U2 = T*qs(:);

fprintf('Max error in Stokeslet traction computation %1.3e\n', max(abs(res-U2)));
fprintf('Timing: mex %1.3e s, reference %1.3e s\n', tmex, tref);
