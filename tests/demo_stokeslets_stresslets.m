% Test evaluation of stokeslet and stresslet source

nt = 10000;
ns = 1000;

%nt = 2;
%ns = 2; 

test_stokeslet = 1; 
test_stresslet = 1; 


xt = rand(nt,3);
xs = rand(ns,3);
nns = rand(ns,3);
qs_3D = rand(ns,3); % Generate random source strengths
qs = qs_3D';

% Evaluate Stokeslet flow at targets via direct sum
if test_stokeslet
    t_mex = tic;
    U = SE0P_Stokeslet_direct_full_ext_mex(xs, qs_3D, struct('eval_ext_x', xt));
    tmex = toc(t_mex);
    
    U = U';              % Transpose to match output stacking
    res = 1/(8*pi)*U(:); % Apply Stokeslet prefactor
    
    % Compute reference flow field 
    t_ref = tic;
    G = stokes_SLP(xs,xt);
    tref = toc(t_ref);
    U2 = G*qs(:); 

    fprintf('Max error in Stokeslet computation %1.3e\n',max(res-U2))
    fprintf('Timing: mex %1.3e s, reference %1.3e s\n', tmex, tref);
end

if test_stresslet
    %Same test for stresslet 
    t_mex = tic;
    U = SE0P_Stresslet_direct_full_ext_mex(xs, qs_3D, nns, struct('eval_ext_x', xt));
    tmex = toc(t_mex);
    U = U';
    res = 1/(8*pi)*U(:); 

    % Compute reference flow field. Swapping the roles of source and target
    % in traction computation results in DLP
    t_ref = tic;
    %T = stokes_traction(xt, xs, nns); %source at xt, target at xs   
    %T = T'; %same thing as DLP
    T = stokes_DLP(xs, xt, nns);   %source at xs, target at xt
    tref = toc(t_ref);
    U2 = T*qs(:);
    fprintf('Max error in stresslet computation %1.3e\n',max(res-U2))
    fprintf('Timing: mex %1.3e s, reference %1.3e s\n', tmex, tref);
end

