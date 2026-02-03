% Test evaluation of stokeslet and stresslet source

nt = 10000;
ns = 1000;

nt = 2;
ns = 2; 

test_stokeslet = 1; 
test_stresslet = 1; 


xt = rand(nt,3);
xs = rand(ns,3);
nns = rand(ns,3);
qs_3D = rand(ns,3); % Generate random source strengths
qs = qs_3D';

% Evaluate Stokeslet flow at targets via direct sum
if test_stokeslet
    U = SE0P_Stokeslet_direct_full_ext_mex(xs, qs_3D, struct('eval_ext_x', xt));
    
    U = U';              % Transpose to match output stacking
    res = 1/(8*pi)*U(:); % Apply Stokeslet prefactor
    
    % Compute reference flow field 
    G = stokes_SLP(xs,xt);
    U2 = G*qs(:); 

    fprintf('Max error in Stokeslet computation %1.3e\n',max(res-U2))
end

if test_stresslet
    %Same test for stresslet 
    U = SE0P_Stresslet_direct_full_ext_mex(xs, qs_3D, nns, struct('eval_ext_x', xt));
    U = U';
    res = 1/(8*pi)*U(:); 

    % Compute reference flow field. Swapping the roles of source and target
    % in traction computation results in DLP
    T = stokes_traction(xt, xs, nns); %source at xt, target at xs   
    T = T'; %same thing as DLP
    T = stokes_DLP(xs, xt, nns);   %source at xs, target at xt
    U2 = T*qs(:);
    fprintf('Max error in stresslet computation %1.3e\n',max(res-U2))
end


function G = stokes_SLP(rin, rout)
%STOKES_SLP Builds dense Stokeslet target-from-source matrix for 3D Stokes flow.
%
%   G = STOKES_SLP(rin, rout)
%
%   Constructs the full dense matrix G representing the (single-layer)
%   Stokeslet kernel in 3D, mapping a vector of point forces at source locations
%   `rin` to velocities at target locations `rout`.
%
%   INPUTS:
%       rin  - N x 3 matrix of source points.
%       rout - M x 3 matrix of target points.
%
%   OUTPUT:
%       G    - 3M x 3N block matrix where each 3x3 block represents the Stokeslet
%              kernel G(x, y) between target x = rout(i,:) and source y = rin(j,:).
%              That is, for each pair (i, j),
%              G(x,y) = (I + r⊗r / |r|^2) / (8π|r|),  where r = x - y.
%
%   NOTES:
%       - This matrix is expensive to construct and store for large N.
%         It is mostly useful for verification, small problems, or dense solvers.
%       - The output format is consistent with applying S * f where f is a
%         stacked 3N x 1 force vector.
%
%   EXAMPLE USAGE:
%       rin = rand(10,3);
%       rout = rand(12,3);
%       S = stokes_SLP(rin, rout);
%
%   Anna Broms June, 12, 2025

d1 = bsxfun(@minus,rin(:,1)',rout(:,1)); % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,rin(:,2)',rout(:,2));
d3 = bsxfun(@minus,rin(:,3)',rout(:,3));
r2 = d1.^2+d2.^2+d3.^2;   % dist^2 mat

ir = 1 ./ sqrt(r2);   % 1/r, 
ir3 = ir ./ r2;        	% 1/r^3, 

d1d2 = d1.*d2.*ir3;
d1d3 = d1.*d3.*ir3;
d2d3 = d2.*d3.*ir3;


D1 = zeros(3); D1(1,1) = 1;
D2 = zeros(3); D2(2,2) = 1;
D3 = zeros(3); D3(3,3) = 1;

E1 = zeros(3); E1(1,2) = 1; E1(2,1) = 1;
E2 = zeros(3); E2(1,3) = 1; E2(3,1) = 1;
E3 = zeros(3); E3(3,2) = 1; E3(2,3) = 1;

dxd_ir3 = kron(d1.^2.*ir3,D1) + kron(d2.^2.*ir3,D2) + kron(d3.^2.*ir3,D3) + kron(d1d2,E1) +...
    kron(d1d3,E2) + kron(d2d3,E3);
G = kron(ir,(1/8/pi)*eye(3)) + (1/8/pi)*dxd_ir3; %  incl prefac


end

function T = stokes_traction(r, rout, nn)
%stokes_traction  Evaluate 3D traction matrix for Stokes single-layer kernel.
%
%   T = stokes_traction(r, rout, nn)
%
%   Constructs the traction matrix mapping source forces at locations `r`
%   to surface tractions (stress dotted with normal) at target points `rout`
%   with outward unit normals `nn`. This corresponds to applying the 
%   Stokeslet traction kernel in 3D free-space.
%
%   INPUTS:
%     r      - N×3 array of source points (e.g., inner surface nodes).
%     rout   - M×3 array of target points (e.g., outer surface nodes).
%     nn     - M×3 array of unit normals at each target point in `rout`.
%
%   OUTPUT:
%     T      - 3M×3N traction matrix, such that T * f approximates the 
%              surface traction field induced by source forces f.
%
%   NOTES:
%     - This routine assumes Stokes flow in 3D, with viscosity μ = 1.
%     - The matrix maps from 3N source force densities to 3M traction values.
%     - The implementation uses the singular fundamental solution and assumes
%       sources and targets are distinct (no regularization).
%
% Anna Broms, Jan 26, 2026

%Note that there is also fmm code to compare to in eg FMM3D (flatiron)

r = r'; 
rout = rout'; 

xsrc = r(1,:); 
ysrc = r(2,:);
zsrc = r(3,:);
xtar = rout(1,:); 
ytar = rout(2,:); 
ztar = rout(3,:); 
nx = nn(:,1);
ny = nn(:,2); 
nz = nn(:,3); 
mu = 1; 

Nsrc = numel(xsrc);
Ntar = numel(xtar);
T = zeros(3*Ntar,3*Nsrc);

prefac = 6 / (8*pi*mu);

for j = 1:Ntar
    %Vector differences r = x - y_j
    rx = xtar(j) - xsrc;
    ry = ytar(j) - ysrc;
    rz = ztar(j) - zsrc;
    r2 = rx.^2 + ry.^2 + rz.^2;
    r5 = sqrt(r2).^5;

    %Dot with normal at target
    rdotn = rx*nx(j) + ry*ny(j) + rz*nz(j);

    %Compute 3x3 block per target
    Txx = prefac * rdotn .* (rx.*rx) ./ r5;
    Txy = prefac * rdotn .* (rx.*ry) ./ r5;
    Txz = prefac * rdotn .* (rx.*rz) ./ r5;
    Tyx = prefac * rdotn .* (ry.*rx) ./ r5;
    Tyy = prefac * rdotn .* (ry.*ry) ./ r5;
    Tyz = prefac * rdotn .* (ry.*rz) ./ r5;
    Tzx = prefac * rdotn .* (rz.*rx) ./ r5;
    Tzy = prefac * rdotn .* (rz.*ry) ./ r5;
    Tzz = prefac * rdotn .* (rz.*rz) ./ r5;

    %Assemble into block structure
    T(3*j-2:3*j,:) = reshape( ...
    [Txx; Txy; Txz;  ...
     Tyx; Tyy; Tyz;  ...
     Tzx; Tzy; Tzz], ...
    3, 3*Nsrc);
end

 
end


function T = stokes_DLP(source, target, dir)
%STOKES_DLP 3D Stokes double-layer potential (stresslet) matrix.
%
%   T = STOKES_DLP(source, target, dir) returns the target-from-source
%   matrix for the free-space 3D Stokes double-layer potential (stresslet
%   kernel). The matrix maps force dipole strengths located at SOURCE points,
%   with dipole directions given by DIR, to velocities evaluated at TARGET
%   points. The viscosity is taken as mu = 1.
%
%   INPUTS:
%     source - Ns-by-3 real array of source locations s_k = (s_xk, s_yk, s_zk)
%     target - Nt-by-3 real array of target locations t_m = (t_xm, t_ym, t_zm)
%     dir    - Ns-by-3 real array of dipole directions d_k at source points
%
%   OUTPUT:
%     T      - (3*Nt)-by-(3*Ns) real matrix mapping dipole strengths
%              f = [f_x1; f_y1; f_z1; f_x2; f_y2; f_z2; ...; f_xNs; f_yNs; f_zNs]
%              to velocities
%              u = [u_x1; u_y1; u_z1; u_x2; u_y2; u_z2; ...; u_xNt; u_yNt; u_zNt]
%              via
%
%                  u = T * f .
%
%              Each 3-by-3 block corresponds to one target/source pair (t_m, s_k).
%              Block ordering is x-, y-, then z-components.
%
%   KERNEL (stresslet / double-layer kernel):
%     Let r = t - s be the separation vector between a target t and a source s,
%     with r = |r|. For a prescribed dipole direction d, the (velocity) stresslet
%     kernel used here is
%
%       T_ij(r, d) = -6/(8*pi*mu) * ( (r·d) r_i r_j ) / r^5 .
%
%     This routine assembles the dense matrix such that
%
%       u_i(t_m) = sum_k T_ij(t_m - s_k, d_k) f_j(s_k).
%
%     (Here mu = 1 is assumed.)
%
%   SINGULARITY / SELF TERMS:
%     The kernel is singular like 1/|r|^2 as r -> 0. No special treatment of
%     diagonal/self interactions is performed; the routine assumes that
%     source and target sets are disjoint.
%
%   NOTES:
%     This matrix is expensive to form and store for large Nt or Ns; it is
%     mainly useful for verification or small problems. Fast evaluation can
%     be achieved using FMM-based implementations.
%
%   See also: STOKES_SLP_3D, STOKES_SLP_2D, STOKES_DIPOLE_2D
%
%   Anna Broms, Feb 3, 2026

r    = source.';
rout = target.';

xsrc = r(1,:);
ysrc = r(2,:);
zsrc = r(3,:);
xtar = rout(1,:);
ytar = rout(2,:);
ztar = rout(3,:);

dx = dir(:,1);
dy = dir(:,2);
dz = dir(:,3);

mu = 1;

Nsrc = numel(xsrc);
Ntar = numel(xtar);
T = zeros(3*Ntar, 3*Nsrc);

prefac = -6 / (8*pi*mu);

for j = 1:Nsrc
    % Differences r = x - y_j (targets minus jth source)
    rx = xtar - xsrc(j);
    ry = ytar - ysrc(j);
    rz = ztar - zsrc(j);

    r2 = rx.^2 + ry.^2 + rz.^2;
    r5 = sqrt(r2).^5;

    % Dot with dipole direction at source j
    rdotd = rx*dx(j) + ry*dy(j) + rz*dz(j);

    % 3x3 block entries (each is Nt-by-1)
    Txx = prefac * rdotd .* (rx.*rx) ./ r5;
    Txy = prefac * rdotd .* (rx.*ry) ./ r5;
    Txz = prefac * rdotd .* (rx.*rz) ./ r5;

    Tyx = prefac * rdotd .* (ry.*rx) ./ r5;
    Tyy = prefac * rdotd .* (ry.*ry) ./ r5;
    Tyz = prefac * rdotd .* (ry.*rz) ./ r5;

    Tzx = prefac * rdotd .* (rz.*rx) ./ r5;
    Tzy = prefac * rdotd .* (rz.*ry) ./ r5;
    Tzz = prefac * rdotd .* (rz.*rz) ./ r5;

    % Assemble into block structure
    T(:, 3*j-2:3*j) = reshape( ...
        [Txx; Txy; Txz; ...
         Tyx; Tyy; Tyz; ...
         Tzx; Tzy; Tzz], ...
        3, 3*Ntar).';
end
end



