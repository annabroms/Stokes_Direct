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
%     - This routine assumes Stokes flow in 3D, with viscosity mu = 1.
%     - The matrix maps from 3N source force densities to 3M traction values.
%     - The implementation uses the singular fundamental solution and assumes
%       sources and targets are distinct (no regularization).
%
%   Anna Broms, Jan 26, 2026

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
