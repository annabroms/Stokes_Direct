function U = SE0P_Stokestraction_direct_full_ext(xs, qs_3D, nnt, opt)
%SE0P_STOKESTRACTION_DIRECT_FULL_EXT Wrapper for Stokeslet traction MEX.
%
%   U = SE0P_Stokestraction_direct_full_ext(xs, qs_3D, nnt, opt)
%
%   This calls SE0P_Stokestraction_direct_full_ext_mex with the same
%   arguments and returns the raw MEX output (num_eval-by-3). Any
%   prefactors (e.g., 1/(8*pi)) should be applied by the caller.

if nargin < 4 || isempty(opt)
    opt = struct();
end

U = SE0P_Stokestraction_direct_full_ext_mex(xs, qs_3D, nnt, opt);
end
