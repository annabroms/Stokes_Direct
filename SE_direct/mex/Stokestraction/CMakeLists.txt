add_definitions(-DMATLAB_MEX_FILE)
# add_definitions(-DVERBOSE)

## 0P
# 0P full_ext
matlab_add_mex(
  NAME SE0P_Stokestraction_direct_full_ext_mex
  SRC SE0P_Stokestraction_direct.c SE_Stokestraction_direct_mex.c
)
target_compile_definitions(SE0P_Stokestraction_direct_full_ext_mex PUBLIC ZERO_PERIODIC EWALD_FULL EXTERNAL)
