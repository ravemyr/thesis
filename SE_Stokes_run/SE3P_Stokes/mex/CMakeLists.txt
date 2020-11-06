
add_definitions(-DHASIMOTO -DMATLAB_MEX_FILE)
# add_definitions(-DVERBOSE)

matlab_add_mex(
  NAME se3p_fast_k_scaling
  SRC se3p_fast_k_scaling.c
)
