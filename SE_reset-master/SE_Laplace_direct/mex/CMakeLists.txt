## FGG MEX libraries
add_library(
  mathint
  mathint.c mathint.h SE_direct.h math_x86.h
  )
target_compile_definitions(mathint PUBLIC POTENTIAL gls gslcblas)
add_library(
  mathint_force
  mathint.c mathint.h SE_direct.h math_x86.h
  )
target_compile_definitions(mathint_force PUBLIC FORCE gls gslcblas)

## 3P
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_mex
  SRC SE3P_Laplace_direct_real.c
  )
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_force_mex
  SRC SE3P_Laplace_direct_real.c
  )
target_compile_definitions(SE3P_Laplace_direct_real_force_mex PUBLIC FORCE)
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_rc_mex
  SRC SE3P_Laplace_direct_real.c
  )
target_compile_definitions(SE3P_Laplace_direct_real_rc_mex PUBLIC CUTOFF)
matlab_add_mex(
  NAME SE3P_Laplace_direct_real_rc_force_mex
  SRC SE3P_Laplace_direct_real.c
  )
target_compile_definitions(SE3P_Laplace_direct_real_rc_force_mex PUBLIC CUTOFF FORCE)
matlab_add_mex(
  NAME SE3P_Laplace_direct_fd_mex
  SRC SE3P_Laplace_direct_fd.c
  )
matlab_add_mex(
  NAME SE3P_Laplace_direct_fd_force_mex
  SRC SE3P_Laplace_direct_fd.c
  )
target_compile_definitions(SE3P_Laplace_direct_fd_force_mex PUBLIC FORCE)
matlab_add_mex(
  NAME SE3P_Laplace_direct_self_mex
  SRC SE3P_Laplace_direct_self.c
  )

## 2P
#  matlab_add_mex(
#    NAME SE2P_direct_real_mex
#    SRC SE2P_direct_real.c
#    )
# FIXME: replace rsrc by real_rc
#  matlab_add_mex(
#    NAME SE2P_direct_rsrc_mex
#    SRC SE2P_direct_rsrc.c
#    )
matlab_add_mex(
  NAME SE2P_Laplace_direct_fd_mex
  SRC SE2P_Laplace_direct_fd.c
  )
#  matlab_add_mex(
#    NAME SE2P_direct_self_mex
#    SRC SE2P_direct_self.c
#    )
matlab_add_mex(
  NAME SE2P_Laplace_direct_k0_mex
  SRC SE2P_Laplace_direct_k0.c
  )

## 1P
matlab_add_mex(
  NAME SE1P_Laplace_direct_fd_mex
  SRC SE1P_Laplace_direct_fd.c
  )
target_link_libraries(SE1P_Laplace_direct_fd_mex mathint gsl gslcblas)
matlab_add_mex(
  NAME SE1P_Laplace_direct_fd_force_mex
  SRC SE1P_Laplace_direct_fd.c
  )
target_link_libraries(SE1P_Laplace_direct_fd_force_mex mathint_force gsl gslcblas)
matlab_add_mex(
  NAME SE1P_Laplace_direct_k0_mex
  SRC SE1P_Laplace_direct_k0.c
  )
target_link_libraries(SE1P_Laplace_direct_k0_mex mathint gsl gslcblas)
matlab_add_mex(
  NAME SE1P_Laplace_direct_k0_force_mex
  SRC SE1P_Laplace_direct_k0.c
  )
target_link_libraries(SE1P_Laplace_direct_k0_force_mex mathint_force gsl gslcblas)
