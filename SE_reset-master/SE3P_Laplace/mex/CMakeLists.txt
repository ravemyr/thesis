if(MKL_FOUND)
  # Build SE3P_Laplace realspace using cell list and MKL
  matlab_add_mex(
    NAME SE3P_Laplace_real_rc_cell_mkl_mex
    SRC SE3P_Laplace_real_rc_cell.c cell_list.c
    )
  # Linking twice, this might have to do with circular dependencies
  target_link_libraries(SE3P_Laplace_real_rc_cell_mkl_mex ${MKL_STATIC})
  target_link_libraries(SE3P_Laplace_real_rc_cell_mkl_mex ${MKL_STATIC})
  target_compile_definitions(SE3P_Laplace_real_rc_cell_mkl_mex PUBLIC INTEL_MKL)
endif()

# Non-MKL version
matlab_add_mex(
  NAME SE3P_Laplace_real_rc_cell_mex
  SRC SE3P_Laplace_real_rc_cell.c cell_list.c
  )
matlab_add_mex(
  NAME SE3P_Laplace_real_rc_cell_force_mex
  SRC SE3P_Laplace_real_rc_cell.c cell_list.c
  )
target_compile_definitions(SE3P_Laplace_real_rc_cell_force_mex PUBLIC FORCE)
