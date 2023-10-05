set(libf2c_SOURCES 
	f77vers.c i77vers.c main.c s_rnge.c abort_.c exit_.c getarg_.c iargc_.c
	getenv_.c signal_.c s_stop.c s_paus.c system_.c cabs.c ctype.c
	derf_.c derfc_.c erf_.c erfc_.c sig_die.c uninit.c
	pow_ci.c pow_dd.c pow_di.c pow_hh.c pow_ii.c pow_ri.c pow_zi.c pow_zz.c
	c_abs.c c_cos.c c_div.c c_exp.c c_log.c c_sin.c c_sqrt.c
	z_abs.c z_cos.c z_div.c z_exp.c z_log.c z_sin.c z_sqrt.c
	r_abs.c r_acos.c r_asin.c r_atan.c r_atn2.c r_cnjg.c r_cos.c
	r_cosh.c r_dim.c r_exp.c r_imag.c r_int.c
	r_lg10.c r_log.c r_mod.c r_nint.c r_sign.c
	r_sin.c r_sinh.c r_sqrt.c r_tan.c r_tanh.c
	d_abs.c d_acos.c d_asin.c d_atan.c d_atn2.c
	d_cnjg.c d_cos.c d_cosh.c d_dim.c d_exp.c
	d_imag.c d_int.c d_lg10.c d_log.c d_mod.c
	d_nint.c d_prod.c d_sign.c d_sin.c d_sinh.c
	d_sqrt.c d_tan.c d_tanh.c
	i_abs.c i_dim.c i_dnnt.c i_indx.c i_len.c i_mod.c i_nint.c i_sign.c
	lbitbits.c lbitshft.c
	h_abs.c h_dim.c h_dnnt.c h_indx.c h_len.c h_mod.c h_nint.c h_sign.c
	l_ge.c l_gt.c l_le.c l_lt.c hl_ge.c hl_gt.c hl_le.c hl_lt.c
	ef1asc_.c ef1cmc_.c
	f77_aloc.c s_cat.c s_cmp.c s_copy.c
	backspac.c close.c dfe.c dolio.c due.c endfile.c err.c
	fmt.c fmtlib.c ftell_.c iio.c ilnw.c inquire.c lread.c lwrite.c
	open.c rdfmt.c rewind.c rsfe.c rsli.c rsne.c sfe.c sue.c
	typesize.c uio.c util.c wref.c wrtfmt.c wsfe.c wsle.c wsne.c xwsne.c
#	pow_qq.c qbitbits.c qbitshft.c ftell64_.c
	dtime_.c etime_.c

)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/signal1.h0 ${CMAKE_CURRENT_BINARY_DIR}/signal1.h @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/sysdep1.h0 ${CMAKE_CURRENT_BINARY_DIR}/sysdep1.h @ONLY)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/f2c.h0"  temp1)
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/f2ch.add"  temp2)
string(APPEND temp1 "${temp2}")
file(WRITE  ${CMAKE_CURRENT_BINARY_DIR}/f2c.h "${temp1}")

add_executable(arithchk arithchk.c)
target_compile_definitions(arithchk PRIVATE NO_FPINIT)
target_link_libraries(arithchk PRIVATE m)
add_custom_command(
      OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/arith.h 
      COMMAND arithchk > arith.h
      DEPENDS arithchk
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)
add_custom_target( arith DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/arith.h)

sherpa_create_git_info(libf2c)
 
add_library(libf2cSherpa SHARED ${libf2c_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/Git_Info.C ${CMAKE_CURRENT_BINARY_DIR}/f2c.h )
target_include_directories(libf2cSherpa PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
add_dependencies(  libf2cSherpa arith )
sherpa_handle_shared_library(libf2cSherpa f2cSherpa)
