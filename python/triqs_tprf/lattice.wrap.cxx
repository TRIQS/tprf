
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>

using c2py::operator""_a;

// ==================== enums =====================

template <> constexpr bool c2py::is_wrapped<triqs_tprf::Channel_t> = true;
template <>
const std::map<triqs_tprf::Channel_t, str_t> c2py::enum_to_string<triqs_tprf::Channel_t> = {{triqs_tprf::Channel_t::PP, "PP"},
                                                                                            {triqs_tprf::Channel_t::PH, "PH"},
                                                                                            {triqs_tprf::Channel_t::PH_bar, "PH_bar"}};

// ==================== module classes =====================

// ==================== module functions ====================

// add_dynamic_and_static
static auto const _c2py_fun_0 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_fk_t g_dyn_fk, triqs_tprf::e_k_t g_stat_k) { return triqs_tprf::add_dynamic_and_static(g_dyn_fk, g_stat_k); },
              "g_dyn_fk", "g_stat_k"),
   c2py::cfun(
      [](triqs_tprf::chi_fk_t chi_dyn_fk, triqs_tprf::chi_k_t chi_stat_k) { return triqs_tprf::add_dynamic_and_static(chi_dyn_fk, chi_stat_k); },
      "chi_dyn_fk", "chi_stat_k"),
   c2py::cfun([](triqs_tprf::g_wk_t g_dyn_wk, triqs_tprf::e_k_t g_stat_k) { return triqs_tprf::add_dynamic_and_static(g_dyn_wk, g_stat_k); },
              "g_dyn_wk", "g_stat_k"),
   c2py::cfun(
      [](triqs_tprf::chi_wk_t chi_dyn_wk, triqs_tprf::chi_k_t chi_stat_k) { return triqs_tprf::add_dynamic_and_static(chi_dyn_wk, chi_stat_k); },
      "chi_dyn_wk", "chi_stat_k"),
   c2py::cfun([](triqs_tprf::g_Dwk_t g_dyn_wk, triqs_tprf::e_k_t g_stat_k) { return triqs_tprf::add_dynamic_and_static(g_dyn_wk, g_stat_k); },
              "g_dyn_wk", "g_stat_k"),
   c2py::cfun(
      [](triqs_tprf::chi_Dwk_t chi_dyn_wk, triqs_tprf::chi_k_t chi_stat_k) { return triqs_tprf::add_dynamic_and_static(chi_dyn_wk, chi_stat_k); },
      "chi_dyn_wk", "chi_stat_k")};

// attatch_tri_vert
static auto const _c2py_fun_1 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](triqs_tprf::chi_nn_cvt L_wn, triqs_tprf::chi_kwnn_cvt chi_kwnn) { return triqs_tprf::attatch_tri_vert(L_wn, chi_kwnn); }, "L_wn", "chi_kwnn")};

// bose
static auto const _c2py_fun_2 = c2py::dispatcher_f_kw_t{c2py::cfun([](double e) { return triqs_tprf::bose(e); }, "e")};

// chi0_Tr_from_g_Tr_PH
static auto const _c2py_fun_3 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_Tr_cvt g_Tr_les, triqs_tprf::g_Tr_cvt g_Tr_gtr) { return triqs_tprf::chi0_Tr_from_g_Tr_PH(g_Tr_les, g_Tr_gtr); },
              "g_Tr_les", "g_Tr_gtr")};

// chi0_nr_from_gr_PH_at_specific_w
static auto const _c2py_fun_4 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](int nw_index, int nn, triqs_tprf::g_wr_cvt g_nr) { return triqs_tprf::chi0_nr_from_gr_PH_at_specific_w(nw_index, nn, g_nr); },
              "nw_index", "nn", "g_nr")};

// chi0_tr_from_grt_PH
static auto const _c2py_fun_5 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_tr_cvt g_tr) { return triqs_tprf::chi0_tr_from_grt_PH(g_tr); }, "g_tr"),
   c2py::cfun([](triqs_tprf::g_tr_cvt g_tr, triqs_tprf::g_tr_cvt g_bwd_tr) { return triqs_tprf::chi0_tr_from_grt_PH(g_tr, g_bwd_tr); }, "g_tr",
              "g_bwd_tr"),
   c2py::cfun([](triqs_tprf::g_Dtr_cvt g_tr) { return triqs_tprf::chi0_tr_from_grt_PH(g_tr); }, "g_tr"),
   c2py::cfun([](triqs_tprf::g_Dtr_cvt g_tr, triqs_tprf::g_Dtr_cvt g_bwd_tr) { return triqs_tprf::chi0_tr_from_grt_PH(g_tr, g_bwd_tr); }, "g_tr",
              "g_bwd_tr")};

// chi0_w0r_from_grt_PH
static auto const _c2py_fun_6 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_tr_cvt g_tr) { return triqs_tprf::chi0_w0r_from_grt_PH(g_tr); }, "g_tr"),
   c2py::cfun([](triqs_tprf::g_tr_cvt g_tr, triqs_tprf::g_tr_cvt g_bwd_tr) { return triqs_tprf::chi0_w0r_from_grt_PH(g_tr, g_bwd_tr); }, "g_tr",
              "g_bwd_tr"),
   c2py::cfun([](triqs_tprf::g_Dtr_cvt g_tr) { return triqs_tprf::chi0_w0r_from_grt_PH(g_tr); }, "g_tr"),
   c2py::cfun([](triqs_tprf::g_Dtr_cvt g_tr, triqs_tprf::g_Dtr_cvt g_bwd_tr) { return triqs_tprf::chi0_w0r_from_grt_PH(g_tr, g_bwd_tr); }, "g_tr",
              "g_bwd_tr")};

// chi0_wr_from_grt_PH
static auto const _c2py_fun_7 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_tr_cvt g_tr, int nw) { return triqs_tprf::chi0_wr_from_grt_PH(g_tr, nw); }, "g_tr", "nw"),
   c2py::cfun([](triqs_tprf::g_tr_cvt g_tr, triqs_tprf::g_tr_cvt g_bwd_tr, int nw) { return triqs_tprf::chi0_wr_from_grt_PH(g_tr, g_bwd_tr, nw); },
              "g_tr", "g_bwd_tr", "nw")};

// chi0q_from_chi0r
static auto const _c2py_fun_8 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_wnr_cvt chi_wnr) { return triqs_tprf::chi0q_from_chi0r(chi_wnr); }, "chi_wnr")};

// chi0q_from_g_wk_PH
static auto const _c2py_fun_9 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](int nw, int nn, triqs_tprf::g_wk_cvt g_wk) { return triqs_tprf::chi0q_from_g_wk_PH(nw, nn, g_wk); }, "nw", "nn", "g_wk")};

// chi0q_sum_nu
static auto const _c2py_fun_10 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_wnk_cvt chi_wnk) { return triqs_tprf::chi0q_sum_nu(chi_wnk); }, "chi_wnk")};

// chi0q_sum_nu_q
static auto const _c2py_fun_11 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_wnk_cvt chi_wnk) { return triqs_tprf::chi0q_sum_nu_q(chi_wnk); }, "chi_wnk")};

// chi0q_sum_nu_tail_corr_PH
static auto const _c2py_fun_12 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_wnk_cvt chi_wnk) { return triqs_tprf::chi0q_sum_nu_tail_corr_PH(chi_wnk); }, "chi_wnk")};

// chi0r_from_chi0q
static auto const _c2py_fun_13 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_wnk_cvt chi_wnk) { return triqs_tprf::chi0r_from_chi0q(chi_wnk); }, "chi_wnk")};

// chi0r_from_gr_PH
static auto const _c2py_fun_14 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](int nw, int nn, triqs_tprf::g_wr_cvt g_nr) { return triqs_tprf::chi0r_from_gr_PH(nw, nn, g_nr); }, "nw", "nn", "g_nr")};

// chi0r_from_gr_PH_nompi
static auto const _c2py_fun_15 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](int nw, int nn, triqs_tprf::g_wr_cvt g_nr) { return triqs_tprf::chi0r_from_gr_PH_nompi(nw, nn, g_nr); }, "nw", "nn", "g_nr")};

// chi_tr_from_chi_wr
static auto const _c2py_fun_16 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wr_cvt chi_wr, int ntau) { return triqs_tprf::chi_tr_from_chi_wr(chi_wr, ntau); }, "chi_wr", "ntau"_a = -1),
   c2py::cfun([](triqs_tprf::chi_Dwr_cvt chi_wr, int ntau) { return triqs_tprf::chi_tr_from_chi_wr(chi_wr, ntau); }, "chi_wr", "ntau"_a = -1)};

// chi_trapz_tau
static auto const _c2py_fun_17 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_t_cvt chi_t) { return triqs_tprf::chi_trapz_tau(chi_t); }, "chi_t")};

// chi_w0r_from_chi_tr
static auto const _c2py_fun_18 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_tr_cvt chi_tr) { return triqs_tprf::chi_w0r_from_chi_tr(chi_tr); }, "chi_tr")};

// chi_wk_from_chi_wr
static auto const _c2py_fun_19 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_wr_cvt chi_wr) { return triqs_tprf::chi_wk_from_chi_wr(chi_wr); }, "chi_wr"),
                           c2py::cfun([](triqs_tprf::chi_Dwr_cvt chi_wr) { return triqs_tprf::chi_wk_from_chi_wr(chi_wr); }, "chi_wr")};

// chi_wr_from_chi_tr
static auto const _c2py_fun_20 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_tr_cvt chi_tr, int nw) { return triqs_tprf::chi_wr_from_chi_tr(chi_tr, nw); }, "chi_tr", "nw"),
   c2py::cfun([](triqs_tprf::chi_Dtr_cvt chi_tr, int nw) { return triqs_tprf::chi_wr_from_chi_tr(chi_tr, nw); }, "chi_tr", "nw")};

// chi_wr_from_chi_wk
static auto const _c2py_fun_21 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_wk_cvt chi_wk) { return triqs_tprf::chi_wr_from_chi_wk(chi_wk); }, "chi_wk"),
                           c2py::cfun([](triqs_tprf::chi_Dwk_cvt chi_wk) { return triqs_tprf::chi_wr_from_chi_wk(chi_wk); }, "chi_wk")};

// chiq_from_chi0q_and_gamma_PH
static auto const _c2py_fun_22 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wnk_cvt chi0_wnk,
                 triqs_tprf::chi_wnn_cvt gamma_ph_wnn) { return triqs_tprf::chiq_from_chi0q_and_gamma_PH(chi0_wnk, gamma_ph_wnn); },
              "chi0_wnk", "gamma_ph_wnn")};

// chiq_sum_nu
static auto const _c2py_fun_23 = c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chiq_t chiq) { return triqs_tprf::chiq_sum_nu(chiq); }, "chiq")};

// chiq_sum_nu_from_chi0q_and_gamma_PH
static auto const _c2py_fun_24 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wnk_cvt chi0_wnk,
                 triqs_tprf::chi_wnn_cvt gamma_ph_wnn) { return triqs_tprf::chiq_sum_nu_from_chi0q_and_gamma_PH(chi0_wnk, gamma_ph_wnn); },
              "chi0_wnk", "gamma_ph_wnn")};

// chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH
static auto const _c2py_fun_25 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wnk_cvt chi0_wnk, triqs_tprf::chi_wnn_cvt gamma_ph_wnn,
                 triqs_tprf::chi_nn_cvt L_wn) { return triqs_tprf::chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH(chi0_wnk, gamma_ph_wnn, L_wn); },
              "chi0_wnk", "gamma_ph_wnn", "L_wn")};

// chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH
static auto const _c2py_fun_26 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](double mu, triqs_tprf::ek_vt e_k, triqs_tprf::g_iw_vt sigma_w, triqs_tprf::g2_iw_vt gamma_ph_wnn,
                 int tail_corr_nwf) { return triqs_tprf::chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH(mu, e_k, sigma_w, gamma_ph_wnn, tail_corr_nwf); },
              "mu", "e_k", "sigma_w", "gamma_ph_wnn", "tail_corr_nwf"_a = -1)};

// chiq_sum_nu_from_g_wk_and_gamma_PH
static auto const _c2py_fun_27 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::gk_iw_t g_wk, triqs_tprf::g2_iw_vt gamma_ph_wnn,
                 int tail_corr_nwf) { return triqs_tprf::chiq_sum_nu_from_g_wk_and_gamma_PH(g_wk, gamma_ph_wnn, tail_corr_nwf); },
              "g_wk", "gamma_ph_wnn", "tail_corr_nwf"_a = -1)};

// chiq_sum_nu_q
static auto const _c2py_fun_28 = c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chiq_t chiq) { return triqs_tprf::chiq_sum_nu_q(chiq); }, "chiq")};

// cluster_mesh_fourier_interpolation
static auto const _c2py_fun_29 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](nda::basic_array<double, 2, nda::C_layout, 'A', nda::heap_basic<nda::mem::mallocator<nda::mem::AddressSpace::Host>>> k_vecs,
                 triqs_tprf::chi_wr_cvt chi) { return triqs_tprf::cluster_mesh_fourier_interpolation(k_vecs, chi); },
              "k_vecs", "chi")};

// construct_phi_wk
static auto const _c2py_fun_30 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](triqs_tprf::chi_wk_vt chi,
      nda::basic_array_view<std::complex<double>, 4, nda::C_layout, 'A', nda::default_accessor, nda::borrowed<nda::mem::AddressSpace::Host>> U) {
     return triqs_tprf::construct_phi_wk(chi, U);
   },
   "chi", "U")};

// dlr_on_imfreq
static auto const _c2py_fun_31 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_Dc_cvt g_c, triqs::mesh::imfreq wmesh) { return triqs_tprf::dlr_on_imfreq(g_c, wmesh); }, "g_c", "wmesh"),
   c2py::cfun([](triqs_tprf::chi_Dc_cvt chi_c, triqs::mesh::imfreq wmesh) { return triqs_tprf::dlr_on_imfreq(chi_c, wmesh); }, "chi_c", "wmesh")};

// dlr_on_imtime
static auto const _c2py_fun_32 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_Dc_cvt g_c, triqs::mesh::imtime tmesh) { return triqs_tprf::dlr_on_imtime(g_c, tmesh); }, "g_c", "tmesh"),
   c2py::cfun([](triqs_tprf::chi_Dc_cvt chi_c, triqs::mesh::imtime tmesh) { return triqs_tprf::dlr_on_imtime(chi_c, tmesh); }, "chi_c", "tmesh")};

// dynamic_and_constant_to_tr
static auto const _c2py_fun_33 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wk_vt Gamma_pp_dyn_wk,
                 triqs_tprf::chi_k_vt Gamma_pp_const_k) { return triqs_tprf::dynamic_and_constant_to_tr(Gamma_pp_dyn_wk, Gamma_pp_const_k); },
              "Gamma_pp_dyn_wk", "Gamma_pp_const_k"),
   c2py::cfun([](triqs_tprf::chi_Dwk_vt Gamma_pp_dyn_wk,
                 triqs_tprf::chi_k_vt Gamma_pp_const_k) { return triqs_tprf::dynamic_and_constant_to_tr(Gamma_pp_dyn_wk, Gamma_pp_const_k); },
              "Gamma_pp_dyn_wk", "Gamma_pp_const_k")};

// dynamical_screened_interaction_W
static auto const _c2py_fun_34 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wk_cvt PI_wk, triqs_tprf::chi_k_cvt V_k) { return triqs_tprf::dynamical_screened_interaction_W(PI_wk, V_k); },
              "PI_wk", "V_k"),
   c2py::cfun([](triqs_tprf::chi_Dwk_cvt PI_wk, triqs_tprf::chi_k_cvt V_k) { return triqs_tprf::dynamical_screened_interaction_W(PI_wk, V_k); },
              "PI_wk", "V_k"),
   c2py::cfun([](triqs_tprf::chi_fk_cvt PI_fk, triqs_tprf::chi_k_cvt V_k) { return triqs_tprf::dynamical_screened_interaction_W(PI_fk, V_k); },
              "PI_fk", "V_k"),
   c2py::cfun([](triqs_tprf::chi_wk_cvt PI_wk, triqs_tprf::chi_wk_cvt V_wk) { return triqs_tprf::dynamical_screened_interaction_W(PI_wk, V_wk); },
              "PI_wk", "V_wk"),
   c2py::cfun([](triqs_tprf::chi_Dwk_cvt PI_wk, triqs_tprf::chi_Dwk_cvt V_wk) { return triqs_tprf::dynamical_screened_interaction_W(PI_wk, V_wk); },
              "PI_wk", "V_wk"),
   c2py::cfun([](triqs_tprf::chi_fk_cvt PI_fk, triqs_tprf::chi_fk_cvt V_fk) { return triqs_tprf::dynamical_screened_interaction_W(PI_fk, V_fk); },
              "PI_fk", "V_fk")};

// dynamical_screened_interaction_W_from_generalized_susceptibility
static auto const _c2py_fun_35 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wk_cvt chi_wk,
                 triqs_tprf::chi_k_cvt V_k) { return triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk, V_k); },
              "chi_wk", "V_k"),
   c2py::cfun([](triqs_tprf::chi_Dwk_cvt chi_wk,
                 triqs_tprf::chi_k_cvt V_k) { return triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk, V_k); },
              "chi_wk", "V_k"),
   c2py::cfun([](triqs_tprf::chi_fk_cvt chi_fk,
                 triqs_tprf::chi_k_cvt V_k) { return triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility(chi_fk, V_k); },
              "chi_fk", "V_k"),
   c2py::cfun([](triqs_tprf::chi_wk_cvt chi_wk,
                 triqs_tprf::chi_wk_cvt V_wk) { return triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk, V_wk); },
              "chi_wk", "V_wk"),
   c2py::cfun([](triqs_tprf::chi_Dwk_cvt chi_wk,
                 triqs_tprf::chi_Dwk_cvt V_wk) { return triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility(chi_wk, V_wk); },
              "chi_wk", "V_wk"),
   c2py::cfun([](triqs_tprf::chi_fk_cvt chi_fk,
                 triqs_tprf::chi_fk_cvt V_fk) { return triqs_tprf::dynamical_screened_interaction_W_from_generalized_susceptibility(chi_fk, V_fk); },
              "chi_fk", "V_fk")};

// eliashberg_constant_gamma_f_product
static auto const _c2py_fun_36 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_r_vt Gamma_pp_const_r,
                                         triqs_tprf::g_tr_t F_tr) { return triqs_tprf::eliashberg_constant_gamma_f_product(Gamma_pp_const_r, F_tr); },
                                      "Gamma_pp_const_r", "F_tr")};

// eliashberg_g_delta_g_product
static auto const _c2py_fun_37 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_wk_vt g_wk, triqs_tprf::g_wk_vt delta_wk) { return triqs_tprf::eliashberg_g_delta_g_product(g_wk, delta_wk); }, "g_wk",
              "delta_wk"),
   c2py::cfun([](triqs_tprf::g_Dwk_vt g_wk, triqs_tprf::g_Dwk_vt delta_wk) { return triqs_tprf::eliashberg_g_delta_g_product(g_wk, delta_wk); },
              "g_wk", "delta_wk")};

// eliashberg_product
static auto const _c2py_fun_38 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::chi_wk_vt Gamma_pp, triqs_tprf::g_wk_vt g_wk,
                                         triqs_tprf::g_wk_vt delta_wk) { return triqs_tprf::eliashberg_product(Gamma_pp, g_wk, delta_wk); },
                                      "Gamma_pp", "g_wk", "delta_wk")};

// eliashberg_product_fft
static auto const _c2py_fun_39 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_tr_vt Gamma_pp_dyn_tr, triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_wk_vt g_wk,
                 triqs_tprf::g_wk_vt delta_wk) { return triqs_tprf::eliashberg_product_fft(Gamma_pp_dyn_tr, Gamma_pp_const_r, g_wk, delta_wk); },
              "Gamma_pp_dyn_tr", "Gamma_pp_const_r", "g_wk", "delta_wk"),
   c2py::cfun([](triqs_tprf::chi_Dtr_vt Gamma_pp_dyn_tr, triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_Dwk_vt g_wk,
                 triqs_tprf::g_Dwk_vt delta_wk) { return triqs_tprf::eliashberg_product_fft(Gamma_pp_dyn_tr, Gamma_pp_const_r, g_wk, delta_wk); },
              "Gamma_pp_dyn_tr", "Gamma_pp_const_r", "g_wk", "delta_wk")};

// eliashberg_product_fft_constant
static auto const _c2py_fun_40 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_wk_vt g_wk,
                 triqs_tprf::g_wk_vt delta_wk) { return triqs_tprf::eliashberg_product_fft_constant(Gamma_pp_const_r, g_wk, delta_wk); },
              "Gamma_pp_const_r", "g_wk", "delta_wk"),
   c2py::cfun([](triqs_tprf::chi_r_vt Gamma_pp_const_r, triqs_tprf::g_Dwk_vt g_wk,
                 triqs_tprf::g_Dwk_vt delta_wk) { return triqs_tprf::eliashberg_product_fft_constant(Gamma_pp_const_r, g_wk, delta_wk); },
              "Gamma_pp_const_r", "g_wk", "delta_wk")};

// fermi
static auto const _c2py_fun_41 = c2py::dispatcher_f_kw_t{c2py::cfun([](double e) { return triqs_tprf::fermi(e); }, "e")};

// fock_sigma
static auto const _c2py_fun_42 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_wk_cvt g_wk) { return triqs_tprf::fock_sigma(v_k, g_wk); }, "v_k", "g_wk"),
   c2py::cfun([](triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_Dwk_cvt g_wk) { return triqs_tprf::fock_sigma(v_k, g_wk); }, "v_k", "g_wk"),
   c2py::cfun([](triqs_tprf::chi_r_cvt v_r, triqs_tprf::e_r_cvt rho_r) { return triqs_tprf::fock_sigma(v_r, rho_r); }, "v_r", "rho_r")};

// fourier_Tk_to_Tr
static auto const _c2py_fun_43 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_Tk_cvt g_Tk) { return triqs_tprf::fourier_Tk_to_Tr(g_Tk); }, "g_Tk")};

// fourier_Tr_to_Tk
static auto const _c2py_fun_44 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_Tr_cvt g_Tr) { return triqs_tprf::fourier_Tr_to_Tk(g_Tr); }, "g_Tr"),
                           c2py::cfun([](triqs_tprf::chi_Tr_cvt chi_Tr) { return triqs_tprf::fourier_Tr_to_Tk(chi_Tr); }, "chi_Tr")};

// fourier_fk_to_fr
static auto const _c2py_fun_45 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_fk_cvt g_fk) { return triqs_tprf::fourier_fk_to_fr(g_fk); }, "g_fk")};

// fourier_fr_to_fk
static auto const _c2py_fun_46 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_fr_cvt g_fr) { return triqs_tprf::fourier_fr_to_fk(g_fr); }, "g_fr")};

// fourier_tr_to_wr
static auto const _c2py_fun_47 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_tr_cvt g_tr, int nw) { return triqs_tprf::fourier_tr_to_wr(g_tr, nw); }, "g_tr", "nw"_a = -1),
   c2py::cfun([](triqs_tprf::g_Dtr_cvt g_tr, int nw) { return triqs_tprf::fourier_tr_to_wr(g_tr, nw); }, "g_tr", "nw"_a = -1),
   c2py::cfun([](triqs_tprf::chi_Dtr_cvt chi_tr, int nw) { return triqs_tprf::fourier_tr_to_wr(chi_tr, nw); }, "chi_tr", "nw"_a = -1)};

// fourier_wk_to_wr
static auto const _c2py_fun_48 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_wk_cvt g_wk) { return triqs_tprf::fourier_wk_to_wr(g_wk); }, "g_wk"),
                           c2py::cfun([](triqs_tprf::g_Dwk_cvt g_wk) { return triqs_tprf::fourier_wk_to_wr(g_wk); }, "g_wk"),
                           c2py::cfun([](triqs_tprf::chi_Dwk_cvt chi_wk) { return triqs_tprf::fourier_wk_to_wr(chi_wk); }, "chi_wk")};

// fourier_wr_to_tr
static auto const _c2py_fun_49 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::g_wr_cvt g_wr, int nt) { return triqs_tprf::fourier_wr_to_tr(g_wr, nt); }, "g_wr", "nt"_a = -1),
   c2py::cfun([](triqs_tprf::g_Dwr_cvt g_wr, int nt) { return triqs_tprf::fourier_wr_to_tr(g_wr, nt); }, "g_wr", "nt"_a = -1),
   c2py::cfun([](triqs_tprf::chi_Dwr_cvt chi_wr, int nt) { return triqs_tprf::fourier_wr_to_tr(chi_wr, nt); }, "chi_wr", "nt"_a = -1)};

// fourier_wr_to_wk
static auto const _c2py_fun_50 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_wr_cvt g_wr) { return triqs_tprf::fourier_wr_to_wk(g_wr); }, "g_wr"),
                           c2py::cfun([](triqs_tprf::g_Dwr_cvt g_wr) { return triqs_tprf::fourier_wr_to_wk(g_wr); }, "g_wr"),
                           c2py::cfun([](triqs_tprf::chi_Dwr_cvt chi_wr) { return triqs_tprf::fourier_wr_to_wk(chi_wr); }, "chi_wr")};

// g0_Tk_les_gtr_from_e_k
static auto const _c2py_fun_51 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::e_k_cvt e_k, triqs::mesh::retime Tmesh, double beta) { return triqs_tprf::g0_Tk_les_gtr_from_e_k(e_k, Tmesh, beta); },
              "e_k", "Tmesh", "beta")};

// g0w_dynamic_sigma
static auto const _c2py_fun_52 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta,
                 triqs::mesh::brzone::value_t kpoint) { return triqs_tprf::g0w_dynamic_sigma(mu, beta, e_k, W_fk, v_k, delta, kpoint); },
              "mu", "beta", "e_k", "W_fk", "v_k", "delta", "kpoint"),
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta,
                 triqs::mesh::brzone kmesh) { return triqs_tprf::g0w_dynamic_sigma(mu, beta, e_k, W_fk, v_k, delta, kmesh); },
              "mu", "beta", "e_k", "W_fk", "v_k", "delta", "kmesh"),
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k,
                 double delta) { return triqs_tprf::g0w_dynamic_sigma(mu, beta, e_k, W_fk, v_k, delta); },
              "mu", "beta", "e_k", "W_fk", "v_k", "delta")};

// g0w_sigma
static auto const _c2py_fun_53 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_k_cvt v_k,
                 triqs::mesh::brzone::value_t kpoint) { return triqs_tprf::g0w_sigma(mu, beta, e_k, v_k, kpoint); },
              "mu", "beta", "e_k", "v_k", "kpoint"),
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_k_cvt v_k,
                 triqs::mesh::brzone kmesh) { return triqs_tprf::g0w_sigma(mu, beta, e_k, v_k, kmesh); },
              "mu", "beta", "e_k", "v_k", "kmesh"),
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_k_cvt v_k) { return triqs_tprf::g0w_sigma(mu, beta, e_k, v_k); },
              "mu", "beta", "e_k", "v_k"),
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta,
                 triqs::mesh::brzone::value_t kpoint) { return triqs_tprf::g0w_sigma(mu, beta, e_k, W_fk, v_k, delta, kpoint); },
              "mu", "beta", "e_k", "W_fk", "v_k", "delta", "kpoint"),
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k, double delta,
                 triqs::mesh::brzone kmesh) { return triqs_tprf::g0w_sigma(mu, beta, e_k, W_fk, v_k, delta, kmesh); },
              "mu", "beta", "e_k", "W_fk", "v_k", "delta", "kmesh"),
   c2py::cfun([](double mu, double beta, triqs_tprf::e_k_cvt e_k, triqs_tprf::chi_fk_cvt W_fk, triqs_tprf::chi_k_cvt v_k,
                 double delta) { return triqs_tprf::g0w_sigma(mu, beta, e_k, W_fk, v_k, delta); },
              "mu", "beta", "e_k", "W_fk", "v_k", "delta")};

// gw_dynamic_sigma
static auto const _c2py_fun_54 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_tr_cvt W_tr, triqs_tprf::g_tr_cvt g_tr) { return triqs_tprf::gw_dynamic_sigma(W_tr, g_tr); }, "W_tr", "g_tr"),
   c2py::cfun([](triqs_tprf::chi_Dtr_cvt W_tr, triqs_tprf::g_Dtr_cvt g_tr) { return triqs_tprf::gw_dynamic_sigma(W_tr, g_tr); }, "W_tr", "g_tr")};

// gw_sigma
static auto const _c2py_fun_55 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wk_cvt W_wk, triqs_tprf::g_wk_cvt g_wk) { return triqs_tprf::gw_sigma(W_wk, g_wk); }, "W_wk", "g_wk"),
   c2py::cfun(
      [](triqs_tprf::chi_Dwk_cvt W_wk, triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_Dwk_cvt g_wk) { return triqs_tprf::gw_sigma(W_wk, v_k, g_wk); },
      "W_wk", "v_k", "g_wk"),
   c2py::cfun([](triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_wk_cvt g_wk) { return triqs_tprf::gw_sigma(v_k, g_wk); }, "v_k", "g_wk"),
   c2py::cfun([](triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_Dwk_cvt g_wk) { return triqs_tprf::gw_sigma(v_k, g_wk); }, "v_k", "g_wk")};

// hartree_sigma
static auto const _c2py_fun_56 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_k_cvt v_k, triqs_tprf::g_wk_cvt g_wk) { return triqs_tprf::hartree_sigma(v_k, g_wk); }, "v_k", "g_wk"),
   c2py::cfun([](triqs_tprf::chi_k_cvt v_k, triqs_tprf::e_r_cvt rho_r) { return triqs_tprf::hartree_sigma(v_k, rho_r); }, "v_k", "rho_r")};

// lattice_dyson_g0_fk
static auto const _c2py_fun_57 = c2py::dispatcher_f_kw_t{c2py::cfun(
   [](double mu, triqs_tprf::e_k_cvt e_k, triqs::mesh::refreq mesh, double delta) { return triqs_tprf::lattice_dyson_g0_fk(mu, e_k, mesh, delta); },
   "mu", "e_k", "mesh", "delta")};

// lattice_dyson_g0_wk
static auto const _c2py_fun_58 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs::mesh::imfreq mesh) { return triqs_tprf::lattice_dyson_g0_wk(mu, e_k, mesh); }, "mu",
              "e_k", "mesh"),
   c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs::mesh::dlr_imfreq mesh) { return triqs_tprf::lattice_dyson_g0_wk(mu, e_k, mesh); }, "mu",
              "e_k", "mesh")};

// lattice_dyson_g_f
static auto const _c2py_fun_59 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_f_cvt sigma_f,
                                         double delta) { return triqs_tprf::lattice_dyson_g_f(mu, e_k, sigma_f, delta); },
                                      "mu", "e_k", "sigma_f", "delta")};

// lattice_dyson_g_fk
static auto const _c2py_fun_60 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_fk_cvt sigma_fk,
                                         double delta) { return triqs_tprf::lattice_dyson_g_fk(mu, e_k, sigma_fk, delta); },
                                      "mu", "e_k", "sigma_fk", "delta"),
                           c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_f_cvt sigma_f,
                                         double delta) { return triqs_tprf::lattice_dyson_g_fk(mu, e_k, sigma_f, delta); },
                                      "mu", "e_k", "sigma_f", "delta")};

// lattice_dyson_g_w
static auto const _c2py_fun_61 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_w_cvt sigma_w) { return triqs_tprf::lattice_dyson_g_w(mu, e_k, sigma_w); }, "mu",
              "e_k", "sigma_w"),
   c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_Dw_cvt sigma_w) { return triqs_tprf::lattice_dyson_g_w(mu, e_k, sigma_w); }, "mu",
              "e_k", "sigma_w")};

// lattice_dyson_g_wk
static auto const _c2py_fun_62 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_w_cvt sigma_w) { return triqs_tprf::lattice_dyson_g_wk(mu, e_k, sigma_w); }, "mu",
              "e_k", "sigma_w"),
   c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_wk_cvt sigma_wk) { return triqs_tprf::lattice_dyson_g_wk(mu, e_k, sigma_wk); },
              "mu", "e_k", "sigma_wk"),
   c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_Dwk_cvt sigma_wk) { return triqs_tprf::lattice_dyson_g_wk(mu, e_k, sigma_wk); },
              "mu", "e_k", "sigma_wk"),
   c2py::cfun([](double mu, triqs_tprf::e_k_cvt e_k, triqs_tprf::g_Dw_cvt sigma_w) { return triqs_tprf::lattice_dyson_g_wk(mu, e_k, sigma_w); }, "mu",
              "e_k", "sigma_w")};

// lindhard_chi00
static auto const _c2py_fun_63 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::e_k_cvt e_k, triqs::mesh::imfreq mesh, double mu) { return triqs_tprf::lindhard_chi00(e_k, mesh, mu); }, "e_k", "mesh",
              "mu"),
   c2py::cfun([](triqs_tprf::e_k_cvt e_k, triqs::mesh::dlr_imfreq mesh, double mu) { return triqs_tprf::lindhard_chi00(e_k, mesh, mu); }, "e_k",
              "mesh", "mu"),
   c2py::cfun([](triqs_tprf::e_k_cvt e_k, triqs::mesh::refreq mesh, double beta, double mu,
                 double delta) { return triqs_tprf::lindhard_chi00(e_k, mesh, beta, mu, delta); },
              "e_k", "mesh", "beta", "mu", "delta")};

// rho_k_from_g_wk
static auto const _c2py_fun_64 =
   c2py::dispatcher_f_kw_t{c2py::cfun([](triqs_tprf::g_wk_cvt g_wk) { return triqs_tprf::rho_k_from_g_wk(g_wk); }, "g_wk"),
                           c2py::cfun([](triqs_tprf::g_Dwk_cvt g_wk) { return triqs_tprf::rho_k_from_g_wk(g_wk); }, "g_wk")};

// solve_rpa_PH
static auto const _c2py_fun_65 = c2py::dispatcher_f_kw_t{
   c2py::cfun(
      [](triqs_tprf::chi_wk_vt chi0,
         nda::basic_array_view<std::complex<double>, 4, nda::C_layout, 'A', nda::default_accessor, nda::borrowed<nda::mem::AddressSpace::Host>> U) {
        return triqs_tprf::solve_rpa_PH(chi0, U);
      },
      "chi0", "U"),
   c2py::cfun(
      [](triqs_tprf::chi_Dwk_vt chi0,
         nda::basic_array_view<std::complex<double>, 4, nda::C_layout, 'A', nda::default_accessor, nda::borrowed<nda::mem::AddressSpace::Host>> U) {
        return triqs_tprf::solve_rpa_PH(chi0, U);
      },
      "chi0", "U"),
   c2py::cfun(
      [](triqs_tprf::chi_fk_vt chi0,
         nda::basic_array_view<std::complex<double>, 4, nda::C_layout, 'A', nda::default_accessor, nda::borrowed<nda::mem::AddressSpace::Host>> U) {
        return triqs_tprf::solve_rpa_PH(chi0, U);
      },
      "chi0", "U")};

// split_into_dynamic_wk_and_constant_k
static auto const _c2py_fun_66 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](triqs_tprf::chi_wk_cvt chi_wk) { return triqs_tprf::split_into_dynamic_wk_and_constant_k(chi_wk); }, "chi_wk")};

static const auto _c2py_doc_0 = _c2py_fun_0.doc(
   R"DOC(
Adds a dynamic and a static Green's function.

Parameters
----------
g_dyn_fk : {par_0}
   : general Green's function :math:`G_{ab}(\omega, \mathbf{k})`.
g_stat_k : {par_1}
   : general Green's function :math:`G_{ab}(\mathbf{k})`.
g_dyn_wk : {par_2}
   : general Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`.

Returns
-------
[1] : {ret_0}
   g_fk : Green's function :math:`G_{ab}(\omega, \mathbf{k}) + G_{ab}(\mathbf{k})`.

[3] : {ret_1}
   g_wk : Green's function :math:`G_{ab}(i\omega_n, \mathbf{k}) + G_{ab}(\mathbf{k})`.
)DOC",
   {{c2py::python_typename<triqs_tprf::g_fk_t>()}, {c2py::python_typename<triqs_tprf::e_k_t>()}, {c2py::python_typename<triqs_tprf::g_wk_t>()}},
   {c2py::python_typename<triqs_tprf::g_fk_t>(), c2py::python_typename<triqs_tprf::g_wk_t>()});
static const auto _c2py_doc_1 = _c2py_fun_1.doc(R"DOC()DOC");
static const auto _c2py_doc_2 = _c2py_fun_2.doc(R"DOC(
Helper function to evaluate the Bose-Einstein distribution function

 .. math ::
     n_B() = {1}{() - 1}

Parameters
----------
e : {par_0}
   : point at which to evaluate :math:`n_B(\epsilon)`.

Returns
-------
{ret_0}
   The value of :math:`n_B(\epsilon)`.
)DOC",
                                                {{c2py::python_typename<double>()}}, {c2py::python_typename<double>()});
static const auto _c2py_doc_3 = _c2py_fun_3.doc(R"DOC(
Generalized susceptibility real time bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(t, \mathbf{r})`

 Computes

 

.. math::

   \chi^{(0)R}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
   i G^<_{d\bar{a}}(t, \mathbf{r}) G^>_{b\bar{c}}(-t, -\mathbf{r})
   - i G^>_{d\bar{a}}(t, \mathbf{r}) G^<_{b\bar{c}}(-t, -\mathbf{r})

Parameters
----------
g_Tr_les : {par_0}
   Lesser real time Green's function in real-space, :math:`G^<_{a\bar{b}}(t, \mathbf{r})`.
g_Tr_gtr : {par_1}
   Greater real time Green's function in real-space, :math:`G^>_{a\bar{b}}(t, \mathbf{r})`.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(t, \mathbf{r})` in real time and real-space.
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g_Tr_cvt>()}, {c2py::python_typename<triqs_tprf::g_Tr_cvt>()}},
                                                {c2py::python_typename<triqs_tprf::chi_Tr_t>()});
static const auto _c2py_doc_4 = _c2py_fun_4.doc(R"DOC()DOC");
static const auto _c2py_doc_5 = _c2py_fun_5.doc(R"DOC(
Generalized susceptibility imaginary time bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
   - G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})

Parameters
----------
g_tr : {par_0}
   Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\tau, \mathbf{r})`.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` in imaginary time and real-space.
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g_tr_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_tr_t>()});
static const auto _c2py_doc_6 = _c2py_fun_6.doc(R"DOC(
Generalized susceptibility zero imaginary frequency bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})`

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{r}) =
   - \int_0^\beta d\tau \,
   G_{d\bar{a}}(\tau, \mathbf{r}) G_{b\bar{c}}(-\tau, -\mathbf{r})

Parameters
----------
g_tr : {par_0}
   Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\tau, \mathbf{r})`.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{r})` in real-space.
)DOC",
                                                {{c2py::python_typename<triqs_tprf::g_tr_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_wr_t>()});
static const auto _c2py_doc_7 = _c2py_fun_7.doc(R"DOC()DOC");
static const auto _c2py_doc_8 =
   _c2py_fun_8.doc(R"DOC(
Fourier transform of the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in real space to :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in momentum space.

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q}) =
   \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{q}} \left\{
   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})
   \right\}

Parameters
----------
chi_wnr : {par_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.
)DOC",
                   {{c2py::python_typename<triqs_tprf::chi_wnr_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_wnk_t>()});
static const auto _c2py_doc_9 =
   _c2py_fun_9.doc(R"DOC(
Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` with convolution in k-space.

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q}) =
   - \frac{\beta}{N_k} \sum_\mathbf{k}
   G_{d\bar{a}}(\nu, \mathbf{k}) \cdot G_{b\bar{c}}(\nu + \omega, \mathbf{k} - \mathbf{q})

Parameters
----------
nw : {par_0}
   Number of bosonic Matsubara freqiencies.
nn : {par_1}
   Number of fermionic Matsubara freqiencies.
g_tr : {par_2}
   Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.
)DOC",
                   {{c2py::python_typename<int>()}, {c2py::python_typename<int>()}, {}}, {c2py::python_typename<triqs_tprf::chi_wnk_t>()});
static const auto _c2py_doc_10 =
   _c2py_fun_10.doc(R"DOC(
Sum over fermionic frequency in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})`. (NB! without tail corrections)

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
   \frac{1}{\beta^2} \sum_{\nu=\nu_{min}}^\nu_{max} \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

Parameters
----------
chi_wnk : {par_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Returns
-------
{ret_0}
   Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic Matsubara frequency and momentum space.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_wnk_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_wk_t>()});
static const auto _c2py_doc_11 =
   _c2py_fun_11.doc(R"DOC(
Sum over fermionic frequency and momentum in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})`. (NB! without tail corrections)

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
   \frac{1}{N_k} \sum_\matbf{k} \frac{1}{\beta^2} \sum_{\nu=\nu_{min}}^\nu_{max}
   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

Parameters
----------
chi_wnk : {par_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Returns
-------
{ret_0}
   Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega)` in one bosonic Matsubara frequency.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_wnk_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_w_t>()});
static const auto _c2py_doc_12 =
   _c2py_fun_12.doc(R"DOC(
Sum over fermionic frequency in the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})` using higher order tail corrections when summing to infinity.

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
   \frac{1}{\beta^2} \sum_{\nu=-\infty}^\infty \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{k})

Parameters
----------
chi_wnk : {par_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Returns
-------
{ret_0}
   Susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` in one bosonic Matsubara frequency and momentum space.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_wnk_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_wk_t>()});
static const auto _c2py_doc_13 =
   _c2py_fun_13.doc(R"DOC(
Fourier transform of the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in momentum-space to :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in real-space.

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
   \mathcal{F}_{\mathbf{q} \rightarrow \mathbf{r}} \left\{
   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})
   \right\}

Parameters
----------
chi_wnk : {par_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{q})` in one bosonic and one fermionic Matsuabara frequency and momentum space.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real space.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_wnk_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_wnr_t>()});
static const auto _c2py_doc_14 =
   _c2py_fun_14.doc(R"DOC(
Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})`.

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
   - \beta G_{d\bar{a}}(\nu, \mathbf{r}) \cdot G_{b\bar{c}}(\nu + \omega, -\mathbf{r})

Parameters
----------
nw : {par_0}
   Number of bosonic Matsubara freqiencies.
nn : {par_1}
   Number of fermionic Matsubara freqiencies.
g_tr : {par_2}
   Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.
)DOC",
                    {{c2py::python_typename<int>()}, {c2py::python_typename<int>()}, {}}, {c2py::python_typename<triqs_tprf::chi_wnr_t>()});
static const auto _c2py_doc_15 =
   _c2py_fun_15.doc(R"DOC(
Generalized susceptibility bubble in the particle-hole channel :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` without MPI parallellization.

 Computes

 

.. math::

   \chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r}) =
   - \beta G_{d\bar{a}}(\nu, \mathbf{r}) \cdot G_{b\bar{c}}(\nu + \omega, -\mathbf{r})

Parameters
----------
nw : {par_0}
   Number of bosonic Matsubara freqiencies.
nn : {par_1}
   Number of fermionic Matsubara freqiencies.
g_tr : {par_2}
   Imaginary time Green's function in real-space, :math:`G_{a\bar{b}}(\nu, \mathbf{r})`.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \mathbf{r})` in one bosonic and one fermionic Matsuabara frequency and real-space.
)DOC",
                    {{c2py::python_typename<int>()}, {c2py::python_typename<int>()}, {}}, {c2py::python_typename<triqs_tprf::chi_wnr_t>()});
static const auto _c2py_doc_16 = _c2py_fun_16.doc(R"DOC(
Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})`

 Computes

 

.. math::

   \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r}) =
   \mathcal{F}_{\omega \rightarrow \tau} \left\{
   \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
   \right\}

Parameters
----------
chi_tr : {par_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` 
                  in imaginary time and real space.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` 
            in Matsubara frequency and real-space.
)DOC",
                                                  {{}}, {c2py::python_typename<triqs_tprf::chi_tr_t>()});
static const auto _c2py_doc_17 = _c2py_fun_17.doc(R"DOC()DOC");
static const auto _c2py_doc_18 =
   _c2py_fun_18.doc(R"DOC(
Static susceptibility calculation :math:`\chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})`
  
 Explicit calculation of the static, zero frequency response, by 2nd order trapetzoidal 
 integration in imaginary time, i.e.

 

.. math::

   \chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r}) =
   \int_0^\beta d\tau \, \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})

Parameters
----------
chi_tr : {par_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` 
                  in imaginary time and real space.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega=0, \mathbf{r})` 
            at zero Matsubara frequency and real-space.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_tr_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_wr_t>()});
static const auto _c2py_doc_19 =
   _c2py_fun_19.doc(R"DOC(
Parallel Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`

 Computes

 

.. math::

   \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
   \mathcal{F}_{\mathbf{r} \rightarrow \mathbf{k}} \left\{
   \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})
   \right\}

Parameters
----------
chi_wr : {par_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` 
                  in Matsubara frequency and real space.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` 
            in Matsubara frequency and momentum space.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_wr_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_wk_t>()});
static const auto _c2py_doc_20 =
   _c2py_fun_20.doc(R"DOC(
Parallel Fourier transform from  :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`

 Computes

 

.. math::

   \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
   \mathcal{F}_{\tau \rightarrow \omega} \left\{
   \chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})
   \right\}

Parameters
----------
chi_tr : {par_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\tau, \mathbf{r})` 
                  in imaginary time and real space.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` 
            in Matsubara frequency and real-space.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_tr_cvt>()}}, {c2py::python_typename<triqs_tprf::chi_wr_t>()});
static const auto _c2py_doc_21 = _c2py_fun_21.doc(R"DOC(
Parallel Fourier transform from :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` to :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})`

 Computes

 

.. math::

   \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r}) =
   \mathcal{F}_{\mathbf{k} \rightarrow \mathbf{r}} \left\{
   \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})
   \right\}

Parameters
----------
chi_wr : {par_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})` 
                  in imaginary time and momentum space.

Returns
-------
{ret_0}
   Generalized susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{r})` 
            in Matsubara frequency and real space.
)DOC",
                                                  {{}}, {c2py::python_typename<triqs_tprf::chi_wr_t>()});
static const auto _c2py_doc_22 =
   _c2py_fun_22.doc(R"DOC(
Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k})`.

 Computes

 

.. math::

   \chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k}) =
   \chi^{(0)} \left[ 1 - \Gamma^{(PH)} \chi^{(0)} \right]^{-1}

Parameters
----------
chi0_wnk : {par_0}
   Generalized lattice bubble susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.
gamma_ph_wnn : {par_1}
   Local particle-hole vertex function :math:`\Gamma^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`.

Returns
-------
{ret_0}
   Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \nu, \nu', \mathbf{k})`.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_wnk_cvt>()}, {c2py::python_typename<triqs_tprf::chi_wnn_cvt>()}},
                    {c2py::python_typename<triqs_tprf::chi_kwnn_t>()});
static const auto _c2py_doc_23 = _c2py_fun_23.doc(R"DOC()DOC");
static const auto _c2py_doc_24 =
   _c2py_fun_24.doc(R"DOC(
Lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

 Computes

 

.. math::

   \chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k}) =
   \chi^{(0)} \left[ 1 - \Gamma^{(PH)} \chi^{(0)} \right]^{-1}

Parameters
----------
chi0_wnk : {par_0}
   Generalized lattice bubble susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.
gamma_ph_wnn : {par_1}
   Local particle-hole vertex function :math:`\Gamma^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`.

Returns
-------
{ret_0}
   Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_wnk_cvt>()}, {c2py::python_typename<triqs_tprf::chi_wnn_cvt>()}},
                    {c2py::python_typename<triqs_tprf::chi_kw_t>()});
static const auto _c2py_doc_25 = _c2py_fun_25.doc(R"DOC(
Dual lattice Bethe-Salpeter equation solver for the generalized susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.

Parameters
----------
chi0_wnk : {par_0}
   Generalized lattice bubble susceptibility :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.
gamma_ph_wnn : {par_1}
   Local particle-hole vertex function :math:`\Gamma^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu, \nu')`.
L_wn : {par_2}
   Local triangular particle-hole vertex function :math:`L^{(PH)}_{\bar{a}b\bar{c}d}(\omega, \nu)`.

Returns
-------
{ret_0}
   Generalized lattice susceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\omega, \mathbf{k})`.
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::chi_wnk_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_wnn_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_nn_cvt>()}},
                                                  {c2py::python_typename<triqs_tprf::chi_kw_t>()});
static const auto _c2py_doc_26 = _c2py_fun_26.doc(R"DOC()DOC");
static const auto _c2py_doc_27 = _c2py_fun_27.doc(R"DOC()DOC");
static const auto _c2py_doc_28 = _c2py_fun_28.doc(R"DOC()DOC");
static const auto _c2py_doc_29 = _c2py_fun_29.doc(R"DOC()DOC");
static const auto _c2py_doc_30 = _c2py_fun_30.doc(
   R"DOC(
Computes reducible ladder vertex for the approximation of a local and static vertex.

   In this approximation the reducible ladder vertex in density/magnetic channel are given by

   

.. math::

   \Phi^{\text{d/m}}_{a\overline{b}c\overline{d}}(Q)
   &\approx
   \frac{1}{(N_\mathbf{k}\beta)^2}
   \sum_{K'', K'''}
   U^{\text{d/m}}\chi^{\text{d/m}}(Q, K'', K''') U^{\text{d/m}}
   \\
   &\approx
   U^{\mathrm{d/m}}
   \chi^{\text{d/m}}(Q) U^{\mathrm{d/m}}\,,

   where all products are particle-hole products.
   The reducible ladder vertex in then only dependent on one bosonic frequency and momentum.
   It can then be used in :meth:`triqs_tprf.eliashberg.construct_gamma_singlet_rpa`
   or :meth:`triqs_tprf.eliashberg.construct_gamma__rpa` to construct the
   irreducible singlet/triplet vertex.

Parameters
----------
chi : {par_0}
   density/magnetic susceptibility  :math:`\chi^{\mathrm{d/m}}_{\bar{a}b\bar{c}d}(i\omega_n,\mathbf{q})`
U : {par_1}
   density/magnetic local and static vertex  :math:`U^{\mathrm{d/m}}_{a\bar{b}c\bar{d}}`

Returns
-------
{ret_0}
   The reducible ladder vertex in the density/magnetic channel :math:`\Phi^{\mathrm{d/m}}(i\omega_n,\mathbf{q})`
)DOC",
   {{c2py::python_typename<triqs_tprf::chi_wk_vt>()},
    {c2py::python_typename<
       nda::basic_array_view<std::complex<double>, 4, nda::C_layout, 'A', nda::default_accessor, nda::borrowed<nda::mem::AddressSpace::Host>>>()}},
   {c2py::python_typename<triqs_tprf::chi_wk_t>()});
static const auto _c2py_doc_31 = _c2py_fun_31.doc(R"DOC()DOC");
static const auto _c2py_doc_32 = _c2py_fun_32.doc(R"DOC()DOC");
static const auto _c2py_doc_33 = _c2py_fun_33.doc(R"DOC(
Fourier transform Gamma parts to imaginary time and real-space

Parameters
----------
Gamma_pp_dyn_wk : {par_0}
   : The dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`.
Gamma_pp_const_k : {par_1}
   : The part of Gamma that is constant in Matsubara frequency space :math:`\Gamma(\mathbf{k})`.

Returns
-------
[1] : {ret_0}
   Tuple of Gamma_pp_dyn_tr,  the dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`, but now in :math:`\tau`-space, Gamma_pp_const_r, the constant part of Gamma in real-space.

[2] : {ret_1}
   Tuple of Gamma_pp_dyn_tr,  the dynamic part of Gamma, which converges to zero for :math:`\omega_n \rightarrow \infty`, but now in :math:`\tau`-space, Gamma_pp_const_r, the constant part of Gamma in real-space.
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::chi_wk_vt>(), c2py::python_typename<triqs_tprf::chi_Dwk_vt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_k_vt>()}},
                                                  {c2py::python_typename<std::tuple<triqs_tprf::chi_tr_t, triqs_tprf::chi_r_t>>(),
                                                   c2py::python_typename<std::tuple<triqs_tprf::chi_Dtr_t, triqs_tprf::chi_r_t>>()});
static const auto _c2py_doc_34 = _c2py_fun_34.doc(R"DOC(
[1] Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})`.

   The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
   is given by

   

.. math::

   W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
   V_{abcd}(\mathbf{k}) +
   \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
   \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
   W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

------

[3] Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})`.

   The full screened interaction :math:`W(\omega, \mathbf{k})`
   is given by

   

.. math::

   W^{(full)}_{abcd}(\omega, \mathbf{k}) =
   V_{abcd}(\mathbf{k}) +
   \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
   \Pi_{fegh}(\omega, \mathbf{k}) \cdot
   W^{(full)}_{hgcd}(\omega, \mathbf{k})

------

[4] Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(i\omega_n, \mathbf{k})`.

   The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
   is given by

   

.. math::

   W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
   V_{abcd}(i\omega_n, \mathbf{k}) +
   \sum_{efgh} V_{abef}(i\omega_n, \mathbf{k}) \cdot
   \Pi_{fegh}(i\omega_n, \mathbf{k}) \cdot
   W^{(full)}_{hgcd}(i\omega_n, \mathbf{k})

------

[6] Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(\omega, \mathbf{k})`.

   The full screened interaction :math:`W(\omega, \mathbf{k})`
   is given by

   

.. math::

   W^{(full)}_{abcd}(\omega, \mathbf{k}) =
   V_{abcd}(\omega, \mathbf{k}) +
   \sum_{efgh} V_{abef}(\omega, \mathbf{k}) \cdot
   \Pi_{fegh}(\omega, \mathbf{k}) \cdot
   W^{(full)}_{hgcd}(\omega, \mathbf{k})

------

Parameters
----------
PI_wk : {par_0}
   polarization bubble :math:`\Pi_{abcd}(i\omega_n, \mathbf{k})`
V_k : {par_1}
   static bare interaction :math:`V_{abcd}(\mathbf{k})`
PI_fk : {par_2}
   polarization bubble :math:`\Pi_{abcd}(\omega, \mathbf{k})`
V_wk : {par_3}
   dynamic bare interaction :math:`V_{abcd}(i\omega_n, \mathbf{k})`
V_fk : {par_4}
   dynamic bare interaction :math:`V_{abcd}(\omega, \mathbf{k})`

Returns
-------
[1, 4] : {ret_0}
   dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`

[3, 6] : {ret_1}
   dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::chi_wk_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_k_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_fk_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_wk_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_fk_cvt>()}},
                                                  {c2py::python_typename<triqs_tprf::chi_wk_t>(), c2py::python_typename<triqs_tprf::chi_fk_t>()});
static const auto _c2py_doc_35 = _c2py_fun_35.doc(R"DOC(
[1] Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})` and known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

   The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
   is given by

   

.. math::

   W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
   V_{abcd}(\mathbf{k}) +
   \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
   \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
   V_{hgcd}(\mathbf{k})

------

[3] Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for static momentum-dependent bare interactions :math:`V(\mathbf{k})` and known generalized susceptibility :math:`\chi(\omega, \mathbf{k})`

   The full screened interaction :math:`W(\omega, \mathbf{k})`
   is given by

   

.. math::

   W^{(full)}_{abcd}(\omega, \mathbf{k}) =
   V_{abcd}(\mathbf{k}) +
   \sum_{efgh} V_{abef}(\mathbf{k}) \cdot
   \chi_{fegh}(\omega, \mathbf{k}) \cdot
   V_{hgcd}(\mathbf{k})

------

[4] Dynamical screened interaction :math:`W(i\omega_n, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(i\omega_n, \mathbf{k})` and known generalized susceptibility :math:`\chi(i\omega_n, \mathbf{k})`

   The full screened interaction :math:`W(i\omega_n, \mathbf{k})`
   is given by

   

.. math::

   W^{(full)}_{abcd}(i\omega_n, \mathbf{k}) =
   V_{abcd}(i\omega_n, \mathbf{k}) +
   \sum_{efgh} V_{abef}(i\omega_n, \mathbf{k}) \cdot
   \chi_{fegh}(i\omega_n, \mathbf{k}) \cdot
   V_{hgcd}(i\omega_n, \mathbf{k})

------

[6] Dynamical screened interaction :math:`W(\omega, \mathbf{k})` calculator for dynamic momentum-dependent bare interactions :math:`V(\omega, \mathbf{k})` and known generalized susceptibility :math:`\chi(\omega, \mathbf{k})`

   The full screened interaction :math:`W(\omega, \mathbf{k})`
   is given by

   

.. math::

   W^{(full)}_{abcd}(\omega, \mathbf{k}) =
   V_{abcd}(\omega, \mathbf{k}) +
   \sum_{efgh} V_{abef}(\omega, \mathbf{k}) \cdot
   \chi_{fegh}(\omega, \mathbf{k}) \cdot
   V_{hgcd}(\omega, \mathbf{k})

------

Parameters
----------
chi_wk : {par_0}
   generalized susceptibility :math:`\chi_{abcd}(i\omega_n, \mathbf{k})`
V_k : {par_1}
   static bare interaction :math:`V_{abcd}(\mathbf{k})`
chi_fk : {par_2}
   generalized susceptibility :math:`\chi_{abcd}(\omega, \mathbf{k})`
V_wk : {par_3}
   dynamic bare interaction :math:`V_{abcd}(i\omega_n, \mathbf{k})`
V_fk : {par_4}
   dynamic bare interaction :math:`V_{abcd}(\omega, \mathbf{k})`

Returns
-------
[1, 4] : {ret_0}
   dynamical screened interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`

[3, 6] : {ret_1}
   dynamical screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::chi_wk_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_k_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_fk_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_wk_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_fk_cvt>()}},
                                                  {c2py::python_typename<triqs_tprf::chi_wk_t>(), c2py::python_typename<triqs_tprf::chi_fk_t>()});
static const auto _c2py_doc_36 = _c2py_fun_36.doc(R"DOC()DOC");
static const auto _c2py_doc_37 = _c2py_fun_37.doc(R"DOC()DOC");
static const auto _c2py_doc_38 = _c2py_fun_38.doc(
   R"DOC(
Linearized Eliashberg product via summation

    Computes the linearized Eliashberg product in the singlet/triplet channel given by

    

.. math::

   \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu,\mathbf{k})
   =
   -\frac{1}{2N_\mathbf{k} \beta}\sum_{i\nu'}\sum_{\mathbf{k}'}
   \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(i\nu - i\nu',\mathbf{k}-\mathbf{k}')
   \\
   \times
   G_{c\bar{e}}(i\nu',\mathbf{k}')
   G_{d\bar{f}}(-i\nu',-\mathbf{k}')
   \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{e}\bar{f}}(i\nu',\mathbf{k}')\,,

    by summation.

Parameters
----------
Gamma_pp : {par_0}
   particle-particle vertex :math:`\Gamma^{\mathrm{s/t}}_{a\bar{b}c\bar{d}}(i\nu_n,\mathbf{k})`
g_wk : {par_1}
   single particle Green's function :math:`G_{a\bar{b}}(i\nu_n,\mathbf{k})`
delta_wk : {par_2}
   superconducting gap :math:`\Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})`

Returns
-------
{ret_0}
   Gives the result of the product :math:`\Delta^{\mathrm{s/t}, \mathrm{out}}`
)DOC",
   {{c2py::python_typename<triqs_tprf::chi_wk_vt>()}, {c2py::python_typename<triqs_tprf::g_wk_vt>()}, {c2py::python_typename<triqs_tprf::g_wk_vt>()}},
   {c2py::python_typename<triqs_tprf::g_wk_t>()});
static const auto _c2py_doc_39 = _c2py_fun_39.doc(R"DOC(
Linearized Eliashberg product via FFT

    Computes the linearized Eliashberg product in the singlet/triplet channel given by

    

.. math::

   \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu,\mathbf{k})
   =
   -\frac{1}{2N_\mathbf{k} \beta}\sum_{i\nu'}\sum_{\mathbf{k}'}
   \Gamma^{\mathrm{s/t}}_{c\bar{a}d\bar{b}}(i\nu - i\nu',\mathbf{k}-\mathbf{k}')
   \\
   \times
   G_{c\bar{e}}(i\nu',\mathbf{k}')
   G_{d\bar{f}}(-i\nu',-\mathbf{k}')
   \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{e}\bar{f}}(i\nu',\mathbf{k}')\,,

    by taking advantage of the convolution theorem.

    We therefore first calculate

    

.. math::

   F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
   =
   G_{a\bar{c}}(i\nu,\mathbf{k})
   G_{b\bar{d}}(-i\nu,-\mathbf{k})
   \Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{c}\bar{d}}(i\nu,\mathbf{k})\,,

    which we then Fourier transform to imaginary time and real-space

    

.. math::

   F^{\mathrm{s/t}}_{ab}(\tau,\mathbf{r})
   =
   \mathcal{F}^2
   \big(
   F^{\mathrm{s/t}}_{ab}(i\nu,\mathbf{k})
   \big)\,.

    We then calculate first the dynamic gap

    

.. math::

   \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(\tau,\mathbf{r})
   =
   -\frac{1}{2}
   \Gamma^{\mathrm{s/t}, \mathrm{dynamic}}_{c\bar{a}d\bar{b}}(\tau, \mathbf{r})
   F^{\mathrm{s/t}}_{cd}(\tau, \mathbf{r})\,,

    and then the static gap

    

.. math::

   \Delta^{\mathrm{s/t}, \mathrm{static}}_{\bar{a}\bar{b}}(\mathbf{r})
   =
   -\frac{1}{2}
   \Gamma^{\mathrm{s/t}, \mathrm{static}}_{c\bar{a}d\bar{b}}(\mathbf{r})
   F^{\mathrm{s/t}}_{cd}(\tau=0, \mathbf{r})\,.

    We then Fourier transform the dynamic gap to imaginary frequencies

    

.. math::

   \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
   =
   \mathcal{F}
   \big(
   \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(\tau,\mathbf{r})
   \big)\,,

    and then add both component together

    

.. math::

   \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
   =
   \Delta^{\mathrm{s/t}, \mathrm{dynamic}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
   +
   \Delta^{\mathrm{s/t}, \mathrm{static}}_{\bar{a}\bar{b}}(\mathbf{r})\,,

   and then finally Fourier transform to :math:`\mathbf{k}`-space

   

.. math::

   \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})
   =
   \mathcal{F}
   \big(
   \Delta^{\mathrm{s/t}, \mathrm{out}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{r})
   \big)\,.

Parameters
----------
Gamma_pp_dyn_tr : {par_0}
   dynamic part of the particle-particle vertex :math:`\Gamma^{\mathrm{s/t}, \mathrm{dynamic}}_{c\bar{a}d\bar{b}}(\tau, \mathbf{r})`
Gamma_pp_const_r : {par_1}
   static part of the particle-particle vertex :math:`\Gamma^{\mathrm{s/t}, \mathrm{static}}_{c\bar{a}d\bar{b}}(\mathbf{r})`
g_wk : {par_2}
   one-particle Green's function :math:`G_{a\bar{b}}(i\nu_n,\mathbf{k})`
delta_wk : {par_3}
   superconducting gap :math:`\Delta^{\mathrm{s/t}, \mathrm{in}}_{\bar{a}\bar{b}}(i\nu_n,\mathbf{k})`

Returns
-------
{ret_0}
   Gives the result of the product :math:`\Delta^{\mathrm{s/t}, \mathrm{out}}`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::chi_tr_vt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_r_vt>()},
                                                   {c2py::python_typename<triqs_tprf::g_wk_vt>()},
                                                   {c2py::python_typename<triqs_tprf::g_wk_vt>()}},
                                                  {c2py::python_typename<triqs_tprf::g_wk_t>()});
static const auto _c2py_doc_40 = _c2py_fun_40.doc(R"DOC()DOC");
static const auto _c2py_doc_41 = _c2py_fun_41.doc(R"DOC(
Helper function to evaluate the Fermi-Dirac distribution function

 .. math ::
     f() = {1}{() + 1}

Parameters
----------
e : {par_0}
   : point at which to evaluate :math:`f(\epsilon)`.

Returns
-------
{ret_0}
   The value of :math:`f(\epsilon)`.
)DOC",
                                                  {{c2py::python_typename<double>()}}, {c2py::python_typename<double>()});
static const auto _c2py_doc_42 =
   _c2py_fun_42.doc(R"DOC(
Fock self energy :math:`\Sigma_{ab}(\mathbf{k})` calculator

   Computes the Fock self-energy of a static interaction as the sum

   

.. math::

   \Sigma_{ab}(\mathbf{k}) = -\frac{1}{N_k}
   \sum_{\mathbf{q},cd} V_{acdb}(\mathbf{q}) \rho_{dc}(\mathbf{k} + \mathbf{q})

   where :math:`\rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k})` is the density matrix of the
   single particle Green's function.

Parameters
----------
V_k : {par_0}
   static interaction :math:`V_{abcd}(\mathbf{q})`
g_wk : {par_1}
   single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
{ret_0}
   Fock self-energy :math:`\Sigma_{ab}(\mathbf{k})`
)DOC",
                    {{}, {c2py::python_typename<triqs_tprf::g_wk_cvt>()}}, {c2py::python_typename<triqs_tprf::e_k_t>()});
static const auto _c2py_doc_43 = _c2py_fun_43.doc(R"DOC(
Inverse fast fourier transform of real time Green's function from k-space to real space

   Computes: :math:`G_{a\bar{b}}(t, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(t, \mathbf{k})\right\}`

Parameters
----------
g_Tk : {par_0}
   k-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{k})`

Returns
-------
{ret_0}
   real-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{r})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_Tk_cvt>()}}, {c2py::python_typename<triqs_tprf::g_Tr_t>()});
static const auto _c2py_doc_44 = _c2py_fun_44.doc(R"DOC(
[1] Fast fourier transform of real time Green's function from real-space to k-space

   Computes: :math:`G_{a\bar{b}}(t, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(t, \mathbf{r}) \right\}`

------

[2] Fast fourier transform of real time Green's function from real-space to k-space

   Computes: :math:`G_{a\bar{b}c\bar{d}}(t, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}c\bar{d}}(t, \mathbf{r}) \right\}`

------

Parameters
----------
g_Tr : {par_0}
   real-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{r})`

Returns
-------
[1] : {ret_0}
   k-space real time Green's function :math:`G_{a\bar{b}}(t, \mathbf{k})`

[2] : {ret_1}
   k-space real time Green's function :math:`G_{a\bar{b}c\bar{d}}(t, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_Tr_cvt>()}},
                                                  {c2py::python_typename<triqs_tprf::g_Tk_t>(), c2py::python_typename<triqs_tprf::chi_Tk_t>()});
static const auto _c2py_doc_45 = _c2py_fun_45.doc(R"DOC(
Inverse fast fourier transform of real frequency Green's function from k-space to real space

   Computes: :math:`G_{a\bar{b}}(\omega, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(\omega, \mathbf{k})\right\}`

Parameters
----------
g_fk : {par_0}
   k-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`

Returns
-------
{ret_0}
   real-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{r})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_fk_cvt>()}}, {c2py::python_typename<triqs_tprf::g_fr_t>()});
static const auto _c2py_doc_46 = _c2py_fun_46.doc(R"DOC(
Fast fourier transform of real frequency Green's function from real-space to k-space

   Computes: :math:`G_{a\bar{b}}(\omega, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(\omega, \mathbf{r}) \right\}`

Parameters
----------
g_fr : {par_0}
   real-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{r})`

Returns
-------
{ret_0}
   k-space real frequency Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_fr_cvt>()}}, {c2py::python_typename<triqs_tprf::g_fk_t>()});
static const auto _c2py_doc_47 = _c2py_fun_47.doc(R"DOC(
Fast fourier transform of real-space Green's function from imaginary time to Matsubara frequency

   Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(\tau, \mathbf{r}) \right\}`

Parameters
----------
g_tr : {par_0}
   real-space imaginary time Green's function :math:`G_{a\bar{b}}(\tau, \mathbf{r})`

Returns
-------
{ret_0}
   real-space Matsubara frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_tr_cvt>()}}, {c2py::python_typename<triqs_tprf::g_wr_t>()});
static const auto _c2py_doc_48 = _c2py_fun_48.doc(R"DOC(
Inverse fast fourier transform of imaginary frequency Green's function from k-space to real space

   Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r}) = \mathcal{F}^{-1} \left\{G_{a\bar{b}}(i\omega_n, \mathbf{k})\right\}`

Parameters
----------
g_wk : {par_0}
   k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`

Returns
-------
{ret_0}
   real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_wk_cvt>()}}, {c2py::python_typename<triqs_tprf::g_wr_t>()});
static const auto _c2py_doc_49 = _c2py_fun_49.doc(R"DOC(
Fast fourier transform of real-space Green's function from Matsubara frequency to imaginary time

   Computes: :math:`G_{a\bar{b}}(\tau, \mathbf{r}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`

Parameters
----------
g_wr : {par_0}
   real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`

Returns
-------
{ret_0}
   real-space imaginary time Green's function :math:`G_{a\bar{b}}(\tau, \mathbf{r})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_wr_cvt>()}}, {c2py::python_typename<triqs_tprf::g_tr_t>()});
static const auto _c2py_doc_50 = _c2py_fun_50.doc(R"DOC(
Fast fourier transform of imaginary frequency Green's function from real-space to k-space

   Computes: :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \mathcal{F} \left\{ G_{a\bar{b}}(i\omega_n, \mathbf{r}) \right\}`

Parameters
----------
g_wr : {par_0}
   real-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{r})`

Returns
-------
{ret_0}
   k-space imaginary frequency Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_wr_cvt>()}}, {c2py::python_typename<triqs_tprf::g_wk_t>()});
static const auto _c2py_doc_51 = _c2py_fun_51.doc(R"DOC()DOC");
static const auto _c2py_doc_52 = _c2py_fun_52.doc(R"DOC(
[1, 2] some documentation

------

[3] Real frequency GW self energy :math:`\Sigma(\omega, \mathbf{k})` calculator via the spectral representation

   Computes the spectral function of the dynamic part of the screened interaction

   

.. math::

   W^{(spec)}_{ab}(\omega, \mathbf{k}) = \frac{-1}{\pi} \text{Im}
   \left( W_{aabb}(\omega, \mathbf{k}) - V_{aabb}(\mathbf{k}) \right)

         
   and constructs the dynamic part of the GW self energy via the spectral representation

   

.. math::

   \Sigma_{ab}(\omega, \mathbf{k}) = \frac{\delta_{\omega}}{N_k} \sum_{\mathbf{q}} \sum_{\omega'}
   U_{al}(\mathbf{k}+\mathbf{q}) U^{\dagger}_{lb}(\mathbf{k}+\mathbf{q})
   W^{(spec)}_{ab}(\omega', \mathbf{q})
   \frac{n_B(\omega') + f(\epsilon_{\mathbf{k}+\mathbf{q}, l})}{\omega + i\delta + \omega' - \epsilon_{\mathbf{k}+\mathbf{q}, l} + \mu}

         
   where :math:`\delta_{\omega}` is the real-frequency mesh spacing and the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued 
   dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

   

.. math::

   \sum_{\bar{a}b} U^\dagger_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U_{bj} (\mathbf{k})
   = \delta_{ij} \epsilon_{\mathbf{k}, i}

------

Parameters
----------
mu : {par_0}
   chemical potential :math:`\mu`
beta : {par_1}
   inverse temperature
e_k : {par_2}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
W_fk : {par_3}
   fully screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
V_k : {par_4}
   bare interaction :math:`V_{abcd}(\mathbf{k})`
delta : {par_5}
   broadening :math:`\delta`

Returns
-------
{ret_0}
   real frequency GW self-energy :math:`\Sigma_{ab}(\omega, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<double>()},
                                                   {c2py::python_typename<double>()},
                                                   {c2py::python_typename<triqs_tprf::e_k_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::chi_fk_cvt>()},
                                                   {},
                                                   {c2py::python_typename<double>()}},
                                                  {c2py::python_typename<triqs_tprf::g_fk_t>()});
static const auto _c2py_doc_53 = _c2py_fun_53.doc(R"DOC(
[1, 2] Some documentation

------

[3] GW self energy :math:`\Sigma(\mathbf{k})` calculator for static interactions

   Computes the GW self-energy of a static interaction as the product

   

.. math::

   \Sigma_{ab}(\mathbf{k}) = \frac{-1}{N_k} \sum_{\mathbf{q}} \sum_{l}
   U_{al}(\mathbf{k}+\mathbf{q}) U^\dagger_{lb}(\mathbf{k}+\mathbf{q})
   V_{aabb}(\mathbf{q}) f(\epsilon_{\mathbf{k}+\mathbf{q}, l})

   where the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued 
   dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

   

.. math::

   \sum_{\bar{a}b} U^\dagger_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U_{bj} (\mathbf{k})
   = \delta_{ij} \epsilon_{\mathbf{k}, i}

------

[4, 5] some documentation

------

[6] Real frequency GW self energy :math:`\Sigma(\omega, \mathbf{k})` calculator via the spectral representation

   Computes the spectral function of the dynamic part of the screened interaction

   

.. math::

   W^{(spec)}_{ab}(\omega, \mathbf{k}) = \frac{-1}{\pi} \text{Im}
   \left( W_{aabb}(\omega, \mathbf{k}) - V_{aabb}(\mathbf{k}) \right)

         
   and constructs the GW self energy via the spectral representation

   

.. math::

   \Sigma_{ab}(\omega, \mathbf{k}) = \frac{-1}{N_k} \sum_{\mathbf{q}} \sum_{l}
   U_{al}(\mathbf{k}+\mathbf{q}) U^{\dagger}_{lb}(\mathbf{k}+\mathbf{q})
   V_{aabb}(\mathbf{q}) f(\epsilon_{\mathbf{k}+\mathbf{q}, l}) \\
   + \frac{\delta_{\omega}}{N_k} \sum_{\mathbf{q}} \sum_{\omega'}
   U_{al}(\mathbf{k}+\mathbf{q}) U^{\dagger}_{lb}(\mathbf{k}+\mathbf{q})
   W^{(spec)}_{ab}(\omega', \mathbf{q})
   \frac{n_B(\omega') + f(\epsilon_{\mathbf{k}+\mathbf{q}, l})}{\omega + i\delta + \omega' - \epsilon_{\mathbf{k}+\mathbf{q}, l} + \mu}

         
   where :math:`\delta_{\omega}` is the real-frequency mesh spacing and the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued 
   dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

   

.. math::

   \sum_{\bar{a}b} U^\dagger_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U_{bj} (\mathbf{k})
   = \delta_{ij} \epsilon_{\mathbf{k}, i}

------

Parameters
----------
mu : {par_0}
   chemical potential :math:`\mu`
beta : {par_1}
   inverse temperature
e_k : {par_2}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
V_k : {par_3}
   bare interaction :math:`V_{abcd}(\mathbf{k})`
W_fk : {par_4}
   fully screened interaction :math:`W_{abcd}(\omega, \mathbf{k})`
delta : {par_5}
   broadening :math:`\delta`

Returns
-------
[3] : {ret_0}
   static GW self-energy :math:`\Sigma_{ab}(\mathbf{k})`

[6] : {ret_1}
   real frequency GW self-energy :math:`\Sigma_{ab}(\omega, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<double>()},
                                                   {c2py::python_typename<double>()},
                                                   {c2py::python_typename<triqs_tprf::e_k_cvt>()},
                                                   {},
                                                   {c2py::python_typename<triqs_tprf::chi_fk_cvt>()},
                                                   {c2py::python_typename<double>()}},
                                                  {c2py::python_typename<triqs_tprf::e_k_t>(), c2py::python_typename<triqs_tprf::g_fk_t>()});
static const auto _c2py_doc_54 =
   _c2py_fun_54.doc(R"DOC(
Dynamic GW self energy :math:`\Sigma(\tau, \mathbf{r})` calculator 

   Computes the GW self-energy as the product

   

.. math::

   \Sigma_{ab}(\tau, \mathbf{r}) =
   - \sum_{cd} W_{abcd}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

Parameters
----------
W_tr : {par_0}
   interaction :math:`W_{abcd}(\tau, \mathbf{r})`
g_tr : {par_1}
   single particle Green's function :math:`G_{ab}(\tau, \mathbf{r})`

Returns
-------
{ret_0}
   Dynamic GW self-energy :math:`\Sigma_{ab}(\tau, \mathbf{r})`
)DOC",
                    {{c2py::python_typename<triqs_tprf::chi_tr_cvt>()}, {c2py::python_typename<triqs_tprf::g_tr_cvt>()}},
                    {c2py::python_typename<triqs_tprf::g_tr_t>()});
static const auto _c2py_doc_55 = _c2py_fun_55.doc(
   R"DOC(
[1] GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator for dynamic interactions

   Splits the interaction into a dynamic and a static part

   .. math ::
       W_{abcd}(i, {k}) = 
           W^{(dyn)}_{abcd}(i, {k})
           + V_{abcd}({k})

   by fitting the high-frequency tail.

   Fourier transforms the dynamic part of the interaction and the 
   single-particle Green's function to imaginary time and real space.

   

.. math::

   G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
   \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

   

.. math::

   W^{(dyn)}_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
   \left\{ W^{(dyn)}_{abcd}(i\omega_n, \mathbf{k}) \right\}

   computes the GW self-energy as the product

   

.. math::

   \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) =
   - \sum_{cd} W^{(dyn)}_{acdb}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

   and transforms back to frequency and momentum

   

.. math::

   \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k}) =
   \mathcal{F} \left\{ \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) \right\}

   The self-energy of the static part of the interaction is calculated
   as the sum

   

.. math::

   \Sigma^{(stat)}_{ab}(\mathbf{k}) = -\frac{1}{N_k}
   \sum_{\mathbf{q},cd} V_{acdb}(\mathbf{k}) \rho_{dc}(\mathbf{k} + \mathbf{q})

   where :math:`\rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k})` is the density matrix of the
   single particle Green's function.

   The total GW self-energy is given by

   

.. math::

   \Sigma_{ab}(i\omega_n, \mathbf{k}) =
   \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k})
   + \Sigma^{(stat)}_{ab}(\mathbf{k})

------

[2] GW self energy :math:`\Sigma(i\omega_n, \mathbf{k})` calculator for dynamic interactions

   Fourier transforms the dynamic part of the interaction and the 
   single-particle Green's function to imaginary time and real space.

   

.. math::

   G_{ab}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
   \left\{ G_{ab}(i\omega_n, \mathbf{k}) \right\}

   

.. math::

   W^{(dyn)}_{abcd}(\tau, \mathbf{r}) = \mathcal{F}^{-1}
   \left\{ W^{(dyn)}_{abcd}(i\omega_n, \mathbf{k}) \right\}

   computes the GW self-energy as the product

   

.. math::

   \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) =
   - \sum_{cd} W^{(dyn)}_{acdb}(\tau, \mathbf{r}) G_{cd}(\tau, \mathbf{r})

   and transforms back to frequency and momentum

   

.. math::

   \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k}) =
   \mathcal{F} \left\{ \Sigma^{(dyn)}_{ab}(\tau, \mathbf{r}) \right\}

   The self-energy of the static part of the interaction is calculated
   as the sum

   

.. math::

   \Sigma^{(stat)}_{ab}(\mathbf{k}) = -\frac{1}{N_k}
   \sum_{\mathbf{q},cd} V_{acdb}(\mathbf{k}) \rho_{dc}(\mathbf{k} + \mathbf{q})

   where :math:`\rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k})` is the density matrix of the
   single particle Green's function.

   The total GW self-energy is given by

   

.. math::

   \Sigma_{ab}(i\omega_n, \mathbf{k}) =
   \Sigma^{(dyn)}_{ab}(i\omega_n, \mathbf{k})
   + \Sigma^{(stat)}_{ab}(\mathbf{k})

------

[3] Static GW self energy :math:`\Sigma_{ab}(\mathbf{k})` calculator

   Computes the static GW self-energy (equivalent to the Fock self-energy)

------

Parameters
----------
W_wk : {par_0}
   interaction :math:`W_{abcd}(i\omega_n, \mathbf{k})`
g_wk : {par_1}
   single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`
V_k : {par_2}
   static interaction :math:`V_{abcd}(\mathbf{q})`

Returns
-------
[1] : {ret_0}
   GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`

[2] : {ret_1}
   GW self-energy :math:`\Sigma_{ab}(i\omega_n, \mathbf{k})`

[3] : {ret_2}
   Static GW self-energy (Fock) :math:`\Sigma_{ab}(\mathbf{k})`
)DOC",
   {{c2py::python_typename<triqs_tprf::chi_wk_cvt>(), c2py::python_typename<triqs_tprf::chi_Dwk_cvt>()},
    {c2py::python_typename<triqs_tprf::g_wk_cvt>(), c2py::python_typename<triqs_tprf::g_Dwk_cvt>()},
    {}},
   {c2py::python_typename<triqs_tprf::g_wk_t>(), c2py::python_typename<triqs_tprf::g_Dwk_t>(), c2py::python_typename<triqs_tprf::e_k_t>()});
static const auto _c2py_doc_56 =
   _c2py_fun_56.doc(R"DOC(
Hartree self energy :math:`\Sigma_{ab}(\mathbf{k})` calculator

   Computes the Hartree self-energy of a static interaction as the sum

   

.. math::

   \Sigma_{ab}(\mathbf{k}) = \frac{1}{N_k}
   \sum_{\mathbf{q},cd} V_{abcd}(\mathbf{q}) \rho_{cd}(\mathbf{k} + \mathbf{q})

   where :math:`\rho_{ab}(\mathbf{k}) = -G_{ba}(\beta, \mathbf{k})` is the density matrix of the
   single particle Green's function.

Parameters
----------
V_k : {par_0}
   static interaction :math:`V_{abcd}(\mathbf{q})`
g_wk : {par_1}
   single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
{ret_0}
   Hartree self-energy :math:`\Sigma_{ab}(\mathbf{k})`
)DOC",
                    {{}, {c2py::python_typename<triqs_tprf::g_wk_cvt>()}}, {c2py::python_typename<triqs_tprf::e_k_t>()});
static const auto _c2py_doc_57 = _c2py_fun_57.doc(R"DOC(
Construct a non-interacting real frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(\omega, \mathbf{k})`

 Computes

 

.. math::

   G^{(0)}_{a\bar{b}}(\omega, \mathbf{k}) = \left[
   (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
   \right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, chemical potential :math:`\mu`,
 broadening :math:`\delta`, and a real frequency Green's function mesh.

Parameters
----------
mu : {par_0}
   chemical potential :math:`\mu`
e_k : {par_1}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
mesh : {par_2}
   real frequency mesh
delta : {par_3}
   broadening :math:`\delta`

Returns
-------
{ret_0}
   Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(\omega, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<double>()},
                                                   {c2py::python_typename<triqs_tprf::e_k_cvt>()},
                                                   {c2py::python_typename<triqs::mesh::refreq>()},
                                                   {c2py::python_typename<double>()}},
                                                  {c2py::python_typename<triqs_tprf::g_fk_t>()});
static const auto _c2py_doc_58 = _c2py_fun_58.doc(R"DOC(
Construct a non-interacting Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})`

 Computes

 

.. math::

   G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
   (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k})
   \right]^{-1}_{a\bar{b}},

 using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, chemical potential :math:`\mu`,
 and a Matsubara frequency Green's function mesh.

Parameters
----------
mu : {par_0}
   chemical potential :math:`\mu`
e_k : {par_1}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
mesh : {par_2}
   imaginary frequency mesh

Returns
-------
[1] : {ret_0}
   Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})`

[2] : {ret_1}
   Matsubara frequency lattice Green's function :math:`G^{(0)}_{a\bar{b}}(i\omega_n, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<double>()},
                                                   {c2py::python_typename<triqs_tprf::e_k_cvt>()},
                                                   {c2py::python_typename<triqs::mesh::imfreq>(), c2py::python_typename<triqs::mesh::dlr_imfreq>()}},
                                                  {c2py::python_typename<triqs_tprf::g_wk_t>(), c2py::python_typename<triqs_tprf::g_Dwk_t>()});
static const auto _c2py_doc_59 = _c2py_fun_59.doc(R"DOC(
Construct an interacting real frequency local (:math:`\mathbf{r}=\mathbf{0}`) lattice Green's function :math:`G_{a\bar{b}}(\omega)`
  
Computes

.. math::

   G_{a\bar{b}}(\omega) = \frac{1}{N_k} \sum_\mathbf{k} \left[
   (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega)
   \right]^{-1}_{a\bar{b}},

using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, 
chemical potential :math:`\mu`, and a momentum independent real frequency 
self energy :math:`\Sigma_{\bar{a}b}(\omega)`.

Parameters
----------
mu : {par_0}
   chemical potential :math:`\mu`
e_k : {par_1}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
sigma_f : {par_2}
   real frequency self-energy :math:`\Sigma_{\bar{a}b}(\omega)`
delta : {par_3}
   broadening :math:`\delta`

Returns
-------
{ret_0}
   Real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<double>()},
                                                   {c2py::python_typename<triqs_tprf::e_k_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::g_f_cvt>()},
                                                   {c2py::python_typename<double>()}},
                                                  {c2py::python_typename<triqs_tprf::g_f_t>()});
static const auto _c2py_doc_60 = _c2py_fun_60.doc(R"DOC(
[1] Construct an interacting real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
  
Computes

.. math::

   G_{a\bar{b}}(\omega, \mathbf{k}) = \left[
   (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega, \mathbf{k})
   \right]^{-1}_{a\bar{b}},

using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, 
chemical potential :math:`\mu`, broadening :math:`\delta`, and a real frequency 
self energy :math:`\Sigma_{\bar{a}b}(\omega, \mathbf{k})`.

------

[2] Construct an interacting real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
  
Computes

.. math::

   G_{a\bar{b}}(\omega, \mathbf{k}) = \left[
   (\omega + i\delta + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(\omega)
   \right]^{-1}_{a\bar{b}},

using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, 
chemical potential :math:`\mu`, broadening :math:`\delta`, and a real frequency 
self energy :math:`\Sigma_{\bar{a}b}(\omega)`.

------

Parameters
----------
mu : {par_0}
   chemical potential :math:`\mu`
e_k : {par_1}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
sigma_fk : {par_2}
   real frequency self-energy :math:`\Sigma_{\bar{a}b}(\omega, \mathbf{k})`
delta : {par_3}
   broadening :math:`\delta`

Returns
-------
{ret_0}
   real frequency lattice Green's function :math:`G_{a\bar{b}}(\omega, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<double>()},
                                                   {c2py::python_typename<triqs_tprf::e_k_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::g_fk_cvt>()},
                                                   {c2py::python_typename<double>()}},
                                                  {c2py::python_typename<triqs_tprf::g_fk_t>()});
static const auto _c2py_doc_61 = _c2py_fun_61.doc(
   R"DOC(
Construct an interacting Matsubara frequency local (:math:`\mathbf{r}=\mathbf{0}`) lattice Green's function :math:`G_{a\bar{b}}(i\omega_n)`
  
Computes

.. math::

   G_{a\bar{b}}(i\omega_n) = \frac{1}{N_k} \sum_\mathbf{k} \left[
   (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
   \right]^{-1}_{a\bar{b}},

using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, 
chemical potential :math:`\mu`, and a momentum independent Matsubara frequency 
self energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`.

Parameters
----------
mu : {par_0}
   chemical potential :math:`\mu`
e_k : {par_1}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
sigma_w : {par_2}
   imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`

Returns
-------
{ret_0}
   Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
)DOC",
   {{c2py::python_typename<double>()}, {c2py::python_typename<triqs_tprf::e_k_cvt>()}, {c2py::python_typename<triqs_tprf::g_w_cvt>()}},
   {c2py::python_typename<triqs_tprf::g_w_t>()});
static const auto _c2py_doc_62 = _c2py_fun_62.doc(R"DOC(
[1] Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
  
Computes

.. math::

   G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
   (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
   \right]^{-1}_{a\bar{b}},

using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, 
chemical potential :math:`\mu`, and a momentum independent Matsubara frequency 
self energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`.

------

[2, 3] Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
  
Computes

.. math::

   G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
   (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n, \mathbf{k})
   \right]^{-1}_{a\bar{b}},

using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, 
chemical potential :math:`\mu`, and a Matsubara frequency 
self energy :math:`\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})`.

------

[4] Construct an interacting Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
  
Computes

.. math::

   G_{a\bar{b}}(i\omega_n, \mathbf{k}) = \left[
   (i\omega_n + \mu ) \cdot \mathbf{1}  - \epsilon(\mathbf{k}) - \Sigma(i\omega_n)
   \right]^{-1}_{a\bar{b}},

using a discretized dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, 
chemical potential :math:`\mu`, and a Matsubara frequency 
self energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`.

------

Parameters
----------
mu : {par_0}
   chemical potential :math:`\mu`
e_k : {par_1}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
sigma_w : {par_2}
   imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n)`
sigma_wk : {par_3}
   imaginary frequency self-energy :math:`\Sigma_{\bar{a}b}(i\omega_n, \mathbf{k})`

Returns
-------
[1, 2] : {ret_0}
   Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`

[3, 4] : {ret_1}
   Matsubara frequency lattice Green's function :math:`G_{a\bar{b}}(i\omega_n, \mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<double>()},
                                                   {c2py::python_typename<triqs_tprf::e_k_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::g_w_cvt>(), c2py::python_typename<triqs_tprf::g_Dw_cvt>()},
                                                   {c2py::python_typename<triqs_tprf::g_wk_cvt>(), c2py::python_typename<triqs_tprf::g_Dwk_cvt>()}},
                                                  {c2py::python_typename<triqs_tprf::g_wk_t>(), c2py::python_typename<triqs_tprf::g_Dwk_t>()});
static const auto _c2py_doc_63 = _c2py_fun_63.doc(R"DOC(
[1] Generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q})`.

   Analytic calculation of the generalized (non-interacting) Lindhard susceptibility 
   in the particle-hole channel. The analytic expression is obtained using residue calculus 
   to explicitly evaluate the matsubara sum of the fourier transformed imaginary time
   bubble product of two non-interacting single-particle Green's functions.

   

.. math::

   G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n) =
   \left[ i\omega_n \cdot \mathbf{1} - \epsilon(\mathbf{k}) \right]^{-1} .

   The analytic evaluation of the bubble diagram gives

   

.. math::

   \chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q}) \equiv
   \mathcal{F} \left\{
   - G^{(0)}_{d\bar{a}}(\tau, \mathbf{r}) G^{(0)}_{b\bar{c}}(-\tau, -\mathbf{r})
   \right\}
   =
   - \frac{1}{N_k} \sum_{\nu} \sum_{\mathbf{k}}
   G^{(0)}_{d\bar{a}}(\nu, \mathbf{k})
   G^{(0)}_{b\bar{c}}(\nu + \omega, \mathbf{k} + \mathbf{q})
   \\ =
   - \frac{1}{N_k} \sum_{\nu} \sum_{\mathbf{k}}
   \left( \sum_{i}
   U^\dagger_{di}(\mathbf{k}) \frac{1}{i\nu - \epsilon_{\mathbf{k}, i}} U_{i\bar{a}}(\mathbf{k})
   \right)
   \left( \sum_j
   U^\dagger_{bj}(\mathbf{k} + \mathbf{q})
   \frac{1}{i\nu + i\omega - \epsilon_{\mathbf{k} + \mathbf{q}, j}}
   U_{j\bar{c}}(\mathbf{k} + \mathbf{q})
   \right)
   \\ =
   \frac{1}{N_k} \sum_{\mathbf{k}} \sum_{ij}
   \left(
   [1 - \delta_{0, \omega_n} \delta_{\epsilon_{\mathbf{k},i},\epsilon_{\mathbf{k}+\mathbf{q}, j}})]
   \frac{ f(\epsilon_{\mathbf{k}, i}) - f(\epsilon_{\mathbf{k}+\mathbf{q}, j}) }
   {i\omega_n + \epsilon_{\mathbf{k} + \mathbf{q}, j} - \epsilon_{\mathbf{k}, i}}
   +
   \delta_{0, \omega_n} \delta_{\epsilon_{\mathbf{k},i},\epsilon_{\mathbf{k}+\mathbf{q}, j}}
   \frac{\beta}{4 \cosh^2 (\beta \epsilon_{\mathbf{k}, i} / 2) }
   \right)
   \\ \times
   U_{\bar{a}i}(\mathbf{k}) U^\dagger_{id}(\mathbf{k})
   U_{\bar{c}j}(\mathbf{k} + \mathbf{q}) U^\dagger_{jb}(\mathbf{k} + \mathbf{q})

   where the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued 
   dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

   

   

.. math::

   \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
   = \delta_{ij} \epsilon_{\mathbf{k}, i}

   .. note::
      The analytic formula is sub-optimal in terms of performance for higher temperatures. The evaluation
      scales as :math:`\mathcal{O}(N_k^2)` which is worse than computing the bubble explicitly in imaginary 
      time, with scaling :math:`\mathcal{O}(N_k N_\tau \log(N_k N_\tau)` for :math:`N_k \gg N_\tau`.

   .. note::
      Care must be taken when evaluating the fermionic Matsubara frequency sum of the
      product of two simple poles. By extending the sum to an integral over the complex 
      plane the standard expression for the Lindhard response is obtained when the 
      poles are non-degenerate. The degenerate case produces an additional frequency independent
      contribution (the last term on the last row).

------

[3] Generalized Lindhard susceptibility in the particle-hole channel and for real frequencies :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{q})`.

   Analytic calculation of the generalized (non-interacting) Lindhard susceptibility 
   in the particle-hole channel in real frequencies. The analytic expression is obtained using 
   residue calculus to explicitly evaluate the matsubara sum of the fourier transformed imaginary
   time bubble product of two non-interacting single-particle Green's functions.

   

.. math::

   G^{(0)}_{a\bar{b}}(\mathbf{k}, i\omega_n) =
   \left[ i\omega_n \cdot \mathbf{1} - \epsilon(\mathbf{k}) \right]^{-1} .

   The analytic continuation of the resulting expression to the real frequency axis gives

   

.. math::

   \chi^{(00)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{q}) =
   \frac{1}{N_k} \sum_{\mathbf{k}} \sum_{ij}
   \frac{ f(\epsilon_{\mathbf{k}, i}) - f(\epsilon_{\mathbf{k}+\mathbf{q}, j}) }
   {\omega + i\delta + \epsilon_{\mathbf{k} + \mathbf{q}, j} - \epsilon_{\mathbf{k}, i}}
   \\ \times
   U_{\bar{a}i}(\mathbf{k}) U^\dagger_{id}(\mathbf{k})
   U_{\bar{c}j}(\mathbf{k} + \mathbf{q}) U^\dagger_{jb}(\mathbf{k} + \mathbf{q})

   where the :math:`U(\mathbf{k})` matrices are the diagonalizing unitary transform of the matrix valued 
   dispersion relation :math:`\epsilon_{\bar{a}b}(\mathbf{k})`, i.e.

   

.. math::

   \sum_{\bar{a}b} U_{i\bar{a}}(\mathbf{k}) \epsilon_{\bar{a}b}(\mathbf{k}) U^\dagger_{bj} (\mathbf{k})
   = \delta_{ij} \epsilon_{\mathbf{k}, i}

------

Parameters
----------
e_k : {par_0}
   discretized lattice dispersion :math:`\epsilon_{\bar{a}b}(\mathbf{k})`
mesh : {par_1}
   bosonic Matsubara frequency mesh
mu : {par_2}
   chemical potential :math:`\mu`
beta : {par_3}
   inverse temperature
delta : {par_4}
   broadening :math:`\delta`

Returns
-------
[1] : {ret_0}
   generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(i\omega_n, \mathbf{q})`

[3] : {ret_1}
   real frequency generalized Lindhard susceptibility in the particle-hole channel :math:`\chi^{(00)}_{\bar{a}b\bar{c}d}(\omega, \mathbf{q})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::e_k_cvt>()},
                                                   {c2py::python_typename<triqs::mesh::imfreq>()},
                                                   {c2py::python_typename<double>()},
                                                   {c2py::python_typename<double>()},
                                                   {c2py::python_typename<double>()}},
                                                  {c2py::python_typename<triqs_tprf::chi_wk_t>(), c2py::python_typename<triqs_tprf::chi_fk_t>()});
static const auto _c2py_doc_64 = _c2py_fun_64.doc(R"DOC(
Density matrix from lattic Green's function

Parameters
----------
g_wk : {par_0}
   single particle Green's function :math:`G_{ab}(i\omega_n, \mathbf{k})`

Returns
-------
{ret_0}
   rho_k density matrix :math:`\rho_{ab}(\mathbf{k})`
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::g_wk_cvt>()}}, {c2py::python_typename<triqs_tprf::e_k_t>()});
static const auto _c2py_doc_65 = _c2py_fun_65.doc(
   R"DOC(
Random Phase Approximation (RPA) in the particle-hole channel

    Computes the equation

    

.. math::

   \chi(\bar{a}b\bar{c}d) = \big(
   \mathbb{1}
   - \chi^{(0)}(\bar{a}b\bar{B}A) U(A\bar{B}D\bar{C})
   \big)^{-1} \chi^{(0)}(\bar{C}D\bar{c}d)\,.

Parameters
----------
chi0 : {par_0}
   bare particle-hole bubble :math:`\chi^{(0)}_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`
U : {par_1}
   RPA static vertex as obtained from triqs_tprf.rpa_tensor.get_rpa_tensor :math:`U_{a\bar{b}c\bar{d}}`

Returns
-------
[1] : {ret_0}
   RPA suceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, i\omega_n)`

[3] : {ret_1}
   RPA suceptibility :math:`\chi_{\bar{a}b\bar{c}d}(\mathbf{k}, \omega)`
)DOC",
   {{c2py::python_typename<triqs_tprf::chi_wk_vt>()},
    {c2py::python_typename<
       nda::basic_array_view<std::complex<double>, 4, nda::C_layout, 'A', nda::default_accessor, nda::borrowed<nda::mem::AddressSpace::Host>>>()}},
   {c2py::python_typename<triqs_tprf::chi_wk_t>(), c2py::python_typename<triqs_tprf::chi_fk_t>()});
static const auto _c2py_doc_66 = _c2py_fun_66.doc(R"DOC(
Splits a rank 4 tensor-valued Green's function into dynamic and constant parts by tail fitting

   Splits a general rank 4 tensor-valued Green's function :math:`\chi_{abcd}(i\omega_n, \mathbf{k})` 
   into a dynamic and a constant part in Matsubara frequency space by fitting
   the high-frequency tail.

   .. math ::
       {abcd}(i, {k}) = 
           ^{(dyn)}_{abcd}(i, {k})
           + ^{(stat)}_{abcd}({k})

Parameters
----------
chi_wk : {par_0}
   : general rank 4 tensor-valued Green's function :math:`\chi_{abcd}(i\omega_n, \mathbf{k})`.

Returns
-------
{ret_0}
   Tuple of chi_dyn_wk, the dynamic part of chi :math:`\chi^{(dyn)}_{abcd}(i\omega_n, \mathbf{k})`, which converges to zero for :math:`\omega_n \rightarrow \infty`, and chi_const_k, the part of chi that is constant in Matsubara frequency space :math:`\chi^{(stat)}_{abcd}(\mathbf{k})`.
)DOC",
                                                  {{c2py::python_typename<triqs_tprf::chi_wk_cvt>()}},
                                                  {c2py::python_typename<std::tuple<triqs_tprf::chi_wk_t, triqs_tprf::chi_k_t>>()});
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"add_dynamic_and_static", (PyCFunction)c2py::pyfkw<_c2py_fun_0>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_0.c_str()},
   {"attatch_tri_vert", (PyCFunction)c2py::pyfkw<_c2py_fun_1>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_1.c_str()},
   {"bose", (PyCFunction)c2py::pyfkw<_c2py_fun_2>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_2.c_str()},
   {"chi0_Tr_from_g_Tr_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_3>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_3.c_str()},
   {"chi0_nr_from_gr_PH_at_specific_w", (PyCFunction)c2py::pyfkw<_c2py_fun_4>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_4.c_str()},
   {"chi0_tr_from_grt_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_5>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_5.c_str()},
   {"chi0_w0r_from_grt_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_6>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_6.c_str()},
   {"chi0_wr_from_grt_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_7>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_7.c_str()},
   {"chi0q_from_chi0r", (PyCFunction)c2py::pyfkw<_c2py_fun_8>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_8.c_str()},
   {"chi0q_from_g_wk_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_9>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_9.c_str()},
   {"chi0q_sum_nu", (PyCFunction)c2py::pyfkw<_c2py_fun_10>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_10.c_str()},
   {"chi0q_sum_nu_q", (PyCFunction)c2py::pyfkw<_c2py_fun_11>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_11.c_str()},
   {"chi0q_sum_nu_tail_corr_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_12>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_12.c_str()},
   {"chi0r_from_chi0q", (PyCFunction)c2py::pyfkw<_c2py_fun_13>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_13.c_str()},
   {"chi0r_from_gr_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_14>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_14.c_str()},
   {"chi0r_from_gr_PH_nompi", (PyCFunction)c2py::pyfkw<_c2py_fun_15>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_15.c_str()},
   {"chi_tr_from_chi_wr", (PyCFunction)c2py::pyfkw<_c2py_fun_16>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_16.c_str()},
   {"chi_trapz_tau", (PyCFunction)c2py::pyfkw<_c2py_fun_17>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_17.c_str()},
   {"chi_w0r_from_chi_tr", (PyCFunction)c2py::pyfkw<_c2py_fun_18>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_18.c_str()},
   {"chi_wk_from_chi_wr", (PyCFunction)c2py::pyfkw<_c2py_fun_19>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_19.c_str()},
   {"chi_wr_from_chi_tr", (PyCFunction)c2py::pyfkw<_c2py_fun_20>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_20.c_str()},
   {"chi_wr_from_chi_wk", (PyCFunction)c2py::pyfkw<_c2py_fun_21>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_21.c_str()},
   {"chiq_from_chi0q_and_gamma_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_22>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_22.c_str()},
   {"chiq_sum_nu", (PyCFunction)c2py::pyfkw<_c2py_fun_23>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_23.c_str()},
   {"chiq_sum_nu_from_chi0q_and_gamma_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_24>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_24.c_str()},
   {"chiq_sum_nu_from_chi0q_and_gamma_and_L_wn_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_25>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_25.c_str()},
   {"chiq_sum_nu_from_e_k_sigma_w_and_gamma_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_26>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_26.c_str()},
   {"chiq_sum_nu_from_g_wk_and_gamma_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_27>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_27.c_str()},
   {"chiq_sum_nu_q", (PyCFunction)c2py::pyfkw<_c2py_fun_28>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_28.c_str()},
   {"cluster_mesh_fourier_interpolation", (PyCFunction)c2py::pyfkw<_c2py_fun_29>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_29.c_str()},
   {"construct_phi_wk", (PyCFunction)c2py::pyfkw<_c2py_fun_30>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_30.c_str()},
   {"dlr_on_imfreq", (PyCFunction)c2py::pyfkw<_c2py_fun_31>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_31.c_str()},
   {"dlr_on_imtime", (PyCFunction)c2py::pyfkw<_c2py_fun_32>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_32.c_str()},
   {"dynamic_and_constant_to_tr", (PyCFunction)c2py::pyfkw<_c2py_fun_33>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_33.c_str()},
   {"dynamical_screened_interaction_W", (PyCFunction)c2py::pyfkw<_c2py_fun_34>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_34.c_str()},
   {"dynamical_screened_interaction_W_from_generalized_susceptibility", (PyCFunction)c2py::pyfkw<_c2py_fun_35>, METH_VARARGS | METH_KEYWORDS,
    _c2py_doc_35.c_str()},
   {"eliashberg_constant_gamma_f_product", (PyCFunction)c2py::pyfkw<_c2py_fun_36>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_36.c_str()},
   {"eliashberg_g_delta_g_product", (PyCFunction)c2py::pyfkw<_c2py_fun_37>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_37.c_str()},
   {"eliashberg_product", (PyCFunction)c2py::pyfkw<_c2py_fun_38>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_38.c_str()},
   {"eliashberg_product_fft", (PyCFunction)c2py::pyfkw<_c2py_fun_39>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_39.c_str()},
   {"eliashberg_product_fft_constant", (PyCFunction)c2py::pyfkw<_c2py_fun_40>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_40.c_str()},
   {"fermi", (PyCFunction)c2py::pyfkw<_c2py_fun_41>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_41.c_str()},
   {"fock_sigma", (PyCFunction)c2py::pyfkw<_c2py_fun_42>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_42.c_str()},
   {"fourier_Tk_to_Tr", (PyCFunction)c2py::pyfkw<_c2py_fun_43>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_43.c_str()},
   {"fourier_Tr_to_Tk", (PyCFunction)c2py::pyfkw<_c2py_fun_44>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_44.c_str()},
   {"fourier_fk_to_fr", (PyCFunction)c2py::pyfkw<_c2py_fun_45>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_45.c_str()},
   {"fourier_fr_to_fk", (PyCFunction)c2py::pyfkw<_c2py_fun_46>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_46.c_str()},
   {"fourier_tr_to_wr", (PyCFunction)c2py::pyfkw<_c2py_fun_47>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_47.c_str()},
   {"fourier_wk_to_wr", (PyCFunction)c2py::pyfkw<_c2py_fun_48>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_48.c_str()},
   {"fourier_wr_to_tr", (PyCFunction)c2py::pyfkw<_c2py_fun_49>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_49.c_str()},
   {"fourier_wr_to_wk", (PyCFunction)c2py::pyfkw<_c2py_fun_50>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_50.c_str()},
   {"g0_Tk_les_gtr_from_e_k", (PyCFunction)c2py::pyfkw<_c2py_fun_51>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_51.c_str()},
   {"g0w_dynamic_sigma", (PyCFunction)c2py::pyfkw<_c2py_fun_52>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_52.c_str()},
   {"g0w_sigma", (PyCFunction)c2py::pyfkw<_c2py_fun_53>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_53.c_str()},
   {"gw_dynamic_sigma", (PyCFunction)c2py::pyfkw<_c2py_fun_54>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_54.c_str()},
   {"gw_sigma", (PyCFunction)c2py::pyfkw<_c2py_fun_55>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_55.c_str()},
   {"hartree_sigma", (PyCFunction)c2py::pyfkw<_c2py_fun_56>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_56.c_str()},
   {"lattice_dyson_g0_fk", (PyCFunction)c2py::pyfkw<_c2py_fun_57>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_57.c_str()},
   {"lattice_dyson_g0_wk", (PyCFunction)c2py::pyfkw<_c2py_fun_58>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_58.c_str()},
   {"lattice_dyson_g_f", (PyCFunction)c2py::pyfkw<_c2py_fun_59>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_59.c_str()},
   {"lattice_dyson_g_fk", (PyCFunction)c2py::pyfkw<_c2py_fun_60>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_60.c_str()},
   {"lattice_dyson_g_w", (PyCFunction)c2py::pyfkw<_c2py_fun_61>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_61.c_str()},
   {"lattice_dyson_g_wk", (PyCFunction)c2py::pyfkw<_c2py_fun_62>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_62.c_str()},
   {"lindhard_chi00", (PyCFunction)c2py::pyfkw<_c2py_fun_63>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_63.c_str()},
   {"rho_k_from_g_wk", (PyCFunction)c2py::pyfkw<_c2py_fun_64>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_64.c_str()},
   {"solve_rpa_PH", (PyCFunction)c2py::pyfkw<_c2py_fun_65>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_65.c_str()},
   {"split_into_dynamic_wk_and_constant_k", (PyCFunction)c2py::pyfkw<_c2py_fun_66>, METH_VARARGS | METH_KEYWORDS, _c2py_doc_66.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {PyModuleDef_HEAD_INIT,
                                        "lattice",                              /* name of module */
                                        R"RAWDOC(Lattice functionality)RAWDOC", /* module documentation, may be NULL */
                                        -1, /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
                                        module_methods,
                                        NULL,
                                        NULL,
                                        NULL,
                                        NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_lattice() {

  if (not c2py::check_python_version("lattice")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;
#define _add_type(T, N) c2py::add_type_object_to_main<T>(N, m, conv_table)

#undef _add_type

  return m;
}
#endif
// CLAIR_WRAP_GEN
