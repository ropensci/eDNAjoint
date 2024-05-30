// Generated by rstantools.  Do not edit by hand.

/*
    eDNAjoint is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    eDNAjoint is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with eDNAjoint.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_joint_binary_negbin_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 48> locations_array__ =
  {" (found before start of program)",
  " (in 'string', line 14, column 4 to column 33)",
  " (in 'string', line 15, column 4 to column 22)",
  " (in 'string', line 16, column 4 to column 23)",
  " (in 'string', line 17, column 4 to column 26)",
  " (in 'string', line 20, column 2 to column 43)",
  " (in 'string', line 21, column 2 to column 40)",
  " (in 'string', line 40, column 2 to column 22)",
  " (in 'string', line 41, column 2 to column 11)",
  " (in 'string', line 23, column 4 to column 41)",
  " (in 'string', line 24, column 4 to column 33)",
  " (in 'string', line 22, column 19 to line 25, column 3)",
  " (in 'string', line 22, column 2 to line 25, column 3)",
  " (in 'string', line 42, column 2 to column 21)",
  " (in 'string', line 44, column 10 to column 65)",
  " (in 'string', line 43, column 17 to line 45, column 7)",
  " (in 'string', line 43, column 4 to line 45, column 7)",
  " (in 'string', line 47, column 10 to column 61)",
  " (in 'string', line 46, column 19 to line 48, column 7)",
  " (in 'string', line 46, column 6 to line 48, column 7)",
  " (in 'string', line 29, column 8 to column 45)",
  " (in 'string', line 28, column 18 to line 30, column 5)",
  " (in 'string', line 28, column 4 to line 30, column 5)",
  " (in 'string', line 32, column 8 to column 39)",
  " (in 'string', line 31, column 18 to line 33, column 5)",
  " (in 'string', line 31, column 4 to line 33, column 5)",
  " (in 'string', line 35, column 2 to column 47)",
  " (in 'string', line 36, column 2 to column 22)",
  " (in 'string', line 37, column 2 to column 42)",
  " (in 'string', line 2, column 4 to column 19)",
  " (in 'string', line 3, column 4 to column 19)",
  " (in 'string', line 4, column 10 to column 11)",
  " (in 'string', line 4, column 4 to column 28)",
  " (in 'string', line 5, column 10 to column 11)",
  " (in 'string', line 5, column 4 to column 28)",
  " (in 'string', line 6, column 4 to column 22)",
  " (in 'string', line 7, column 10 to column 11)",
  " (in 'string', line 7, column 4 to column 28)",
  " (in 'string', line 8, column 10 to column 11)",
  " (in 'string', line 8, column 4 to column 28)",
  " (in 'string', line 9, column 10 to column 11)",
  " (in 'string', line 9, column 4 to column 28)",
  " (in 'string', line 10, column 4 to column 28)",
  " (in 'string', line 11, column 4 to column 28)",
  " (in 'string', line 14, column 10 to column 14)",
  " (in 'string', line 20, column 8 to column 12)",
  " (in 'string', line 21, column 8 to column 12)",
  " (in 'string', line 40, column 9 to column 12)"};
#include <stan_meta_header.hpp>
class model_joint_binary_negbin final : public model_base_crtp<model_joint_binary_negbin> {
private:
  int S;
  int C;
  std::vector<int> L;
  std::vector<int> R;
  int Nloc;
  std::vector<int> E;
  std::vector<int> N;
  std::vector<int> K;
  std::vector<double> p10priors;
  std::vector<double> phipriors;
  int log_lik_1dim__;
public:
  ~model_joint_binary_negbin() {}
  model_joint_binary_negbin(stan::io::var_context& context__, unsigned int
                            random_seed__ = 0, std::ostream*
                            pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_joint_binary_negbin_namespace::model_joint_binary_negbin";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 29;
      context__.validate_dims("data initialization", "S", "int",
        std::vector<size_t>{});
      S = std::numeric_limits<int>::min();
      current_statement__ = 29;
      S = context__.vals_i("S")[(1 - 1)];
      current_statement__ = 29;
      stan::math::check_greater_or_equal(function__, "S", S, 1);
      current_statement__ = 30;
      context__.validate_dims("data initialization", "C", "int",
        std::vector<size_t>{});
      C = std::numeric_limits<int>::min();
      current_statement__ = 30;
      C = context__.vals_i("C")[(1 - 1)];
      current_statement__ = 30;
      stan::math::check_greater_or_equal(function__, "C", C, 1);
      current_statement__ = 31;
      stan::math::validate_non_negative_index("L", "S", S);
      current_statement__ = 32;
      context__.validate_dims("data initialization", "L", "int",
        std::vector<size_t>{static_cast<size_t>(S)});
      L = std::vector<int>(S, std::numeric_limits<int>::min());
      current_statement__ = 32;
      L = context__.vals_i("L");
      current_statement__ = 32;
      stan::math::check_greater_or_equal(function__, "L", L, 1);
      current_statement__ = 33;
      stan::math::validate_non_negative_index("R", "C", C);
      current_statement__ = 34;
      context__.validate_dims("data initialization", "R", "int",
        std::vector<size_t>{static_cast<size_t>(C)});
      R = std::vector<int>(C, std::numeric_limits<int>::min());
      current_statement__ = 34;
      R = context__.vals_i("R");
      current_statement__ = 34;
      stan::math::check_greater_or_equal(function__, "R", R, 1);
      current_statement__ = 35;
      context__.validate_dims("data initialization", "Nloc", "int",
        std::vector<size_t>{});
      Nloc = std::numeric_limits<int>::min();
      current_statement__ = 35;
      Nloc = context__.vals_i("Nloc")[(1 - 1)];
      current_statement__ = 35;
      stan::math::check_greater_or_equal(function__, "Nloc", Nloc, 1);
      current_statement__ = 36;
      stan::math::validate_non_negative_index("E", "C", C);
      current_statement__ = 37;
      context__.validate_dims("data initialization", "E", "int",
        std::vector<size_t>{static_cast<size_t>(C)});
      E = std::vector<int>(C, std::numeric_limits<int>::min());
      current_statement__ = 37;
      E = context__.vals_i("E");
      current_statement__ = 37;
      stan::math::check_greater_or_equal(function__, "E", E, 0);
      current_statement__ = 38;
      stan::math::validate_non_negative_index("N", "S", S);
      current_statement__ = 39;
      context__.validate_dims("data initialization", "N", "int",
        std::vector<size_t>{static_cast<size_t>(S)});
      N = std::vector<int>(S, std::numeric_limits<int>::min());
      current_statement__ = 39;
      N = context__.vals_i("N");
      current_statement__ = 39;
      stan::math::check_greater_or_equal(function__, "N", N, 1);
      current_statement__ = 40;
      stan::math::validate_non_negative_index("K", "S", S);
      current_statement__ = 41;
      context__.validate_dims("data initialization", "K", "int",
        std::vector<size_t>{static_cast<size_t>(S)});
      K = std::vector<int>(S, std::numeric_limits<int>::min());
      current_statement__ = 41;
      K = context__.vals_i("K");
      current_statement__ = 41;
      stan::math::check_greater_or_equal(function__, "K", K, 0);
      current_statement__ = 42;
      context__.validate_dims("data initialization", "p10priors", "double",
        std::vector<size_t>{static_cast<size_t>(2)});
      p10priors = std::vector<double>(2,
                    std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 42;
      p10priors = context__.vals_r("p10priors");
      current_statement__ = 43;
      context__.validate_dims("data initialization", "phipriors", "double",
        std::vector<size_t>{static_cast<size_t>(2)});
      phipriors = std::vector<double>(2,
                    std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 43;
      phipriors = context__.vals_r("phipriors");
      current_statement__ = 44;
      stan::math::validate_non_negative_index("mu", "Nloc", Nloc);
      current_statement__ = 45;
      stan::math::validate_non_negative_index("p11", "Nloc", Nloc);
      current_statement__ = 46;
      stan::math::validate_non_negative_index("p", "Nloc", Nloc);
      current_statement__ = 47;
      log_lik_1dim__ = std::numeric_limits<int>::min();
      current_statement__ = 47;
      log_lik_1dim__ = (C + S);
      current_statement__ = 47;
      stan::math::validate_non_negative_index("log_lik", "C + S",
        log_lik_1dim__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = Nloc + 1 + 1 + 1;
  }
  inline std::string model_name() const final {
    return "model_joint_binary_negbin";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_joint_binary_negbin_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      std::vector<local_scalar_t__> mu =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 1;
      mu = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
             jacobian__>(0, lp__, Nloc);
      local_scalar_t__ phi = DUMMY_VAR__;
      current_statement__ = 2;
      phi = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
              lp__);
      local_scalar_t__ beta = DUMMY_VAR__;
      current_statement__ = 3;
      beta = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
               lp__);
      local_scalar_t__ log_p10 = DUMMY_VAR__;
      current_statement__ = 4;
      log_p10 = in__.template read_constrain_ub<local_scalar_t__,
                  jacobian__>(0, lp__);
      std::vector<local_scalar_t__> p11 =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      std::vector<local_scalar_t__> p =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 12;
      for (int i = 1; i <= Nloc; ++i) {
        current_statement__ = 9;
        stan::model::assign(p11,
          (stan::model::rvalue(mu, "mu", stan::model::index_uni(i)) /
          (stan::model::rvalue(mu, "mu", stan::model::index_uni(i)) +
          stan::math::exp(beta))), "assigning variable p11",
          stan::model::index_uni(i));
        current_statement__ = 10;
        stan::model::assign(p,
          (stan::model::rvalue(p11, "p11", stan::model::index_uni(i)) +
          stan::math::exp(log_p10)), "assigning variable p",
          stan::model::index_uni(i));
      }
      current_statement__ = 5;
      stan::math::check_greater_or_equal(function__, "p11", p11, 0);
      current_statement__ = 5;
      stan::math::check_less_or_equal(function__, "p11", p11, 1);
      current_statement__ = 6;
      stan::math::check_greater_or_equal(function__, "p", p, 0);
      current_statement__ = 6;
      stan::math::check_less_or_equal(function__, "p", p, 1);
      {
        current_statement__ = 22;
        for (int j = 1; j <= C; ++j) {
          current_statement__ = 20;
          lp_accum__.add(stan::math::neg_binomial_2_lpmf<propto__>(
                           stan::model::rvalue(E, "E",
                             stan::model::index_uni(j)),
                           stan::model::rvalue(mu, "mu",
                             stan::model::index_uni(
                               stan::model::rvalue(R, "R",
                                 stan::model::index_uni(j)))), phi));
        }
        current_statement__ = 25;
        for (int i = 1; i <= S; ++i) {
          current_statement__ = 23;
          lp_accum__.add(stan::math::binomial_lpmf<propto__>(
                           stan::model::rvalue(K, "K",
                             stan::model::index_uni(i)),
                           stan::model::rvalue(N, "N",
                             stan::model::index_uni(i)),
                           stan::model::rvalue(p, "p",
                             stan::model::index_uni(
                               stan::model::rvalue(L, "L",
                                 stan::model::index_uni(i))))));
        }
        current_statement__ = 26;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(log_p10,
                         stan::model::rvalue(p10priors, "p10priors",
                           stan::model::index_uni(1)),
                         stan::model::rvalue(p10priors, "p10priors",
                           stan::model::index_uni(2))));
        current_statement__ = 27;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(beta, 0, 10));
        current_statement__ = 28;
        lp_accum__.add(stan::math::gamma_lpdf<propto__>(phi,
                         stan::model::rvalue(phipriors, "phipriors",
                           stan::model::index_uni(1)),
                         stan::model::rvalue(phipriors, "phipriors",
                           stan::model::index_uni(2))));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_joint_binary_negbin_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      std::vector<double> mu =
        std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      mu = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
             jacobian__>(0, lp__, Nloc);
      double phi = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      phi = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
              lp__);
      double beta = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 3;
      beta = in__.template read_constrain_lb<local_scalar_t__, jacobian__>(0,
               lp__);
      double log_p10 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 4;
      log_p10 = in__.template read_constrain_ub<local_scalar_t__,
                  jacobian__>(0, lp__);
      std::vector<double> p11 =
        std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      std::vector<double> p =
        std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      out__.write(mu);
      out__.write(phi);
      out__.write(beta);
      out__.write(log_p10);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 12;
      for (int i = 1; i <= Nloc; ++i) {
        current_statement__ = 9;
        stan::model::assign(p11,
          (stan::model::rvalue(mu, "mu", stan::model::index_uni(i)) /
          (stan::model::rvalue(mu, "mu", stan::model::index_uni(i)) +
          stan::math::exp(beta))), "assigning variable p11",
          stan::model::index_uni(i));
        current_statement__ = 10;
        stan::model::assign(p,
          (stan::model::rvalue(p11, "p11", stan::model::index_uni(i)) +
          stan::math::exp(log_p10)), "assigning variable p",
          stan::model::index_uni(i));
      }
      current_statement__ = 5;
      stan::math::check_greater_or_equal(function__, "p11", p11, 0);
      current_statement__ = 5;
      stan::math::check_less_or_equal(function__, "p11", p11, 1);
      current_statement__ = 6;
      stan::math::check_greater_or_equal(function__, "p", p, 0);
      current_statement__ = 6;
      stan::math::check_less_or_equal(function__, "p", p, 1);
      if (emit_transformed_parameters__) {
        out__.write(p11);
        out__.write(p);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
      Eigen::Matrix<double,-1,1> log_lik =
        Eigen::Matrix<double,-1,1>::Constant(log_lik_1dim__,
          std::numeric_limits<double>::quiet_NaN());
      double p10 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 13;
      p10 = stan::math::exp(log_p10);
      current_statement__ = 16;
      for (int j = 1; j <= C; ++j) {
        current_statement__ = 14;
        stan::model::assign(log_lik,
          stan::math::neg_binomial_2_lpmf<false>(
            stan::model::rvalue(E, "E", stan::model::index_uni(j)),
            stan::model::rvalue(mu, "mu",
              stan::model::index_uni(
                stan::model::rvalue(R, "R", stan::model::index_uni(j)))), phi),
          "assigning variable log_lik", stan::model::index_uni(j));
      }
      current_statement__ = 19;
      for (int i = 1; i <= S; ++i) {
        current_statement__ = 17;
        stan::model::assign(log_lik,
          stan::math::binomial_lpmf<false>(
            stan::model::rvalue(K, "K", stan::model::index_uni(i)),
            stan::model::rvalue(N, "N", stan::model::index_uni(i)),
            stan::model::rvalue(p, "p",
              stan::model::index_uni(
                stan::model::rvalue(L, "L", stan::model::index_uni(i))))),
          "assigning variable log_lik", stan::model::index_uni((C + i)));
      }
      out__.write(log_lik);
      out__.write(p10);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      std::vector<local_scalar_t__> mu =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(mu, in__.read<std::vector<local_scalar_t__>>(Nloc),
        "assigning variable mu");
      out__.write_free_lb(0, mu);
      local_scalar_t__ phi = DUMMY_VAR__;
      current_statement__ = 2;
      phi = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, phi);
      local_scalar_t__ beta = DUMMY_VAR__;
      current_statement__ = 3;
      beta = in__.read<local_scalar_t__>();
      out__.write_free_lb(0, beta);
      local_scalar_t__ log_p10 = DUMMY_VAR__;
      current_statement__ = 4;
      log_p10 = in__.read<local_scalar_t__>();
      out__.write_free_ub(0, log_p10);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "mu", "double",
        std::vector<size_t>{static_cast<size_t>(Nloc)});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "phi", "double",
        std::vector<size_t>{});
      current_statement__ = 3;
      context__.validate_dims("parameter initialization", "beta", "double",
        std::vector<size_t>{});
      current_statement__ = 4;
      context__.validate_dims("parameter initialization", "log_p10",
        "double", std::vector<size_t>{});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      std::vector<local_scalar_t__> mu =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 1;
      mu = context__.vals_r("mu");
      out__.write_free_lb(0, mu);
      local_scalar_t__ phi = DUMMY_VAR__;
      current_statement__ = 2;
      phi = context__.vals_r("phi")[(1 - 1)];
      out__.write_free_lb(0, phi);
      local_scalar_t__ beta = DUMMY_VAR__;
      current_statement__ = 3;
      beta = context__.vals_r("beta")[(1 - 1)];
      out__.write_free_lb(0, beta);
      local_scalar_t__ log_p10 = DUMMY_VAR__;
      current_statement__ = 4;
      log_p10 = context__.vals_r("log_p10")[(1 - 1)];
      out__.write_free_ub(0, log_p10);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"mu", "phi", "beta", "log_p10"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"p11", "p"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::string> temp{"log_lik", "p10"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(
                                                                    Nloc)},
                std::vector<size_t>{}, std::vector<size_t>{},
                std::vector<size_t>{}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(Nloc)},
             std::vector<size_t>{static_cast<size_t>(Nloc)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(log_lik_1dim__)},
             std::vector<size_t>{}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
      param_names__.emplace_back(std::string() + "mu" + '.' +
        std::to_string(sym1__));
    }
    param_names__.emplace_back(std::string() + "phi");
    param_names__.emplace_back(std::string() + "beta");
    param_names__.emplace_back(std::string() + "log_p10");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        param_names__.emplace_back(std::string() + "p11" + '.' +
          std::to_string(sym1__));
      }
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        param_names__.emplace_back(std::string() + "p" + '.' +
          std::to_string(sym1__));
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        param_names__.emplace_back(std::string() + "log_lik" + '.' +
          std::to_string(sym1__));
      }
      param_names__.emplace_back(std::string() + "p10");
    }
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
      param_names__.emplace_back(std::string() + "mu" + '.' +
        std::to_string(sym1__));
    }
    param_names__.emplace_back(std::string() + "phi");
    param_names__.emplace_back(std::string() + "beta");
    param_names__.emplace_back(std::string() + "log_p10");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        param_names__.emplace_back(std::string() + "p11" + '.' +
          std::to_string(sym1__));
      }
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        param_names__.emplace_back(std::string() + "p" + '.' +
          std::to_string(sym1__));
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        param_names__.emplace_back(std::string() + "log_lik" + '.' +
          std::to_string(sym1__));
      }
      param_names__.emplace_back(std::string() + "p10");
    }
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"mu\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"phi\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"log_p10\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"p11\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"p\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(log_lik_1dim__) + "},\"block\":\"generated_quantities\"},{\"name\":\"p10\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"mu\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"phi\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"log_p10\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"p11\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"p\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(log_lik_1dim__) + "},\"block\":\"generated_quantities\"},{\"name\":\"p10\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (((Nloc + 1) + 1) + 1);
    const size_t num_transformed = emit_transformed_parameters * ((Nloc +
      Nloc));
    const size_t num_gen_quantities = emit_generated_quantities *
      ((log_lik_1dim__ + 1));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (((Nloc + 1) + 1) + 1);
    const size_t num_transformed = emit_transformed_parameters * ((Nloc +
      Nloc));
    const size_t num_gen_quantities = emit_generated_quantities *
      ((log_lik_1dim__ + 1));
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_joint_binary_negbin_namespace::model_joint_binary_negbin;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_joint_binary_negbin_namespace::profiles__;
}
#endif
#endif
