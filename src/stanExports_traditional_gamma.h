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
namespace model_traditional_gamma_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 31> locations_array__ =
  {" (found before start of program)",
  " (in 'string', line 8, column 4 to column 36)",
  " (in 'string', line 9, column 4 to column 38)",
  " (in 'string', line 12, column 4 to column 35)",
  " (in 'string', line 25, column 2 to column 20)",
  " (in 'string', line 26, column 2 to column 18)",
  " (in 'string', line 14, column 6 to column 42)",
  " (in 'string', line 13, column 17 to line 15, column 5)",
  " (in 'string', line 13, column 4 to line 15, column 5)",
  " (in 'string', line 28, column 4 to column 29)",
  " (in 'string', line 27, column 18 to line 29, column 3)",
  " (in 'string', line 27, column 2 to line 29, column 3)",
  " (in 'string', line 31, column 10 to column 72)",
  " (in 'string', line 30, column 17 to line 32, column 7)",
  " (in 'string', line 30, column 4 to line 32, column 7)",
  " (in 'string', line 19, column 6 to column 49)",
  " (in 'string', line 18, column 17 to line 20, column 5)",
  " (in 'string', line 18, column 4 to line 20, column 5)",
  " (in 'string', line 21, column 3 to column 27)",
  " (in 'string', line 22, column 3 to column 28)",
  " (in 'string', line 2, column 4 to column 19)",
  " (in 'string', line 3, column 10 to column 11)",
  " (in 'string', line 3, column 4 to column 28)",
  " (in 'string', line 4, column 4 to column 22)",
  " (in 'string', line 5, column 10 to column 11)",
  " (in 'string', line 5, column 4 to column 29)",
  " (in 'string', line 8, column 10 to column 14)",
  " (in 'string', line 9, column 10 to column 14)",
  " (in 'string', line 12, column 10 to column 11)",
  " (in 'string', line 25, column 9 to column 10)",
  " (in 'string', line 26, column 9 to column 13)"};
#include <stan_meta_header.hpp>
class model_traditional_gamma final : public model_base_crtp<model_traditional_gamma> {
private:
  int C;
  std::vector<int> R;
  int Nloc;
  std::vector<double> E;
public:
  ~model_traditional_gamma() {}
  model_traditional_gamma(stan::io::var_context& context__, unsigned int
                          random_seed__ = 0, std::ostream*
                          pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_traditional_gamma_namespace::model_traditional_gamma";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 20;
      context__.validate_dims("data initialization", "C", "int",
        std::vector<size_t>{});
      C = std::numeric_limits<int>::min();
      current_statement__ = 20;
      C = context__.vals_i("C")[(1 - 1)];
      current_statement__ = 20;
      stan::math::check_greater_or_equal(function__, "C", C, 1);
      current_statement__ = 21;
      stan::math::validate_non_negative_index("R", "C", C);
      current_statement__ = 22;
      context__.validate_dims("data initialization", "R", "int",
        std::vector<size_t>{static_cast<size_t>(C)});
      R = std::vector<int>(C, std::numeric_limits<int>::min());
      current_statement__ = 22;
      R = context__.vals_i("R");
      current_statement__ = 22;
      stan::math::check_greater_or_equal(function__, "R", R, 1);
      current_statement__ = 23;
      context__.validate_dims("data initialization", "Nloc", "int",
        std::vector<size_t>{});
      Nloc = std::numeric_limits<int>::min();
      current_statement__ = 23;
      Nloc = context__.vals_i("Nloc")[(1 - 1)];
      current_statement__ = 23;
      stan::math::check_greater_or_equal(function__, "Nloc", Nloc, 1);
      current_statement__ = 24;
      stan::math::validate_non_negative_index("E", "C", C);
      current_statement__ = 25;
      context__.validate_dims("data initialization", "E", "double",
        std::vector<size_t>{static_cast<size_t>(C)});
      E = std::vector<double>(C, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 25;
      E = context__.vals_r("E");
      current_statement__ = 25;
      stan::math::check_greater_or_equal(function__, "E", E, 0);
      current_statement__ = 26;
      stan::math::validate_non_negative_index("alpha", "Nloc", Nloc);
      current_statement__ = 27;
      stan::math::validate_non_negative_index("beta", "Nloc", Nloc);
      current_statement__ = 28;
      stan::math::validate_non_negative_index("E_trans", "C", C);
      current_statement__ = 29;
      stan::math::validate_non_negative_index("log_lik", "C", C);
      current_statement__ = 30;
      stan::math::validate_non_negative_index("mu", "Nloc", Nloc);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = Nloc + Nloc;
  }
  inline std::string model_name() const final {
    return "model_traditional_gamma";
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
      "model_traditional_gamma_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      std::vector<local_scalar_t__> alpha =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 1;
      alpha = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
                jacobian__>(0, lp__, Nloc);
      std::vector<local_scalar_t__> beta =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 2;
      beta = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
               jacobian__>(0.01, lp__, Nloc);
      std::vector<local_scalar_t__> E_trans =
        std::vector<local_scalar_t__>(C, DUMMY_VAR__);
      current_statement__ = 8;
      for (int j = 1; j <= C; ++j) {
        current_statement__ = 6;
        stan::model::assign(E_trans,
          (stan::model::rvalue(E, "E", stan::model::index_uni(j)) +
          0.0000000000001), "assigning variable E_trans",
          stan::model::index_uni(j));
      }
      current_statement__ = 3;
      stan::math::check_greater_or_equal(function__, "E_trans", E_trans, 0);
      {
        current_statement__ = 17;
        for (int j = 1; j <= C; ++j) {
          current_statement__ = 15;
          lp_accum__.add(stan::math::gamma_lpdf<propto__>(
                           stan::model::rvalue(E_trans, "E_trans",
                             stan::model::index_uni(j)),
                           stan::model::rvalue(alpha, "alpha",
                             stan::model::index_uni(
                               stan::model::rvalue(R, "R",
                                 stan::model::index_uni(j)))),
                           stan::model::rvalue(beta, "beta",
                             stan::model::index_uni(
                               stan::model::rvalue(R, "R",
                                 stan::model::index_uni(j))))));
        }
        current_statement__ = 18;
        lp_accum__.add(stan::math::gamma_lpdf<propto__>(beta, 0.25, 0.25));
        current_statement__ = 19;
        lp_accum__.add(stan::math::gamma_lpdf<propto__>(alpha, 0.01, 0.01));
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
      "model_traditional_gamma_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      std::vector<double> alpha =
        std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      alpha = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
                jacobian__>(0, lp__, Nloc);
      std::vector<double> beta =
        std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 2;
      beta = in__.template read_constrain_lb<std::vector<local_scalar_t__>,
               jacobian__>(0.01, lp__, Nloc);
      std::vector<double> E_trans =
        std::vector<double>(C, std::numeric_limits<double>::quiet_NaN());
      out__.write(alpha);
      out__.write(beta);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 8;
      for (int j = 1; j <= C; ++j) {
        current_statement__ = 6;
        stan::model::assign(E_trans,
          (stan::model::rvalue(E, "E", stan::model::index_uni(j)) +
          0.0000000000001), "assigning variable E_trans",
          stan::model::index_uni(j));
      }
      current_statement__ = 3;
      stan::math::check_greater_or_equal(function__, "E_trans", E_trans, 0);
      if (emit_transformed_parameters__) {
        out__.write(E_trans);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
      Eigen::Matrix<double,-1,1> log_lik =
        Eigen::Matrix<double,-1,1>::Constant(C,
          std::numeric_limits<double>::quiet_NaN());
      Eigen::Matrix<double,-1,1> mu =
        Eigen::Matrix<double,-1,1>::Constant(Nloc,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 11;
      for (int j = 1; j <= Nloc; ++j) {
        current_statement__ = 9;
        stan::model::assign(mu,
          (stan::model::rvalue(alpha, "alpha", stan::model::index_uni(j)) /
          stan::model::rvalue(beta, "beta", stan::model::index_uni(j))),
          "assigning variable mu", stan::model::index_uni(j));
      }
      current_statement__ = 14;
      for (int j = 1; j <= C; ++j) {
        current_statement__ = 12;
        stan::model::assign(log_lik,
          stan::math::gamma_lpdf<false>(
            stan::model::rvalue(E_trans, "E_trans", stan::model::index_uni(j)),
            stan::model::rvalue(alpha, "alpha",
              stan::model::index_uni(
                stan::model::rvalue(R, "R", stan::model::index_uni(j)))),
            stan::model::rvalue(beta, "beta",
              stan::model::index_uni(
                stan::model::rvalue(R, "R", stan::model::index_uni(j))))),
          "assigning variable log_lik", stan::model::index_uni(j));
      }
      out__.write(log_lik);
      out__.write(mu);
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
      std::vector<local_scalar_t__> alpha =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(alpha,
        in__.read<std::vector<local_scalar_t__>>(Nloc),
        "assigning variable alpha");
      out__.write_free_lb(0, alpha);
      std::vector<local_scalar_t__> beta =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 2;
      stan::model::assign(beta,
        in__.read<std::vector<local_scalar_t__>>(Nloc),
        "assigning variable beta");
      out__.write_free_lb(0.01, beta);
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
      context__.validate_dims("parameter initialization", "alpha", "double",
        std::vector<size_t>{static_cast<size_t>(Nloc)});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "beta", "double",
        std::vector<size_t>{static_cast<size_t>(Nloc)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      std::vector<local_scalar_t__> alpha =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 1;
      alpha = context__.vals_r("alpha");
      out__.write_free_lb(0, alpha);
      std::vector<local_scalar_t__> beta =
        std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      current_statement__ = 2;
      beta = context__.vals_r("beta");
      out__.write_free_lb(0.01, beta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"alpha", "beta"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"E_trans"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::string> temp{"log_lik", "mu"};
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
                std::vector<size_t>{static_cast<size_t>(Nloc)}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(C)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(C)},
             std::vector<size_t>{static_cast<size_t>(Nloc)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
      param_names__.emplace_back(std::string() + "alpha" + '.' +
        std::to_string(sym1__));
    }
    for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
      param_names__.emplace_back(std::string() + "beta" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        param_names__.emplace_back(std::string() + "E_trans" + '.' +
          std::to_string(sym1__));
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        param_names__.emplace_back(std::string() + "log_lik" + '.' +
          std::to_string(sym1__));
      }
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        param_names__.emplace_back(std::string() + "mu" + '.' +
          std::to_string(sym1__));
      }
    }
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
      param_names__.emplace_back(std::string() + "alpha" + '.' +
        std::to_string(sym1__));
    }
    for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
      param_names__.emplace_back(std::string() + "beta" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        param_names__.emplace_back(std::string() + "E_trans" + '.' +
          std::to_string(sym1__));
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        param_names__.emplace_back(std::string() + "log_lik" + '.' +
          std::to_string(sym1__));
      }
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        param_names__.emplace_back(std::string() + "mu" + '.' +
          std::to_string(sym1__));
      }
    }
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"alpha\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"E_trans\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(C) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(C) + "},\"block\":\"generated_quantities\"},{\"name\":\"mu\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Nloc) + "},\"block\":\"generated_quantities\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"alpha\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nloc) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"E_trans\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(C) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(C) + "},\"block\":\"generated_quantities\"},{\"name\":\"mu\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(Nloc) + "},\"block\":\"generated_quantities\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = (Nloc + Nloc);
    const size_t num_transformed = emit_transformed_parameters * (C);
    const size_t num_gen_quantities = emit_generated_quantities * ((C +
      Nloc));
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
    const size_t num_params__ = (Nloc + Nloc);
    const size_t num_transformed = emit_transformed_parameters * (C);
    const size_t num_gen_quantities = emit_generated_quantities * ((C +
      Nloc));
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
using stan_model = model_traditional_gamma_namespace::model_traditional_gamma;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_traditional_gamma_namespace::profiles__;
}
#endif
#endif
