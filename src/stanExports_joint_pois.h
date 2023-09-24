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
// Code generated by stanc v2.26.1-4-gd72b68b7-dirty
#include <stan/model/model_header.hpp>
namespace model_joint_pois_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'joint_pois', line 13, column 4 to column 33)",
                                                      " (in 'joint_pois', line 14, column 4 to column 23)",
                                                      " (in 'joint_pois', line 15, column 4 to column 26)",
                                                      " (in 'joint_pois', line 18, column 2 to column 43)",
                                                      " (in 'joint_pois', line 19, column 2 to column 41)",
                                                      " (in 'joint_pois', line 21, column 4 to column 41)",
                                                      " (in 'joint_pois', line 22, column 4 to column 33)",
                                                      " (in 'joint_pois', line 20, column 19 to line 23, column 3)",
                                                      " (in 'joint_pois', line 20, column 2 to line 23, column 3)",
                                                      " (in 'joint_pois', line 37, column 2 to column 22)",
                                                      " (in 'joint_pois', line 38, column 2 to column 11)",
                                                      " (in 'joint_pois', line 39, column 2 to column 21)",
                                                      " (in 'joint_pois', line 41, column 10 to column 53)",
                                                      " (in 'joint_pois', line 40, column 17 to line 42, column 7)",
                                                      " (in 'joint_pois', line 40, column 4 to line 42, column 7)",
                                                      " (in 'joint_pois', line 44, column 10 to column 61)",
                                                      " (in 'joint_pois', line 43, column 19 to line 45, column 7)",
                                                      " (in 'joint_pois', line 43, column 6 to line 45, column 7)",
                                                      " (in 'joint_pois', line 27, column 8 to column 33)",
                                                      " (in 'joint_pois', line 26, column 17 to line 28, column 5)",
                                                      " (in 'joint_pois', line 26, column 4 to line 28, column 5)",
                                                      " (in 'joint_pois', line 30, column 8 to column 39)",
                                                      " (in 'joint_pois', line 29, column 18 to line 31, column 5)",
                                                      " (in 'joint_pois', line 29, column 4 to line 31, column 5)",
                                                      " (in 'joint_pois', line 33, column 2 to column 47)",
                                                      " (in 'joint_pois', line 34, column 2 to column 22)",
                                                      " (in 'joint_pois', line 2, column 4 to column 19)",
                                                      " (in 'joint_pois', line 3, column 4 to column 19)",
                                                      " (in 'joint_pois', line 4, column 10 to column 11)",
                                                      " (in 'joint_pois', line 4, column 4 to column 28)",
                                                      " (in 'joint_pois', line 5, column 10 to column 11)",
                                                      " (in 'joint_pois', line 5, column 4 to column 28)",
                                                      " (in 'joint_pois', line 6, column 4 to column 22)",
                                                      " (in 'joint_pois', line 7, column 10 to column 11)",
                                                      " (in 'joint_pois', line 7, column 4 to column 28)",
                                                      " (in 'joint_pois', line 8, column 10 to column 11)",
                                                      " (in 'joint_pois', line 8, column 4 to column 28)",
                                                      " (in 'joint_pois', line 9, column 10 to column 11)",
                                                      " (in 'joint_pois', line 9, column 4 to column 28)",
                                                      " (in 'joint_pois', line 10, column 4 to column 28)",
                                                      " (in 'joint_pois', line 13, column 10 to column 14)",
                                                      " (in 'joint_pois', line 18, column 8 to column 12)",
                                                      " (in 'joint_pois', line 19, column 8 to column 12)",
                                                      " (in 'joint_pois', line 37, column 9 to column 12)"};
#include <stan_meta_header.hpp>
class model_joint_pois final : public model_base_crtp<model_joint_pois> {
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
  int log_lik_1dim__;
 
public:
  ~model_joint_pois() { }
  
  inline std::string model_name() const final { return "model_joint_pois"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1-4-gd72b68b7-dirty", "stancflags = "};
  }
  
  
  model_joint_pois(stan::io::var_context& context__,
                   unsigned int random_seed__ = 0,
                   std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_joint_pois_namespace::model_joint_pois";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 27;
      context__.validate_dims("data initialization","S","int",
          context__.to_vec());
      S = std::numeric_limits<int>::min();
      
      current_statement__ = 27;
      S = context__.vals_i("S")[(1 - 1)];
      current_statement__ = 27;
      current_statement__ = 27;
      check_greater_or_equal(function__, "S", S, 1);
      current_statement__ = 28;
      context__.validate_dims("data initialization","C","int",
          context__.to_vec());
      C = std::numeric_limits<int>::min();
      
      current_statement__ = 28;
      C = context__.vals_i("C")[(1 - 1)];
      current_statement__ = 28;
      current_statement__ = 28;
      check_greater_or_equal(function__, "C", C, 1);
      current_statement__ = 29;
      validate_non_negative_index("L", "S", S);
      current_statement__ = 30;
      context__.validate_dims("data initialization","L","int",
          context__.to_vec(S));
      L = std::vector<int>(S, std::numeric_limits<int>::min());
      
      current_statement__ = 30;
      assign(L, nil_index_list(), context__.vals_i("L"),
        "assigning variable L");
      current_statement__ = 30;
      for (int sym1__ = 1; sym1__ <= S; ++sym1__) {
        current_statement__ = 30;
        current_statement__ = 30;
        check_greater_or_equal(function__, "L[sym1__]", L[(sym1__ - 1)], 1);}
      current_statement__ = 31;
      validate_non_negative_index("R", "C", C);
      current_statement__ = 32;
      context__.validate_dims("data initialization","R","int",
          context__.to_vec(C));
      R = std::vector<int>(C, std::numeric_limits<int>::min());
      
      current_statement__ = 32;
      assign(R, nil_index_list(), context__.vals_i("R"),
        "assigning variable R");
      current_statement__ = 32;
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        current_statement__ = 32;
        current_statement__ = 32;
        check_greater_or_equal(function__, "R[sym1__]", R[(sym1__ - 1)], 1);}
      current_statement__ = 33;
      context__.validate_dims("data initialization","Nloc","int",
          context__.to_vec());
      Nloc = std::numeric_limits<int>::min();
      
      current_statement__ = 33;
      Nloc = context__.vals_i("Nloc")[(1 - 1)];
      current_statement__ = 33;
      current_statement__ = 33;
      check_greater_or_equal(function__, "Nloc", Nloc, 1);
      current_statement__ = 34;
      validate_non_negative_index("E", "C", C);
      current_statement__ = 35;
      context__.validate_dims("data initialization","E","int",
          context__.to_vec(C));
      E = std::vector<int>(C, std::numeric_limits<int>::min());
      
      current_statement__ = 35;
      assign(E, nil_index_list(), context__.vals_i("E"),
        "assigning variable E");
      current_statement__ = 35;
      for (int sym1__ = 1; sym1__ <= C; ++sym1__) {
        current_statement__ = 35;
        current_statement__ = 35;
        check_greater_or_equal(function__, "E[sym1__]", E[(sym1__ - 1)], 0);}
      current_statement__ = 36;
      validate_non_negative_index("N", "S", S);
      current_statement__ = 37;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec(S));
      N = std::vector<int>(S, std::numeric_limits<int>::min());
      
      current_statement__ = 37;
      assign(N, nil_index_list(), context__.vals_i("N"),
        "assigning variable N");
      current_statement__ = 37;
      for (int sym1__ = 1; sym1__ <= S; ++sym1__) {
        current_statement__ = 37;
        current_statement__ = 37;
        check_greater_or_equal(function__, "N[sym1__]", N[(sym1__ - 1)], 1);}
      current_statement__ = 38;
      validate_non_negative_index("K", "S", S);
      current_statement__ = 39;
      context__.validate_dims("data initialization","K","int",
          context__.to_vec(S));
      K = std::vector<int>(S, std::numeric_limits<int>::min());
      
      current_statement__ = 39;
      assign(K, nil_index_list(), context__.vals_i("K"),
        "assigning variable K");
      current_statement__ = 39;
      for (int sym1__ = 1; sym1__ <= S; ++sym1__) {
        current_statement__ = 39;
        current_statement__ = 39;
        check_greater_or_equal(function__, "K[sym1__]", K[(sym1__ - 1)], 0);}
      current_statement__ = 40;
      context__.validate_dims("data initialization","p10priors","double",
          context__.to_vec(2));
      p10priors = std::vector<double>(2, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 40;
      assign(p10priors, nil_index_list(), context__.vals_r("p10priors"),
        "assigning variable p10priors");
      current_statement__ = 41;
      validate_non_negative_index("mu", "Nloc", Nloc);
      current_statement__ = 42;
      validate_non_negative_index("p11", "Nloc", Nloc);
      current_statement__ = 43;
      validate_non_negative_index("p", "Nloc", Nloc);
      current_statement__ = 44;
      log_lik_1dim__ = std::numeric_limits<int>::min();
      
      current_statement__ = 44;
      log_lik_1dim__ = (C + S);
      current_statement__ = 44;
      validate_non_negative_index("log_lik", "C + S", log_lik_1dim__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += Nloc;
      num_params_r__ += 1;
      num_params_r__ += 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_joint_pois_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      std::vector<local_scalar_t__> mu;
      mu = std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 1;
        assign(mu, cons_list(index_uni(sym1__), nil_index_list()),
          in__.scalar(), "assigning variable mu");}
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 1;
        if (jacobian__) {
          current_statement__ = 1;
          assign(mu, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(mu[(sym1__ - 1)], 0, lp__),
            "assigning variable mu");
        } else {
          current_statement__ = 1;
          assign(mu, cons_list(index_uni(sym1__), nil_index_list()),
            stan::math::lb_constrain(mu[(sym1__ - 1)], 0),
            "assigning variable mu");
        }}
      local_scalar_t__ beta;
      beta = DUMMY_VAR__;
      
      current_statement__ = 2;
      beta = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        beta = stan::math::lb_constrain(beta, 0, lp__);
      } else {
        current_statement__ = 2;
        beta = stan::math::lb_constrain(beta, 0);
      }
      local_scalar_t__ log_p10;
      log_p10 = DUMMY_VAR__;
      
      current_statement__ = 3;
      log_p10 = in__.scalar();
      current_statement__ = 3;
      if (jacobian__) {
        current_statement__ = 3;
        log_p10 = stan::math::ub_constrain(log_p10, 0, lp__);
      } else {
        current_statement__ = 3;
        log_p10 = stan::math::ub_constrain(log_p10, 0);
      }
      std::vector<local_scalar_t__> p11;
      p11 = std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      
      std::vector<local_scalar_t__> p;
      p = std::vector<local_scalar_t__>(Nloc, DUMMY_VAR__);
      
      current_statement__ = 9;
      for (int i = 1; i <= Nloc; ++i) {
        current_statement__ = 6;
        assign(p11, cons_list(index_uni(i), nil_index_list()),
          (mu[(i - 1)] / (mu[(i - 1)] + stan::math::exp(beta))),
          "assigning variable p11");
        current_statement__ = 7;
        assign(p, cons_list(index_uni(i), nil_index_list()),
          (p11[(i - 1)] + stan::math::exp(log_p10)), "assigning variable p");
      }
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 4;
        current_statement__ = 4;
        check_greater_or_equal(function__, "p11[sym1__]", p11[(sym1__ - 1)],
                               0);}
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 4;
        current_statement__ = 4;
        check_less_or_equal(function__, "p11[sym1__]", p11[(sym1__ - 1)], 1);
      }
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 5;
        current_statement__ = 5;
        check_greater_or_equal(function__, "p[sym1__]", p[(sym1__ - 1)], 0);}
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 5;
        current_statement__ = 5;
        check_less_or_equal(function__, "p[sym1__]", p[(sym1__ - 1)], 1);}
      {
        current_statement__ = 21;
        for (int j = 1; j <= C; ++j) {
          current_statement__ = 19;
          lp_accum__.add(
            poisson_lpmf<propto__>(E[(j - 1)], mu[(R[(j - 1)] - 1)]));}
        current_statement__ = 24;
        for (int i = 1; i <= S; ++i) {
          current_statement__ = 22;
          lp_accum__.add(
            binomial_lpmf<propto__>(K[(i - 1)], N[(i - 1)],
              p[(L[(i - 1)] - 1)]));}
        current_statement__ = 25;
        lp_accum__.add(
          normal_lpdf<propto__>(log_p10, p10priors[(1 - 1)],
            p10priors[(2 - 1)]));
        current_statement__ = 26;
        lp_accum__.add(normal_lpdf<propto__>(beta, 0, 10));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_joint_pois_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      std::vector<double> mu;
      mu = std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 1;
        assign(mu, cons_list(index_uni(sym1__), nil_index_list()),
          in__.scalar(), "assigning variable mu");}
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 1;
        assign(mu, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_constrain(mu[(sym1__ - 1)], 0),
          "assigning variable mu");}
      double beta;
      beta = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      beta = in__.scalar();
      current_statement__ = 2;
      beta = stan::math::lb_constrain(beta, 0);
      double log_p10;
      log_p10 = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      log_p10 = in__.scalar();
      current_statement__ = 3;
      log_p10 = stan::math::ub_constrain(log_p10, 0);
      std::vector<double> p11;
      p11 = std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      
      std::vector<double> p;
      p = std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        vars__.emplace_back(mu[(sym1__ - 1)]);}
      vars__.emplace_back(beta);
      vars__.emplace_back(log_p10);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 9;
      for (int i = 1; i <= Nloc; ++i) {
        current_statement__ = 6;
        assign(p11, cons_list(index_uni(i), nil_index_list()),
          (mu[(i - 1)] / (mu[(i - 1)] + stan::math::exp(beta))),
          "assigning variable p11");
        current_statement__ = 7;
        assign(p, cons_list(index_uni(i), nil_index_list()),
          (p11[(i - 1)] + stan::math::exp(log_p10)), "assigning variable p");
      }
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 4;
        current_statement__ = 4;
        check_greater_or_equal(function__, "p11[sym1__]", p11[(sym1__ - 1)],
                               0);}
      current_statement__ = 4;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 4;
        current_statement__ = 4;
        check_less_or_equal(function__, "p11[sym1__]", p11[(sym1__ - 1)], 1);
      }
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 5;
        current_statement__ = 5;
        check_greater_or_equal(function__, "p[sym1__]", p[(sym1__ - 1)], 0);}
      current_statement__ = 5;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 5;
        current_statement__ = 5;
        check_less_or_equal(function__, "p[sym1__]", p[(sym1__ - 1)], 1);}
      if (emit_transformed_parameters__) {
        for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
          vars__.emplace_back(p11[(sym1__ - 1)]);}
        for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
          vars__.emplace_back(p[(sym1__ - 1)]);}
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
      Eigen::Matrix<double, -1, 1> log_lik;
      log_lik = Eigen::Matrix<double, -1, 1>(log_lik_1dim__);
      stan::math::fill(log_lik, std::numeric_limits<double>::quiet_NaN());
      
      double p10;
      p10 = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 12;
      p10 = stan::math::exp(log_p10);
      current_statement__ = 15;
      for (int j = 1; j <= C; ++j) {
        current_statement__ = 13;
        assign(log_lik, cons_list(index_uni(j), nil_index_list()),
          poisson_lpmf<false>(E[(j - 1)], mu[(R[(j - 1)] - 1)]),
          "assigning variable log_lik");}
      current_statement__ = 18;
      for (int i = 1; i <= S; ++i) {
        current_statement__ = 16;
        assign(log_lik, cons_list(index_uni((C + i)), nil_index_list()),
          binomial_lpmf<false>(K[(i - 1)], N[(i - 1)], p[(L[(i - 1)] - 1)]),
          "assigning variable log_lik");}
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        vars__.emplace_back(log_lik[(sym1__ - 1)]);}
      vars__.emplace_back(p10);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      std::vector<double> mu;
      mu = std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 1;
      assign(mu, nil_index_list(), context__.vals_r("mu"),
        "assigning variable mu");
      std::vector<double> mu_free__;
      mu_free__ = std::vector<double>(Nloc, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 1;
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        current_statement__ = 1;
        assign(mu_free__, cons_list(index_uni(sym1__), nil_index_list()),
          stan::math::lb_free(mu[(sym1__ - 1)], 0),
          "assigning variable mu_free__");}
      double beta;
      beta = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      beta = context__.vals_r("beta")[(1 - 1)];
      double beta_free__;
      beta_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      beta_free__ = stan::math::lb_free(beta, 0);
      double log_p10;
      log_p10 = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      log_p10 = context__.vals_r("log_p10")[(1 - 1)];
      double log_p10_free__;
      log_p10_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      log_p10_free__ = stan::math::ub_free(log_p10, 0);
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        vars__.emplace_back(mu_free__[(sym1__ - 1)]);}
      vars__.emplace_back(beta_free__);
      vars__.emplace_back(log_p10_free__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("mu");
    names__.emplace_back("beta");
    names__.emplace_back("log_p10");
    names__.emplace_back("p11");
    names__.emplace_back("p");
    names__.emplace_back("log_lik");
    names__.emplace_back("p10");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(Nloc)});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(Nloc)});
    
    dimss__.emplace_back(std::vector<size_t>{static_cast<size_t>(Nloc)});
    
    dimss__.emplace_back(std::vector<size_t>{
                                             static_cast<size_t>(log_lik_1dim__)
                                             });
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "mu" + '.' + std::to_string(sym1__));
      }}
    param_names__.emplace_back(std::string() + "beta");
    param_names__.emplace_back(std::string() + "log_p10");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "p11" + '.' + std::to_string(sym1__));
        }}
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "p" + '.' + std::to_string(sym1__));
        }}
    }
    
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "log_lik" + '.' + std::to_string(sym1__));
        }}
      param_names__.emplace_back(std::string() + "p10");
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
      {
        param_names__.emplace_back(std::string() + "mu" + '.' + std::to_string(sym1__));
      }}
    param_names__.emplace_back(std::string() + "beta");
    param_names__.emplace_back(std::string() + "log_p10");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "p11" + '.' + std::to_string(sym1__));
        }}
      for (int sym1__ = 1; sym1__ <= Nloc; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "p" + '.' + std::to_string(sym1__));
        }}
    }
    
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "log_lik" + '.' + std::to_string(sym1__));
        }}
      param_names__.emplace_back(std::string() + "p10");
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"mu\",\"type\":{\"name\":\"array\",\"length\":" << Nloc << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"log_p10\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"p11\",\"type\":{\"name\":\"array\",\"length\":" << Nloc << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"p\",\"type\":{\"name\":\"array\",\"length\":" << Nloc << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" << log_lik_1dim__ << "},\"block\":\"generated_quantities\"},{\"name\":\"p10\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"mu\",\"type\":{\"name\":\"array\",\"length\":" << Nloc << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"},{\"name\":\"beta\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"log_p10\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"p11\",\"type\":{\"name\":\"array\",\"length\":" << Nloc << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"p\",\"type\":{\"name\":\"array\",\"length\":" << Nloc << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" << log_lik_1dim__ << "},\"block\":\"generated_quantities\"},{\"name\":\"p10\",\"type\":{\"name\":\"real\"},\"block\":\"generated_quantities\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_joint_pois_namespace::model_joint_pois;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_joint_pois_namespace::profiles__;
}
#endif
#endif
