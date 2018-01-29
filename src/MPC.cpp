#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

// TODO: Set the timestep length and duration
size_t N = 15;
double dt = 0.1;

// Reference values
double ref_v = 35;

// Weight of each factor in cost function
int wgt_cte = 1000;
int wgt_epsi = 1000;
int wgt_ve = 1;
int wgt_delta = 100;
int wgt_delta_chg = 50000;
int wgt_a = 10;
int wgt_a_change = 100;

// Other reusable variables
size_t x_start = 0;
size_t y_start = N;
size_t psi_start = N * 2;
size_t v_start = N * 3;
size_t cte_start = N * 4;
size_t epsi_start = N * 5;
size_t delta_start = (4 + 2) * N;
size_t a_start = (4 + 2) * N + 1 * (N - 1);

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

class FG_eval {
public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;

  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  void operator()(ADvector &fg, const ADvector &vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    //###############################################################
    // Cost function
    //###############################################################
    fg[0] = 0;
    size_t t;
    // Cost function
    // TODO: Define the cost related the reference state and
    // any anything you think may be beneficial.

    // The part of the cost based on the reference state.
    for (t = 0; t < N; t++) {
      fg[0] += wgt_cte * CppAD::pow(vars[cte_start + t], 2);
      fg[0] += wgt_epsi * CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += wgt_ve * CppAD::pow(vars[v_start + t] - ref_v, 2);
    }

    // Minimize the use of actuators.
    for (t = 0; t < N - 1; t++) {
      fg[0] += wgt_delta * CppAD::pow(vars[delta_start + t], 2);
      fg[0] += wgt_a * CppAD::pow(vars[a_start + t], 2);
    }

    // Minimize the value gap between sequential actuations.
    for (t = 0; t < N - 2; t++) {
      fg[0] += wgt_delta_chg * CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += wgt_a_change * CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    //###############################################################
    // Constraints
    //###############################################################
    // At time t=0
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // From time t=1 on
    for (t = 1; t < N; t++) {
      // The state at time t+1
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * pow(x0, 2) + coeffs[3] * pow(x0, 3);
      AD<double> psides0 = CppAD::atan(coeffs[1] + (2 * coeffs[2] * x0) + (3 * coeffs[3] * pow(x0, 2)));

      // Here's `x` to get you started.
      // The idea here is to constraint this value to be 0.
      //
      // Recall the equations for the model:
      // x_[t] = x[t-1] + v[t-1] * cos(psi[t-1]) * dt
      // y_[t] = y[t-1] + v[t-1] * sin(psi[t-1]) * dt
      // psi_[t] = psi[t-1] + v[t-1] / Lf * delta[t-1] * dt
      // v_[t] = v[t-1] + a[t-1] * dt
      // cte[t] = f(x[t-1]) - y[t-1] + v[t-1] * sin(epsi[t-1]) * dt
      // epsi[t] = psi[t] - psides[t-1] + v[t-1] * delta[t-1] / Lf * dt
      fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t] = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t] = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}

MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // TODO: Set the number of model variables (includes both states and inputs).
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  //
  // 4 * 10 + 2 * 9
  // 4+2 states (4 vehicle states, 2 error states) and 2 controls
  size_t n_vars = (4 + 2) * N + 2 * (N - 1);

  // TODO: Set the number of constraints
  // 4+2 constraints between two updates, a total of N-1 updates plus (4 + 2) additional constraints where at time t
  // xt=xt, yt=yt and so on.
  size_t n_constraints = (4 + 2) * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  // Init states
  double x0 = state[0];
  double y0 = state[1];
  double psi0 = state[2];
  double v0 = state[3];
  double cte0 = state[4];
  double epsi0 = state[5];

  vars[x_start] = x0;
  vars[y_start] = y0;
  vars[psi_start] = psi0;
  vars[v_start] = v0;
  vars[cte_start] = cte0;
  vars[epsi_start] = epsi0;

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  // TODO: Set lower and upper limits for variables.
  // For variables with no explicit limit, numerical limits (-inf, +inf) still need to be applied.
  // Otherwise, ipopt will not run properly
  for (i = 0; i < delta_start; i++) {
    vars_upperbound[i] = 1.0e19;
    vars_lowerbound[i] = -1.0e19;
  }

  // Set delta limits to [-25deg, 25deg] or [-0.436332rad, 0.436332rad]
  for (i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Set a limits to [-1, 1]
  for (i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  // Set init state lower and upper limits, i.e. constraint fg[1 + x_start] = vars[x_start], which is equal to x0.
  constraints_lowerbound[x_start] = x0;
  constraints_lowerbound[y_start] = y0;
  constraints_lowerbound[psi_start] = psi0;
  constraints_lowerbound[v_start] = v0;
  constraints_lowerbound[cte_start] = cte0;
  constraints_lowerbound[epsi_start] = epsi0;

  constraints_upperbound[x_start] = x0;
  constraints_upperbound[y_start] = y0;
  constraints_upperbound[psi_start] = psi0;
  constraints_upperbound[v_start] = v0;
  constraints_upperbound[cte_start] = cte0;
  constraints_upperbound[epsi_start] = epsi0;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  // Return first delta and a and predicted x and y
//  vector<double> returned_vals(2 + 2 * N);
  vector<double> returned_vals;

  returned_vals.push_back(solution.x[delta_start]);
  returned_vals.push_back(solution.x[a_start]);

  cout << "optml values: " <<
       "steering_angle:" << solution.x[delta_start] << " " <<
       "throttle:" << solution.x[a_start] << endl;

  // Predicted trajectory
  for (i = 0; i < N - 1; i++) {
    returned_vals.push_back(solution.x[x_start + i + 1]);
    returned_vals.push_back(solution.x[y_start + i + 1]);
  }

  return returned_vals;
}
