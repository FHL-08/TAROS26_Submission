clear
close all
clc
import casadi.*

t0 = 0; % start time
tf = 40; % end time
dt = 0.1; % step size
k = t0:dt:tf-dt;
ACTIVATE_CONSENSUS = 'active';

% NMPC finite horizon lengths
N = 10; % trajectory following horizon length
T = length(k);

U_0 = 10; % Source Voltage
R = 0.5; % Armature Resistance
R_z = 2e-4; % Internal Resistance of voltage source
J_m = 1e-6; % Inertia of Motor Shaft
L = 100e-3; % Armature Inductance
N_g = 1; % Gearbox Transmission Ratio
l_L = 0.125; % Distance between left wheel and centre of wheel axle
l_R = l_L; % Distance between right wheel and centre of wheel axle
b_val = 0.5*(l_L + l_R);
r = 0.05; % wheel d_safe
r_G = r/N_g;
m_p = 4.385;
m_w = 0.0575;
m = m_p + 2*m_w; % Mass of robot
l_T = 0.11; % distance between centre of wheel shaft and centre of gravity
J_zp = 0.05; % Intertia of the platform about the centroid
J_zw = 0.002; % Intertia of the wheels about the centroid
J_B = J_zp + 2*J_zw + m_w*(2*l_T^2 + l_L^2 + l_R^2); % Intertia about the centre of rotation
K_t = 5.305e-2; % Torque Constant
K_e = K_t; % Torque Constant
k_r = 4e-5; % Damping coefficient of Motor Shaft
k_v = 0.25; % Damping coefficient of linear motion
k_w = 0.15; % Dampin coefficient of rotational motion

% Robot geometric parameters
X_val = m*b_val*l_T;
Y_val = J_B + m*l_T^2;
sigma_val = norm([X_val; Y_val]);
alpha_val = atan2(Y_val, X_val);

alpha = (r_G^2)/(l_L + l_R);
scaled_radius = (alpha/r_G);
a_L = k_r + k_v*l_R*alpha;
a_R = k_r + k_v*l_L*alpha;
c_L = k_r*l_L + k_w*alpha;
c_R = k_r*l_R + k_w*alpha;

b_L = J_m + m*l_R*alpha;
b_R = J_m + m*l_L*alpha;
d_L = J_m*l_L + J_B*alpha;
d_R = J_m*l_R + J_B*alpha;

RL = -(1/L)*[R+R_z R_z;
             R_z R+R_z];
KL = -(K_e/L)*eye(2);
ac_mat = [a_L a_R;
          -c_L c_R];
bd_mat = [d_R -b_R;
          d_L b_L];
LR = [1 1;
      -l_L l_R];
beta = 1/det(bd_mat);

%%% Matrix Definition for low-level controller
A = [RL KL;
     (K_t*beta)*bd_mat*LR -beta*bd_mat*ac_mat];
B = (U_0/L)*eye(4, 2);
C = r_G*[zeros(2) (LR'\eye(2))];
electr_sys = ss(A, B, C, zeros(2));

%%% System Parameters
x = MX.sym('x');
y = MX.sym('y');
theta = MX.sym('theta');

omega = MX.sym('omega');
vel = MX.sym('vel', 2);
Phi = MX.sym('Phi', 2);
tau = MX.sym('tau', 2);

dim_vel = length(Phi);
dim_tau = length(tau);

%%% Robot Kinematics part 1
q = [x; y; theta]; % states
dim_state = length(q); % dimension of state vector

R_theta = [cos(q(3)), -sin(q(3));
           sin(q(3)),  cos(q(3))];
R_theta_func = Function('R_theta_func', {theta}, {R_theta});

R_q = [R_theta*diag([1, l_T]);
       0 1];
R_q_func = Function('R_q_func', {q}, {R_q});
dq = R_q*vel;
dq_func = Function('dq_func', {q, vel}, {dq});

% Nonlinear solver. Produces q(k + 1) given q(k) and u
q_ode = struct('x', q, 'p', vel, 'ode', dq);
q_solver = integrator('q_solver', 'cvodes', q_ode, t0, dt);

% Defines the ith standard basis vector e
function e = std_basis(n, i)
    e = zeros(n,1);
    e(i) = 1;
end

%%% Robot Kinematics part 2
M_bD = [l_R l_L; l_T -l_T];
bar_M_bD = [l_T l_R; l_T -l_L];
S_y_q = scaled_radius*R_theta*M_bD;
inv_S_y_q = (1/(r*l_T))*bar_M_bD*(R_theta');
S_q = [S_y_q; scaled_radius*[1 -1]];

% convert wheel velocities to robot cartesian velocities
S_y_q_func = Function('S_y_q_func', {theta}, {S_y_q});
% inverse transformation: converts cartesian velocities to wheel velocities
inv_S_y_q_func = Function('inv_S_y_q_func', {theta}, {inv_S_y_q});

skew_mat = [0 -1;
            1  0];

% angular velocity
omega_q = scaled_radius*[1 -1]*Phi;

% skew symmetric matrix S(w)
S_w = omega_q*skew_mat;
S_omega = omega*skew_mat;

dS_y_q = S_w*S_y_q;
S_w_func = Function('S_w_func', {Phi, theta}, {S_w});
S_omega_func = Function('S_omega_func', {omega}, {S_omega});

%%% Defintion of unprojected system matrices
% M_q = diag([m; m; J_B]);
% C_q = zeros(3, 1);
% B_q = (N_g/r_G)*[cos(q(3)), cos(q(3));
%                  sin(q(3)), sin(q(3));
%                  l_R      , -l_L];

det_S_y = -(scaled_radius^2)*(l_R + l_L)*(l_T);
bar_M_q = m*(scaled_radius^2)*[(l_R^2 + l_T^2), (l_L*l_R - l_T^2);
                             (l_R*l_L - l_T^2), (l_L^2 + l_T^2)]...
          + J_B*(scaled_radius^2)*[1, -1; -1, 1];
inv_bar_M_q = bar_M_q\eye(2);
bar_C_q = m*det_S_y*S_w;
bar_B_q_coeff = N_g;
bar_B_q = bar_B_q_coeff*eye(dim_tau);

dphi = inv_bar_M_q*(-bar_C_q*Phi + bar_B_q*tau);
dphi_func = Function('dphi_func', {Phi, tau}, {dphi});

x_state = [q; Phi];
dx_state = [S_q*Phi; dphi];
global_robot_acceleration = dS_y_q*Phi + S_y_q*dphi;

% Nonlinear solver. Produces x_state(k + 1) given x_state(k) and tau
state_ode = struct('x', x_state, 'p', tau, 'ode', dx_state);
state_solver = integrator('state_solver', 'cvodes', state_ode, t0, dt);
dx_state_func = Function('dx_state_func', {x_state, tau}, {dx_state});
global_robot_acceleration_func = Function('global_robot_acceleration_func',...
                                            {x_state, tau}, {global_robot_acceleration});

%%% Define Saturation limits
volts_ub = U_0; % Volts
volts_lb = -U_0; % Volts
forward_vel = -C*(A\B)*[volts_ub; volts_ub];
wheel_vel_ub = full(inv_S_y_q_func(0)*forward_vel);
wheel_vel_lb = -wheel_vel_ub;
max_torque = 0.25;
torque_ub = max_torque*ones(2, 1);
torque_lb = -max_torque*ones(2, 1);

%%% Define Communication Network
mrs_size = 3; % number of robots
adjacency = [0 1 1;
             1 0 1;
             1 1 0]; % adjacency matrix (consensus)
vertex_degrees = sum(adjacency, 2);
leader_neighbour = [1 1 1]; % leader neighbourhood vector (consensus)
switch ACTIVATE_CONSENSUS
    case 'inactive'
        adjacency = zeros(mrs_size);  % adjacency matrix (leader-follower)
        vertex_degrees = sum(adjacency, 2); 
        leader_neighbour = ones(1, mrs_size); % leader neighbourhood (leader-follower)
    case 'active'
        % do nothing
end
degree = diag(vertex_degrees);
laplacian = degree - adjacency;
leader_adjacency = diag(leader_neighbour);
modified_laplacian = laplacian + leader_adjacency;
g_dis = modified_laplacian + adjacency;

%%% trajectory error variables
leader_state = MX.sym('leader_state', dim_state);
dleader_state = MX.sym('dleader_state', dim_state);
ddleader_state = MX.sym('ddleader_state', dim_state);
g_i_dis_var = MX.sym('g_i_dis_var');
basis_var = MX.sym('xi_var', mrs_size);
xi_var = MX.sym('xi_var', dim_state-1);
eta_i = MX.sym('eta_i', 2*(dim_state-1));
nominal_eta_var = MX.sym('nominal_eta_var', 2*(dim_state-1));
e_i = MX.sym('e_i', dim_state-1);
de_i = MX.sym('de_i', dim_state-1);
h_i_dis_var = MX.sym('h_i_dis_var', dim_state-1);
dh_i_dis_var = MX.sym('dh_i_dis_var', dim_state-1);
ddh_i_dis_var = MX.sym('ddh_i_dis_var', dim_state-1);
nu_var = MX.sym('nu_var', dim_tau);
P_j_var = MX.sym('P_j_var', dim_state-1, mrs_size);
dP_j_var = MX.sym('dP_j_var', dim_state-1, mrs_size);
ddP_j_var = MX.sym('ddP_j_var', dim_state-1, mrs_size);
Delta_var = MX.sym('Delta_var', dim_state-1, mrs_size);
dDelta_var = MX.sym('dDelta_var', dim_state-1, mrs_size);
ddDelta_var = MX.sym('dDelta_var', dim_state-1, mrs_size);

%%% trajectory error matrices
calculate_torque = bar_B_q\(bar_M_q*nu_var + bar_C_q*Phi);
calculate_torque_func = Function('calculate_torque_func', {x_state, nu_var}, {calculate_torque});

To_cartesian = eye(2, 3);
copy_over_network = ones(1, mrs_size);
Leader_position_matrix = To_cartesian*leader_state*copy_over_network;
Leader_velocity_matrix = To_cartesian*dleader_state*copy_over_network;
Leader_acceleration_matrix = To_cartesian*ddleader_state*copy_over_network;

tracking_error = ((P_j_var - Delta_var)*modified_laplacian...
                      - Leader_position_matrix*leader_adjacency)*basis_var;
tracking_error_rate = -(dP_j_var*adjacency + dDelta_var*modified_laplacian...
                           + Leader_velocity_matrix*leader_adjacency)*basis_var...
                           + g_i_dis_var*S_y_q*Phi;
eta = [tracking_error; tracking_error_rate];
eta_func = Function('eta_func', {P_j_var, dP_j_var, Delta_var, dDelta_var,...
                                 leader_state, dleader_state,...
                                 x_state, g_i_dis_var, basis_var},...
                                {eta});

tilde_n_i = g_i_dis_var*S_y_q;
inv_tilde_n_i = inv_S_y_q/g_i_dis_var;
tilde_m_i = g_i_dis_var*dS_y_q*Phi;
tilde_l_i = -(ddP_j_var*adjacency + ddDelta_var*modified_laplacian...
            + Leader_acceleration_matrix*leader_adjacency)*basis_var;

tracking_error_accel = tilde_l_i + tilde_m_i + tilde_n_i*nu_var;
deta = [tracking_error_rate; tracking_error_accel];
deta_func = Function('deta_func', {dP_j_var, ddP_j_var, dDelta_var, ddDelta_var,...
                                   dleader_state, ddleader_state,...
                                   x_state, g_i_dis_var, basis_var, nu_var},...
                                  {deta});

convert_xi_to_nu = inv_tilde_n_i*(xi_var - tilde_l_i - tilde_m_i);
convert_xi_to_nu_func = Function('convert_xi_to_nu_func',...
                                 {ddP_j_var, ddDelta_var, ddleader_state,...
                                  x_state, g_i_dis_var, basis_var, xi_var},...
                                 {convert_xi_to_nu});

%%% formation parameters
Delta = [0.5, -0.5*cosd(60), -0.5*cosd(60);
           0, -0.5*sind(60),  0.5*sind(60)];

%%% position of obstacles
robot_radius = 0.3;
obstacle_radius = 0.5;
d_safe = obstacle_radius + robot_radius;
p_obs_1 = [-1; 1.5];
p_obs_2 = [5.5; 0];
p_obs_3 = [5.5; 3];
obstacle_array = [p_obs_1 p_obs_2 p_obs_3];

%%% leader trajectory (virtual signal generation)
initial_position = [-0.2; -1.5];
goal_position = [2.2; 4.8];
heading_vector = goal_position - initial_position;
initial_angle = atan2(heading_vector(2), heading_vector(1));
start = [initial_position; initial_angle];

p_l = start;
dp_l = zeros(dim_state, 1);
ddp_l = zeros(dim_state, 1);
Delta_history = Delta;
dDelta_history = zeros(dim_state-1, mrs_size);
ddDelta_history = zeros(dim_state-1, mrs_size);

q_l = start; % iniital state of virtual robot

%%% define velocities of virtual leader over trajectory
v_l = zeros(2, T);

%%% initial conditions
q_k = zeros(dim_state, mrs_size); % empty vector for states

% Initial Robot possitions and orientations
p1_0 = [-0.75; -1.3975; 0.327];
p2_0 = [-0.498; -2.802; 0.799];
p3_0 = [-1.685; -1.397; 0.718];

% zero starting error
% p1_0 = [-0.022; -1.0328; initial_angle];
% p2_0 = [0.1156; -1.8878; initial_angle];
% p3_0 = [-0.6936; -1.5795; initial_angle];

q_next_k = [p1_0 p2_0 p3_0]; % initial system state
dq_next_k = zeros(dim_state, mrs_size); % initial system velocity
robot_accelerations = zeros(dim_state-1, mrs_size);

%%% cost weighting matrices (tune)
S_matrix = [zeros(2), eye(2)];
Q = diag([50, 50, 35, 35]);
R = 10*eye(dim_state-1);
R_horizon = sparse(kron(eye(N), R));
H = [zeros(dim_state-1) eye(dim_state-1);
     zeros((dim_state-1), 2*(dim_state-1))];
G = [zeros(dim_state-1); eye(dim_state-1)];

% discrete consensus error model
auxiliary_sys = ss(H, G, [], []);
auxiliary_sysd = c2d(auxiliary_sys, dt);
dt_H = auxiliary_sysd.A;
dt_G = auxiliary_sysd.B;

% Construct Phi (n x Nm)
Phi_matrix = zeros(2*(dim_state-1), N*(dim_state-1));
dt_Hk = eye(2*(dim_state-1));
for i = N:-1:1
    Phi_matrix(:, (N - i) * (dim_state-1) + (1:dim_state-1)) = dt_Hk * dt_G;
    dt_Hk = dt_Hk * dt_H;
end
dt_H_N = dt_H^N;

[P,~,K] = dare(dt_H, dt_G, Q, R);

% Define Optimisation Problem Bounds
z_var = 1/(g_i_dis_var*l_T);
c_var = scaled_radius*z_var;
lambda = 0.5*J_B*c_var*z_var;

H_theta = lambda*[-sin(2*theta) cos(2*theta); 
                   cos(2*theta) sin(2*theta)];
f_theta = c_var*sigma_val*[cos(theta + alpha_val) sin(theta + alpha_val);
                           cos(theta - alpha_val) sin(theta - alpha_val)];
d_E_theta = ((de_i - dh_i_dis_var)')*H_theta*(de_i - dh_i_dis_var)*[-1; 1];
torque_bound = d_E_theta + f_theta*(xi_var - ddh_i_dis_var);
torque_bound_func = Function('torque_bound_func', {de_i, xi_var, theta, dh_i_dis_var,...
                                                    ddh_i_dis_var, g_i_dis_var},...
                                                  {torque_bound});
wheel_vel_bound = inv_tilde_n_i*(de_i - dh_i_dis_var);
wheel_vel_bound_func = Function('wheel_vel_bound_func', {de_i, theta, dh_i_dis_var, g_i_dis_var},...
                                                  {wheel_vel_bound});

Phi_e_i = inv_tilde_n_i*(de_i - dh_i_dis_var);
next_theta = theta + S_q(3, :)*Phi_e_i*dt;
next_theta_func = Function('next_theta_func', {de_i, theta, dh_i_dis_var, g_i_dis_var},...
                                              {next_theta});

%%% memory allocation
consensus_opt_output = repmat(repmat(zeros(dim_state-1, 1), N, 1), 1, mrs_size);
consensus_opt_params = MX.sym('consensus_opt_params', length(consensus_opt_output));

v_i = repmat(zeros(dim_vel, 1), 1, mrs_size);
wheel_vel_i = (1/(r*l_T))*[l_T l_R; l_T -l_L]*v_i;
tau_i = repmat(zeros(dim_tau, 1), 1, mrs_size);
q_i_library = repmat(zeros(dim_state, length(k)), 1, 1, mrs_size);
formation_error = repmat(zeros(dim_state-1, length(k)), 1, 1, mrs_size);
vel_library = repmat(zeros(dim_vel, length(k)), 1, 1, mrs_size);
tau_library = repmat(zeros(dim_tau, length(k)), 1, 1, mrs_size);
wheeled_vel_library = repmat(zeros(dim_vel, length(k)), 1, 1, mrs_size);

%%% Define IPOPT options
opts = struct;
opts.ipopt.tol = 1e-8; % Tolerance for optimization
opts.ipopt.constr_viol_tol = 1e-9; % Constraint violation robot_radius
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_constr_viol_tol = 1e-8;
opts.ipopt.max_iter = 120;  % Max number of iterations
opts.ipopt.print_level = 0; % Suppress output (0 for no output)
opts.print_time = false;
opts.record_time = true;
opts.ipopt.sb = 'yes'; % Suppress the IPOPT startup banner

p_diffeomorphism_k_init = MX.sym('p_diffeomorphism_k_init', 4);
p_h_i_dis = MX.sym('p_h_i_dis', 2);
p_dh_i_dis = MX.sym('p_dh_i_dis', 2);
p_ddh_i_dis = MX.sym('p_ddh_i_dis', 2);
p_angle_init = MX.sym('p_angle_init');
p_g_i_dis = MX.sym('p_g_i_dis');

% Initialise symbolic states and cost using the parameters
p_diffeomorphism_k_next = p_diffeomorphism_k_init;
h_i_dis_next = p_h_i_dis;
h_vel = p_dh_i_dis;
angle = p_angle_init;
consensus_stage_cost = 0;

torque_constraint_vec = MX((N+1)*(dim_state-1), 1);
wheel_vel_constraint_vec = MX((N+1)*(dim_state-1), 1);

for z=1:N
    if (z == 1)
        h_accel  = p_ddh_i_dis;  
    else
        h_accel  = MX.zeros(2, 1);
    end

    opt_var_z = consensus_opt_params((z-1)*(dim_state-1)+1:z*(dim_state-1));
    angle_predict = next_theta_func(p_diffeomorphism_k_next(3:4), angle, h_vel, p_g_i_dis);
    p_diffeomorphism_predict = dt_H*p_diffeomorphism_k_next + dt_G*opt_var_z;

    torque_var = torque_bound_func(p_diffeomorphism_k_next(3:4), opt_var_z,...
                                   angle, h_vel, h_accel, p_g_i_dis);
    wheel_bound_var = wheel_vel_bound_func(p_diffeomorphism_predict(3:4),...
                                           angle_predict, h_vel, p_g_i_dis);

    torque_constraint_vec((z-1)*(dim_state-1)+1:z*(dim_state-1)) = torque_var;
    wheel_vel_constraint_vec((z-1)*(dim_state-1)+1:z*(dim_state-1)) = wheel_bound_var;
    
    p_diffeomorphism_k_next = p_diffeomorphism_predict;
    h_i_dis_next = h_i_dis_next + h_vel*dt;
    h_vel = h_vel + h_accel*dt;
    angle = angle_predict;
    
    if (z == N)
        consensus_stage_cost = consensus_stage_cost + (p_diffeomorphism_k_next')*P...
                                *p_diffeomorphism_k_next;
        terminal_control_law = -K*p_diffeomorphism_k_next;
        terminal_state = dt_H*p_diffeomorphism_k_next + dt_G*terminal_control_law;
        angle_predict = next_theta_func(p_diffeomorphism_k_next(3:4), angle, h_vel, p_g_i_dis);

        torque_var = torque_bound_func(p_diffeomorphism_k_next(3:4), terminal_control_law,...
                                         angle, h_vel, h_accel, p_g_i_dis);
        wheel_bound_var = wheel_vel_bound_func(terminal_state(3:4),...
                                               angle_predict, h_vel, p_g_i_dis);
        torque_constraint_vec(z*(dim_state-1)+1:(z+1)*(dim_state-1)) = torque_var;
        wheel_vel_constraint_vec(z*(dim_state-1)+1:(z+1)*(dim_state-1)) = wheel_bound_var;
    else
        consensus_stage_cost = consensus_stage_cost + (p_diffeomorphism_k_next')*Q...
                                *p_diffeomorphism_k_next;
    end
end

%%% Define NLP problem
formation_sys_objective = consensus_stage_cost...
                          + (consensus_opt_params(1:N*(dim_state-1))')*R_horizon ...
                             *consensus_opt_params(1:N*(dim_state-1));

formation_sys_constraint = [torque_constraint_vec;
                            wheel_vel_constraint_vec];
%%% optimisation parameter bounds
formation_lbg = [repmat(torque_lb, N+1, 1); repmat(wheel_vel_lb, N+1, 1)];
formation_ubg = [repmat(torque_ub, N+1, 1); repmat(wheel_vel_ub, N+1, 1)];

% Assemble all symbolic parameters into a single vector
all_params_sym = [p_diffeomorphism_k_init; p_h_i_dis; p_dh_i_dis; p_ddh_i_dis; 
                  p_angle_init; p_g_i_dis];

% Create the NLP solver object ONCE
consensus_nlp = struct('x', consensus_opt_params, 'p', all_params_sym, ...
                       'f', formation_sys_objective, 'g', formation_sys_constraint);
consensus_nlpsolve = nlpsol('solver', 'ipopt', consensus_nlp, opts);
solver_times = zeros(mrs_size, length(k));

h = waitbar(0, 'Please wait... Running Simulation');

%% Simulate Distributed Multirobot System Control
for j=1:length(k)
    waitbar(j / length(k), h, sprintf('Simulation Progress: %d%%', floor((j/length(k))*100)));
    
    p_j_array_database = q_next_k(1:dim_state-1, :);
    p_jdot_database = dq_next_k(1:dim_state-1, :);
    P_jddot_database = robot_accelerations;

    formation_radius = max(vecnorm(p_j_array_database - To_cartesian*q_l, 2, 1));
    rho_0 = formation_radius + d_safe;
    goal_obs_dist = min(vecnorm(obstacle_array - goal_position, 2, 1));

    k_rep = 1;
    term_a = 2*goal_obs_dist/(27*rho_0^3);
    term_b = 2/(3*rho_0^2);
    k_att = ceil(k_rep*(term_a - term_b + (term_b/3 + term_a)*sqrt(1 + 3*(rho_0/goal_obs_dist))));
    F_att = k_att*(goal_position - To_cartesian*q_l);

    [min_robot_obs_dist, obs_index] = min(vecnorm(obstacle_array - To_cartesian*q_l, 2, 1));
    if (min_robot_obs_dist <= rho_0)
        rho_vc_goal = norm(goal_position - To_cartesian*q_l);
        F_rep_1 = k_rep*(1/min_robot_obs_dist - 1/rho_0)*((rho_vc_goal^2)/(min_robot_obs_dist^3));
        F_rep_2 = -k_rep*((1/min_robot_obs_dist - 1/rho_0)^2);
        F_rep = F_rep_1*(To_cartesian*q_l - obstacle_array(:, obs_index))...
                + F_rep_2*(To_cartesian*q_l - goal_position);
    else
        F_rep = zeros(dim_state-1, 1);
    end

    F_total = 0.02*(F_att + F_rep);
    velocity_vector = full((R_theta_func(q_l(3))')*F_total);
    v_l(1, j) = velocity_vector(1);
    v_l(2, j) = velocity_vector(2)/l_T;

    q_l_solution = q_solver('x0', q_l, 'p', v_l(:, j));
    q_l = full(q_l_solution.xf);
    dq_l = R_q_func(q_l)*v_l(:, j);
    p_l(:, j) = q_l;
    dp_l(:, j) = full(dq_l);
    ddp_l(:, j) = zeros(dim_state, 1);
    Delta_history(:, :, j) = full((R_theta_func(q_l(3)))*Delta);
    dDelta_history(:, :, j) = full(S_omega_func(dq_l(3))*(R_theta_func(q_l(3)))*Delta);
    ddDelta_history(:, :, j) = full(S_omega_func(0)*(R_theta_func(q_l(3)))...
                                    + (S_omega_func(dq_l(3))^2)*(R_theta_func(q_l(3))))*Delta;

    leader_distribution = To_cartesian*p_l(:, j)*copy_over_network;
    leader_velocity_distribution = To_cartesian*dp_l(:, j)*copy_over_network;
    leader_acceleration_distribution = To_cartesian*ddp_l(:, j)*copy_over_network;
    Delta_j = Delta_history(:, :, j);
    dDelta_j = dDelta_history(:, :, j);
    ddDelta_j = ddDelta_history(:, :, j);

    global_formation_errors = (p_j_array_database - Delta_j)*modified_laplacian...
                                - leader_distribution*leader_adjacency;
    velocity_sums = p_jdot_database*adjacency + leader_velocity_distribution*leader_adjacency;

    for i=1:mrs_size
        g_i_dis_val = g_dis(i,i);
        q_k(:, i) = q_next_k(:, i);
        measured_q = q_k(:, i);
        x_state_k = [measured_q; wheel_vel_i(:, i)];
        angle_val = x_state_k(3);
        p_i = p_j_array_database(:, i);

        q_i_library(:, j, i) = measured_q;
        vel_library(:, j, i) = v_i(:, i);
        wheeled_vel_library(:, j, i) = wheel_vel_i(:, i);

        leader_error = leader_neighbour(i)*((p_i - p_l(1:dim_state-1, j)) - Delta_j(:, i));
        global_formation_error = global_formation_errors(:, i);
        % Evaluate local formation error excluding trajectory error
        formation_error(:, j, i) = abs(global_formation_error);  % - leader_error
        ith_basis_vector = std_basis(mrs_size, i);
        h_i_dis_val = -(p_j_array_database*adjacency + Delta_j*modified_laplacian...
                        + leader_distribution*leader_adjacency)*ith_basis_vector;
        dh_i_dis_val = -(p_jdot_database*adjacency + dDelta_j*modified_laplacian...
                        + leader_velocity_distribution*leader_adjacency)*ith_basis_vector;
        ddh_i_dis_val = -(P_jddot_database*adjacency + ddDelta_j*modified_laplacian...
                        + leader_acceleration_distribution*leader_adjacency)*ith_basis_vector;
        diffeomorphism_k_val = eta_func(p_j_array_database, p_jdot_database,...
                                     Delta_j, dDelta_j, p_l(:, j), dp_l(:, j),...
                                     x_state_k, g_i_dis_val, ith_basis_vector);
        relative_average_velocity_val = velocity_sums*ith_basis_vector/g_i_dis_val;

        %%% Consensus Errors
        % Pack numerical data into a vector in the SAME ORDER as 'all_params_sym'
        numerical_param_values = [diffeomorphism_k_val; h_i_dis_val; dh_i_dis_val;...
                                  ddh_i_dis_val; angle_val; g_i_dis_val];
        
        % Call the pre-built solver with the numerical data
        consensus_horizon_variables = consensus_nlpsolve('x0', consensus_opt_output(:, i), ...
                                                         'p', numerical_param_values, ...
                                                         'lbg', formation_lbg,...
                                                         'ubg', formation_ubg);

        consensus_opt_output(:, i) = full(consensus_horizon_variables.x);
        optimal_Xi = full(consensus_opt_output(1:N*(dim_state-1), i));
        optimal_xi = full(consensus_opt_output(1:dim_state-1, i));
        optimal_nu = convert_xi_to_nu_func(P_jddot_database, ddDelta_j, ddp_l(:, j),...
                                           x_state_k, g_i_dis_val, ith_basis_vector,...
                                           optimal_xi);
        tau_i(:, i) = full(calculate_torque_func(x_state_k, optimal_nu));

        %%% update states
        state_solution_f = state_solver('x0', x_state_k, 'p', tau_i(:, i));
        x_state_next = full(state_solution_f.xf);
        dx_state_next = full(dx_state_func(x_state_k, tau_i(:, i)));
        robot_i_acceleration = full(global_robot_acceleration_func(x_state_k, tau_i(:, i)));

        q_next_k(:, i) = x_state_next(1:dim_state);
        dq_next_k(:, i) = dx_state_next(1:dim_state);
        robot_accelerations(:, i) = robot_i_acceleration;
        wheel_vel_i(:, i) = x_state_next(end-dim_vel+1:end);
        v_i(:, i) = [norm(dq_next_k(1:2, i)); dq_next_k(3, i)];
        eta_i_N  = full(dt_H_N*diffeomorphism_k_val + Phi_matrix*optimal_Xi);
        next_optimal_Xi = [optimal_Xi(dim_state:end); -K*eta_i_N];
        consensus_opt_output(1:N*(dim_state-1), i) = next_optimal_Xi;

        tau_library(:, j, i) = tau_i(:, i);
        solver_times(i, j) = consensus_nlpsolve.stats().t_wall_total;
    end
end

close(h);

save('APF_Trajectory_Values.mat', 'formation_error', 'q_i_library', ...
      'p_obs_1', 'p_obs_2', 'p_obs_3', 'obstacle_radius', 'k', 'p_l',...
      'vel_library', 'v_l', 'goal_position', 'initial_position');

%% Calculation of Performance Metrics
square_errors = formation_error.^2;
square_norm_errors = sum(square_errors, 1);
integral_of_norm_errors = sum(square_norm_errors, 2);
avg_consensus_error = mean(squeeze(sqrt(integral_of_norm_errors)));

avg_solver_time = mean(solver_times*1000, 'all');

square_torques = tau_library.^2;
square_norm_torque = sum(square_torques, 1);
integral_of_norm_torques = sum(square_norm_torque, 2);
avg_energy_consumed = mean(squeeze(sqrt(integral_of_norm_torques)));

TV = sum(abs(diff(tau_library))); % Total Variation
avg_TV = mean(squeeze(TV));

disp('Controller Metrics:')
disp([avg_consensus_error, avg_solver_time, avg_energy_consumed, avg_TV]);

%% Plot Results
%%% trajectories
figure(1);
hold on;
plot(q_i_library(1, :, 1), q_i_library(2, :, 1), 'r', 'LineWidth', 1.5);
plot(q_i_library(1, :, 2), q_i_library(2, :, 2), 'b', 'LineWidth', 1.5);
plot(q_i_library(1, :, 3), q_i_library(2, :, 3), 'g', 'LineWidth', 1.5);
plot(p_l(1, 1:length(k)), p_l(2, 1:length(k)), 'c--', 'LineWidth', 1.5);
rectangle('Position', [p_obs_1(1)-obstacle_radius,...
                       p_obs_1(2)-obstacle_radius,...
                       2*obstacle_radius, 2*obstacle_radius],...
                      'Curvature', [1, 1], 'FaceColor',[0 0 0]);
rectangle('Position', [p_obs_2(1)-obstacle_radius,...
                       p_obs_2(2)-obstacle_radius,...
                       2*obstacle_radius, 2*obstacle_radius],...
                      'Curvature', [1, 1], 'FaceColor',[0 0 0]);
rectangle('Position', [p_obs_3(1)-obstacle_radius,...
                       p_obs_3(2)-obstacle_radius,...
                       2*obstacle_radius, 2*obstacle_radius],...
                      'Curvature', [1, 1], 'FaceColor',[0 0 0]);

cord1x = reshape(q_i_library(1, 1, 1:3), [], 1);
cord1y = reshape(q_i_library(2, 1, 1:3), [], 1);
plot([cord1x; cord1x(1); p_l(1, 1)], [cord1y; cord1y(1); p_l(2, 1)], '.-',...
    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);

plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);
plot(nan, nan, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10);

plot(goal_position(1), goal_position(2),'rp', 'MarkerSize', 8);
plot(initial_position(1), initial_position(2),'bp', 'MarkerSize', 8);

time = int32(1*length(k)/5);
cord1x = reshape(q_i_library(1, time, 1:3), [], 1);
cord1y = reshape(q_i_library(2, time, 1:3), [], 1);
plot([cord1x; cord1x(1); p_l(1, time)], [cord1y; cord1y(1); p_l(2, time)], '.-',...
    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);

time = int32(3.5*length(k)/5);
cord1x = reshape(q_i_library(1, time, 1:3), [], 1);
cord1y = reshape(q_i_library(2, time, 1:3), [], 1);
plot([cord1x; cord1x(1); p_l(1, time)], [cord1y; cord1y(1); p_l(2, time)], '.-',...
    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);

cord1x = reshape(q_i_library(1, length(k), 1:3), [], 1);
cord1y = reshape(q_i_library(2, length(k), 1:3), [], 1);
plot([cord1x; cord1x(1); p_l(1, length(k))], [cord1y; cord1y(1); p_l(2, length(k))], '.-',...
    'Color', [0.5, 0.5, 0.5, 0.5],'MarkerSize', 10, 'LineWidth', 2);
hold off;
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3', 'Leader', 'formation tracker', ...
       'Obstacle1', 'Obstacle2', 'Obstacle3', 'Start', 'Goal');
title("Robots' Trajectories", 'fontsize', 15);
xlabel('X coordinate [m]', 'fontSize', 12);
xlim([-3 8]);
ylabel('Y coordinate [m]', 'fontSize', 12);
ylim([-3 6]);
pbaspect([diff(xlim) diff(ylim) 1]);
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% linear velocities
figure(2);
hold on;
plot(k, vel_library(1, :, 1), 'r', 'LineWidth', 1.5);
plot(k, vel_library(1, :, 2), 'b', 'LineWidth', 1.5);
plot(k, vel_library(1, :, 3), 'g', 'LineWidth', 1.5);
plot(k, v_l(1, 1:length(k)), 'c--', 'LineWidth', 1.5);
hold off;
title('Linear Velocities', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('Linear velocities v_i [m/s]', 'fontSize', 12);
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3', 'Leader');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% angular velocities
figure(3);
hold on;
plot(k, vel_library(2, :, 1), 'r', 'LineWidth', 1.5);
plot(k, vel_library(2, :, 2), 'b', 'LineWidth', 1.5);
plot(k, vel_library(2, :, 3), 'g', 'LineWidth', 1.5);
plot(k, v_l(2, 1:length(k)), 'c--', 'LineWidth', 1.5);
hold off;
title('Angular Velocities', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('Angular velocities $\omega_i$ [rad/s]', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3', 'Leader');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% X formation errors
figure(4);
hold on;
plot(k, formation_error(1, :, 1), 'r', 'LineWidth', 1.5);
plot(k, formation_error(1, :, 2), 'b', 'LineWidth', 1.5);
plot(k, formation_error(1, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Absolute longitudinal axis formation error (x)', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('x-position error', 'fontSize', 12);
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% Y formation errors
figure(5);
hold on;
plot(k, formation_error(2, :, 1), 'r', 'LineWidth', 1.5);
plot(k, formation_error(2, :, 2), 'b', 'LineWidth', 1.5);
plot(k, formation_error(2, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Absolute lateral axis formation error (y)', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('y-position error', 'fontSize', 12);
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% right wheel angular velocities
figure(6);
hold on;
plot(k, wheeled_vel_library(1, :, 1), 'r', 'LineWidth', 1.5);
plot(k, wheeled_vel_library(1, :, 2), 'b', 'LineWidth', 1.5);
plot(k, wheeled_vel_library(1, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Right Wheel Angular Velocities', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('$\dot{\Phi}_r$ [rad/s]', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% left wheel angular velocities
figure(7);
hold on;
plot(k, wheeled_vel_library(2, :, 1), 'r', 'LineWidth', 1.5);
plot(k, wheeled_vel_library(2, :, 2), 'b', 'LineWidth', 1.5);
plot(k, wheeled_vel_library(2, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Left Wheel Angular Velocities', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('$\dot{\Phi}_l$ [rad/s]', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');

%%% right wheel torque
figure(8);
hold on;
plot(k, tau_library(1, :, 1), 'r', 'LineWidth', 1.5);
plot(k, tau_library(1, :, 2), 'b', 'LineWidth', 1.5);
plot(k, tau_library(1, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Right Wheel Torque', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('$\tau_r\ [N\cdot m]$', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);

%%% left wheel torque
figure(9);
hold on;
plot(k, tau_library(2, :, 1), 'r', 'LineWidth', 1.5);
plot(k, tau_library(2, :, 2), 'b', 'LineWidth', 1.5);
plot(k, tau_library(2, :, 3), 'g', 'LineWidth', 1.5);
hold off;
title('Left Wheel Torque', 'fontsize', 15);
xlabel('time [s]', 'fontSize', 12);
ylabel('$\tau_l\ [N\cdot m]$', 'FontSize', 12, 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
legend('Robot 1', 'Robot 2', 'Robot 3');
box on
grid on
grid minor
set(gca, 'fontsize', 12);
set(gca, 'TickDir', 'out');