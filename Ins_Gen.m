%% Problem Instances Generation
clear all;
num_s = 100;   % number of scenarios
num_i = 10;  % for each scenario, run a number of instances
num_sub = 3; % number of subproblems, set to 3

%% Randomize Nonzero Parameters
para_1 = [];
para_2 = [];
para_3 = [];
for i = 1:num_s
    para_1(i) = 2+randi(100)/100;
    para_2(i) = 2+randi(200)/100;
    para_3(i) = -10-randi(200)/10;
end

% para_1 = [3,2.5,2];
% para_2 = [3.6,3,2.4];
% para_3 = [-24,-20,-16];

%% LP Instance with Unbounded Subproblems
% linking constraints
L1 = [1,zeros(1,num_s*3)];
L2 = [1,zeros(1,num_s*3)];
L3 = [1,zeros(1,num_s*4)];

% master problem parameters
master_u.b = 500;
master_u.L{1} = L1;
master_u.L{2} = L2;
master_u.L{3} = L3;

% objective coefficient vector: c
c1 = [150];
c2 = [230];
c3 = [260];
for i = 1:num_s
    c1 = [c1;238/num_s;-170/num_s];
    c2 = [c2;210/num_s;-150/num_s];
    c3 = [c3;-36/num_s;-10/num_s];
end
c1 = [c1;zeros(num_s,1)];
c2 = [c2;zeros(num_s,1)];
c3 = [c3;zeros(2*num_s,1)];

% RHS vector: b
b1 = [200*ones(num_s,1)];
b2 = [240*ones(num_s,1)];
b3 = [zeros(num_s,1);6000*ones(num_s,1)];

% constraint matrix of subproblems: A
A1 = zeros(num_s,num_s*3+1);
A2 = zeros(num_s,num_s*3+1);
A3 = zeros(2*num_s,num_s*4+1);
A1(:,1) = para_1';
A2(:,1) = para_2';
A3(:,1) = [para_3';zeros(num_s,1)];
for i = 1:num_s
    A1(i,i*2:i*2+1) = [1,-1];
    A2(i,i*2:i*2+1) = [1,-1];
    A3(i,i*2:i*2+1) = [1,1];
end
A1(:,num_s*2+2:end) = -1*eye(num_s);
A2(:,num_s*2+2:end) = -1*eye(num_s);
for i = 1:num_s
    A3(num_s+i,2*i) = 1;
end
A3(:,num_s*2+2:end) = eye(num_s*2);

% basic variables
basic_v1 = [];
basic_v2 = [];
basic_v3 = [];
for i = 1:num_s
    basic_v1(i) = 2*i;
    basic_v2(i) = 2*i;
    basic_v3(i) = 2*i+1;
    basic_v3(num_s+i) = num_s*3+1+i;
end

% initial extreme points
v1_init = zeros(length(c1),1);
B1 = A1(:,basic_v1);
v1_init_t = B1\b1;
v1_init(basic_v1) = v1_init_t;

v2_init = zeros(length(c2),1);
B2 = A2(:,basic_v2);
v2_init_t = B2\b2;
v2_init(basic_v2) = v2_init_t;

v3_init = zeros(length(c3),1);
B3 = A3(:,basic_v3);
v3_init_t = B3\b3;
v3_init(basic_v3) = v3_init_t;

% master problem parameters
sub_u.c{1} = c1;
sub_u.c{2} = c2;
sub_u.c{3} = c3;
sub_u.b{1} = b1;
sub_u.b{2} = b2;
sub_u.b{3} = b3;
sub_u.A{1} = A1;
sub_u.A{2} = A2;
sub_u.A{3} = A3;
sub_u.basis{1} = basic_v1;
sub_u.basis{2} = basic_v2;
sub_u.basis{3} = basic_v3;
sub_u.v0{1} = v1_init;
sub_u.v0{2} = v2_init;
sub_u.v0{3} = v3_init;

%% Use DW-Decomposition with SimplexMethod as its LP Solver
t1 = clock;
[soln_u, fval_u, flag_u] = dw_dec(master_u, sub_u, num_sub);
t2 = clock;
t(2) = etime(t2,t1)

%% Directly Use SimplexMethod to Solve the LP
% need one more slack variable
c_ = [c1;c2;c3;0];
b_ = [500;b1;b2;b3];
A_sub = blkdiag(A1,A2,A3);
A_link = [L1,L2,L3];
A_ = [A_link;A_sub];
slack_link = zeros(length(b_),1);
slack_link(1) = 1;
A_ = [A_,slack_link];
B_bar = [];
for i = 1:num_s
    B_bar = [B_bar;2*i];
end
for i = 1:num_s
    B_bar = [B_bar;1+3*num_s+2*i];
end
for i = 1:num_s*2
    B_bar = [B_bar;2+6*num_s+2*num_s+1+i];
end
B_bar = [B_bar;length(c_)];

t1 = clock;
[soln_s, fval_s, flag_s] = SimplexMethod(A_, b_, c_, B_bar);
t2 = clock;
t(1) = etime(t2,t1)

%% Directly Use linprog (interior-point-legacy) to Solve the LP
options = optimoptions('linprog','Algorithm','interior-point-legacy');
t1 = clock;
[soln_lp1,fval_lp1,flag_lp1,output_lp1]=linprog(c_,[],[],A_,b_,zeros(length(c_),1),[],options);
t2 = clock;
t(4) = etime(t2,t1)

%% Directly Use linprog (interior-point) to Solve the LP
options = optimoptions('linprog','Algorithm','interior-point');
t1 = clock;
[soln_lp2,fval_lp2,flag_lp2,output_lp2]=linprog(c_,[],[],A_,b_,zeros(length(c_),1),[],options);
t2 = clock;
t(5) = etime(t2,t1)

%% Directly Use linprog (dual-simplex) to Solve the LP
options = optimoptions('linprog','Algorithm','dual-simplex');
t1 = clock;
[soln_lp3,fval_lp3,flag_lp3,output_lp3]=linprog(c_,[],[],A_,b_,zeros(length(c_),1),[],options);
t2 = clock;
t(6) = etime(t2,t1)

%% The Same LP Instance with Bounded Subproblems
% linking constraints
L1 = [1,zeros(1,num_s*3+1)];
L2 = [1,zeros(1,num_s*3+1)];
L3 = [1,zeros(1,num_s*4+1)];

% master problem parameters
master_b.b = 500;
master_b.L{1} = L1;
master_b.L{2} = L2;
master_b.L{3} = L3;

% objective coefficient vector: c
c1 = [150];
c2 = [230];
c3 = [260];
for i = 1:num_s
    c1 = [c1;238/num_s;-170/num_s];
    c2 = [c2;210/num_s;-150/num_s];
    c3 = [c3;-36/num_s;-10/num_s];
end
c1 = [c1;zeros(num_s+1,1)];
c2 = [c2;zeros(num_s+1,1)];
c3 = [c3;zeros(2*num_s+1,1)];

% RHS vector: b
b1 = [200*ones(num_s,1);500];
b2 = [240*ones(num_s,1);500];
b3 = [zeros(num_s,1);6000*ones(num_s,1);500];

% constraint matrix of subproblems: A
A1 = zeros(num_s+1,num_s*3+2);
A2 = zeros(num_s+1,num_s*3+2);
A3 = zeros(2*num_s+1,num_s*4+2);
A1(1:num_s,1) = para_1';
A2(1:num_s,1) = para_2';
A3(1:2*num_s,1) = [para_3';zeros(num_s,1)];
for i = 1:num_s
    A1(i,i*2:i*2+1) = [1,-1];
    A2(i,i*2:i*2+1) = [1,-1];
    A3(i,i*2:i*2+1) = [1,1];
end
A1(1:num_s,num_s*2+2:end-1) = -1*eye(num_s);
A2(1:num_s,num_s*2+2:end-1) = -1*eye(num_s);
A1(num_s+1,1) = 1;
A1(num_s+1,end) = 1;
A2(num_s+1,1) = 1;
A2(num_s+1,end) = 1;
for i = 1:num_s
    A3(num_s+i,2*i) = 1;
end
A3(1:2*num_s,num_s*2+2:end-1) = eye(num_s*2);
A3(2*num_s+1,1) = 1;
A3(2*num_s+1,end) = 1;

% basic variables
basic_v1 = [];
basic_v2 = [];
basic_v3 = [];
for i = 1:num_s
    basic_v1(i) = 2*i;
    basic_v2(i) = 2*i;
    basic_v3(i) = 2*i+1;
    basic_v3(num_s+i) = num_s*3+1+i;
end
basic_v1 = [basic_v1,num_s*3+2];
basic_v2 = [basic_v2,num_s*3+2];
basic_v3 = [basic_v3,num_s*4+2];

% initial extreme points
v1_init = zeros(length(c1),1);
B1 = A1(:,basic_v1);
v1_init_t = B1\b1;
v1_init(basic_v1) = v1_init_t;

v2_init = zeros(length(c2),1);
B2 = A2(:,basic_v2);
v2_init_t = B2\b2;
v2_init(basic_v2) = v2_init_t;

v3_init = zeros(length(c3),1);
B3 = A3(:,basic_v3);
v3_init_t = B3\b3;
v3_init(basic_v3) = v3_init_t;

% master problem parameters
sub_b.c{1} = c1;
sub_b.c{2} = c2;
sub_b.c{3} = c3;
sub_b.b{1} = b1;
sub_b.b{2} = b2;
sub_b.b{3} = b3;
sub_b.A{1} = A1;
sub_b.A{2} = A2;
sub_b.A{3} = A3;
sub_b.basis{1} = basic_v1;
sub_b.basis{2} = basic_v2;
sub_b.basis{3} = basic_v3;
sub_b.v0{1} = v1_init;
sub_b.v0{2} = v2_init;
sub_b.v0{3} = v3_init;

%% Use DW-Decomposition with SimplexMethod as its LP Solver
t1 = clock;
[soln_b, fval_b, flag_b] = dw_dec(master_b, sub_b, num_sub);
t2 = clock;
t(3) = etime(t2,t1)

% Timer Index:
% 1: Simplex
% 2: DW with Unbounded Subproblems
% 3: DW with Bounded Subproblems
% 4: linprog with interior-point-legacy
% 5: linprog with interior-point
% 6: linprog with dual-simplex
