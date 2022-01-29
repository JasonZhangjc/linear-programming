function [soln obj status]=SimplexMethod(A, b, c, B_bar)
% The function SimplexMethod uses the simplex method to 
% solve an LP in standard form:
%       min c' * x
%  s.t. A * x == b
%           x >= 0
%
% Inputs:
% c = n*1 vector of objective coefficients
% A = m*n matrix with m < n
% b = m*1 vector of RHS coefficients
% B_bar = m*1 vector that contains indices of basic variables
%
% Output:
% soln is an n*1 vector, the final solution of LP
% obj is a number, the final objective value of LP
% status describes the termination condition of the problem as follows:
% 1 - optimal solution
% -3 - unbounded problem
% If the LP is unbounded, soln returns the correspondig extreme direction

soln=[]; obj=[]; status=[];

%% Step 0: Initialization
% Generate an initial basic feasible solution and partition x, A, and c
% so that x=[x_B | x_N], A=[B | N], and c=[c_B | c_N]
indices = [1:length(c)]';
indices(find(ismember(indices, B_bar)==1)) = [];
N_bar = indices;
B = A(:,B_bar);           %basis matrix B
N = A(:,N_bar);           %non-basis matrix N
x_B = B\b;                %basic variables 
x_N = zeros(length(N_bar),1); %non-basic variables
x = [x_B; x_N];           %partition x according to basis
c_B = c(B_bar);           %obj coefficients of basic variables
c_N = c(N_bar);           %obj coefficients of non-basic variables
obj_init = [c_B; c_N]'*x; %initial objective function value
k = 0;

while k >= 0
    %% Step 1: Checking Optimality
    % Compute the reduced costs r_q=c_q-c_B¡¯*B^(-1)*N_q for q in N_bar
    % if r_q >= 0, STOP, current solution is optimal, else go to STEP 2
    pi = B'\c_B;        %solve the system B^T*pi=c_B for simplex multipliers
    r_N = c_N' - pi'*N; %compute reduced cost for non-basic variables
    ratio = find(r_N<0);
    if isempty(ratio)   %if r_q >= 0, then STOP. Optimal
        disp('probelm solved')
        status = 1;
        obj = [c_B; c_N]'*x;
        %indices of x are in increasing order
        indices_temp = [B_bar; N_bar];
        for a = 1:length(c)
            soln(a,1) = x(find(indices_temp==a));
        end
        break
    else % if r_q < 0, GO TO Step 2
        %% Step 2: Generating Direction Vector
        % Construct d_q=[-B^(-1)*N; e_q].
        % If d_q >= 0, then LP is unbounded, STOP, else go to STEP 3.
        enter_idx = ratio(1);  %choosing entering variable
        e = zeros(length(N_bar),1);
        e(enter_idx) = 1;      %construct vector e_q
        d = -B\(N(:,enter_idx)); %solve Bd=-N_q
        direction = [d; e];    %improved direction d
        d_idx = find(direction < 0);
        if isempty(d_idx)     %if direction > 0, then STOP.(unbounded)
            disp('unbounded problem')
            status=-3;
            indices_temp = [B_bar; N_bar];
            for a = 1:length(c)
                soln(a,1) = direction(find(indices_temp==a));
            end
            break
        else %if d_q < 0, GO TO Step 3
            %% Step 3: Generating Step Length
            % Compute step length by the minimum ratio test. Go to STEP 4.
            step_indices = -x(d_idx)./direction(d_idx);
            step = min(step_indices);
            for i = 1:length(d_idx)
                if step == -x(d_idx(i))/direction(d_idx(i))
                    leave_idx = d_idx(i);
                end
            end
            %% Step 4: Generating Improved Solution
            % Let x(k+1) = x(k) + alpha*d_q. Go to Step 5.
            x_d = x + step*direction;
            %leave_indices = find(x_d(1:length(B_bar))==0);
            %leave_idx = leave_indices(1); %determining leaving variable
            %% Step 5: Updating Basis
            % Generate the new basis B for next iteration,
            % Update c=[c_B|c_N],x=[x_B|x_N],& Aeq=[B|N]. Go to STEP 1.
            B_bar_temp = B_bar;
            N_bar_temp = N_bar;
            x_B = x_d(1:length(B_bar));
            x_N = x_d(length(B_bar)+1:end);
            x_B_temp = x_d(1:length(B_bar));
            x_N_temp = x_d(length(B_bar)+1:end);
            %exchange the entering and leaving variables in B_bar
            B_bar(leave_idx) = N_bar_temp(enter_idx);
            N_bar(enter_idx) = B_bar_temp(leave_idx);
            x_B(leave_idx) = x_N_temp(enter_idx);
            x_N(enter_idx) = x_B_temp(leave_idx);
            x = [x_B; x_N];           %update x = [x_B | x_N]
            B = A(:,B_bar);           %update basis B
            N = A(:,N_bar);           %update non-basis N
            c_B = c(B_bar);           %update c_B
            c_N = c(N_bar);           %update c_N
            obj_init = [c_B; c_N]'*x; %update objective value
            k = k+1; %GO TO Step 1
        end
    end
end