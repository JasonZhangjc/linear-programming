function [soln, fval, flag] = dw_dec(master, sub, num_sub)
% This function implements the Dantzig-Wolfe Decomposition Algorithm
% This function is for farming problem only
% With slight modification, this function should apply to other problems
% This function works for both BOUNDED and UNBOUNDED subproblems
% This function uses SimplexMethod as its LP solver for subproblems
% Inputs:
% master: parameters of the master problem
% sub: parameters of the subproblems
% num_sub: number of subproblems
% num_sub = 3 in the farming problem
% Outputs:
% soln: solution vector if the LP is bounded
% fval: objective value if the LP is bounded
% flag: status of termination
% flag == 1: optimal solution found
% flag == -3: problem is unbounded
%% Step 0: Initialization
    soln = {};                   % solution
    fval = 0;                    % objective value
    flag = 0;                    % exit flag
    flag_optimal = 0;            % optimal flag
    flag_unbounded = 0;          % unbounded flag
    epsilon = 1e-16;             % termination threshold     
    B = eye(4);                  % basis matrix of restricted MP
    b = [500;ones(num_sub,1)];   % RHS vector of restricted MP
    sol.c{1} = [0];              % objective coefficients
    sol.x{1} = [500];            % extreme points/rays
    for i = 2:num_sub+1
        sol.c{i} = sub.c{i-1};
        sol.x{i} = sub.v0{i-1};
    end
    sol.p = [0;1;2;3];           % problem indices of extreme points
    
    f_B = zeros(num_sub+1,1);    % f_B
    x_simp = {};                 % solution from SimplexMethod
    fval_simp = {};              % objective value from SimplexMethod
    flag_simp = {};          % exitflag from SimplexMethod
    
    % DW-Decomposition Main loop
    while ((1-flag_optimal) && (1-flag_unbounded))
        unbnd_sub = 0; % The flag indicating the existence of an unbounded subproblem
        %% Step 1: Simplex Multiplier Generation
        for i = 1:num_sub+1
            f_B(i) = sol.c{i}'*sol.x{i};
        end
        x_B = B\b;               % x_B
        pie = B'\f_B;            % pie
        pie_1 = pie(1);          % Recall pie_1 and pie_2 in pie
        
        %% Step 2: Optimality Check
        % Use my SimplexMethod instead of linprog
        enter_var = -1;          % entering variable
        for i = 1:num_sub
            A_sub = sub.A{i};    % A matrix of a subproblem
            b_sub = sub.b{i};    % b vector of a subproblem
            c_sub = (sub.c{i}' - pie_1*master.L{i})'; % c vector of a subproblem
            basis_sub = sub.basis{i}';   % basis indices of a subproblem
            [x_simp{i},fval_simp{i},flag_simp{i}]=SimplexMethod(A_sub,b_sub,c_sub,basis_sub);
            sol.x_temp{i} = x_simp{i};   % temporary extreme points
            if flag_simp{i} == -3
                enter_var = i+1;
                r_N = [0;-1;-1;-1];      % reduced costs
                unbnd_sub = 1;
                break;
            end
        end

        if enter_var == -1
            r_N = [0;fval_simp{1}-pie(2);fval_simp{2}-pie(3);fval_simp{3}-pie(4)];
            enter_var = find(r_N==min(r_N));
        end
        
        if sum(r_N<-1e-3==1)==0    % if r_N >= 0, stop, optimal solution found
            disp('Optimal Solution Found !')
            flag_optimal = 1;
            flag = 1;
            % output objective value
            fval = f_B' * (B\b);
            display(fval)
            % output varaibles
            soln{1}=zeros(length(sub.c{1}),1);
            soln{2}=zeros(length(sub.c{2}),1);
            soln{3}=zeros(length(sub.c{3}),1);
            x_B_opt=B\b;
            for i = 1:num_sub+1
                prob = sol.p(i);
                soln{prob} = soln{prob} + x_B_opt(i)*sol.x{i};
            end  
            solution = [soln{1};soln{2};soln{3}];
            disp(solution)
            break;
        end
        
        %% Step 3: Column Generation
        a_bar=zeros(num_sub+1,1);         % a_bar
        enter_prob = sol.p(enter_var);    % subproblem index of entering var
        if enter_prob>0
            if unbnd_sub == 0
                a_bar(enter_prob+1)=1;
            else
                a_bar(enter_prob+1)=0;
            end
            a_bar(1)=master.L{enter_prob}*sol.x_temp{enter_prob};
        end
        if enter_prob==0
           a_bar(1)=1;
        end
        
        %% Step 4: Descent Direction Generation
        d = B\(-a_bar);                   % descent direction
        if sum(d<-epsilon==1)==0          % if # of all d>>0, stop, problem unbounded
            flag_unbounded=1;
            flag = -3;
            disp('Problem Unbounded !')
            break;
        end
        
        %% Step 5: Step Length Generation
        leave_cand = find(d<-epsilon);    % leaving variable candidates
        ratio=-x_B(leave_cand)./d(leave_cand);
        alpha=min(ratio);                 % alpha: min ratio test
        leave_cand = leave_cand(find(ratio==alpha));
        leave_var = leave_cand(1);        % leaving variable
        
        %% Step 7: Basis Update (Step 6 is done in another way)
        B(:,leave_var)=a_bar;
        if enter_prob>0
           f_B(leave_var) = sub.c{enter_prob}'*sol.x_temp{enter_prob};
           sol.x{leave_var} = sol.x_temp{enter_prob};
           sol.c{leave_var} = sub.c{enter_prob};
           sol.p(leave_var) = enter_prob;
        end

        if enter_prob==0
           f_B(leave_var)=B\b(leave_var);
           sol.c{leave_var} = [0];
           sol.p(leave_var) = 0;
        end
    end    
end