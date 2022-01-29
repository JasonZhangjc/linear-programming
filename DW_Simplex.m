clear all
%% Instance with Bounded Sub-Problems
compute_matrices = 1;
compute_initial_extreme_points=1;
construct_structure = 1;

% while (compute_matrices)
%     scene_num = 3;
%     epsilon = 1e-16;
%     flag_optimal=0;
%     flag_unbounded=0;
%     
%     % Initialization of Matrices and Vectors
%     % Subproblem 1 variables: x_1, y_1i, w_1i
%     % Subproblem 1 slack variables: z_1i
%     % Subproblem 1 x_1 upper bound slack variable z1
%     % numbers = 3*scene_num + 2
%     var_num_s1=3*scene_num+2; % with slack variables and x_1 <= 500
%     var_num_s2=3*scene_num+2;
%     var_num_s3=4*scene_num+2; % without x_4
%     
%     A1 = zeros(scene_num+1,var_num_s1);
%     A2 = zeros(scene_num+1,var_num_s2);
%     A3 = zeros(2*scene_num+1,var_num_s3);
%     
%     b1 = zeros(scene_num+1,1);
%     b2 = zeros(scene_num+1,1);
%     b3 = zeros(2*scene_num+1,1);
%     
%     c1 = zeros(var_num_s1,1);
%     c2 = zeros(var_num_s2,1);
%     c3 = zeros(var_num_s3,1);
%     
%     L1 = zeros(1,var_num_s1);
%     L2 = zeros(1,var_num_s2);
%     L3 = zeros(1,var_num_s3);
% 
%     % Initiate Yields
%     yield = zeros(scene_num,3);
%     yield(1,:)=[3.0, 3.6, 24];
%     yield(2,:)=[2.5, 3.0, 20];
%     yield(3,:)=[2.0, 2.4, 16];
%     
%     for i = 1:scene_num
%     A1(i,1)=yield(i,1);
%     A1(i,2*i:(2*i+1))=[1,-1];
%     A1(i,2*scene_num+1+i)=-1;
%     end
%     A1(scene_num+1,1)=1;
%     A1(scene_num+1,3*scene_num+2)=1;
%     
%     b1(1:scene_num) = 200;
%     b1(scene_num+1) = 500;
%     c1(1) = 150;
%     for i = 1:scene_num
%     c1(2*i:(2*i+1))=1/scene_num*[238;-170];
%     end
%     
%     % Subproblem 2 variables: x_2, y_2i, w_2i 
%     % Subproblem 2 slack variables: z_2i
%     % numbers = 3*scene_num + 1
%     for i = 1:scene_num
%     A2(i,1)=yield(i,2);
%     A2(i,2*i:(2*i+1))=[1,-1];
%     A2(i,2*scene_num+1+i)=-1;
%     end
%     A2(scene_num+1,1)=1;
%     A2(scene_num+1,3*scene_num+2)=1;
%     
%     b2(1:scene_num) = 240;
%     b2(scene_num+1) = 500;
%     c2(1) = 230;
%     for i = 1:scene_num
%     c2(2*i:(2*i+1))=1/scene_num*[210;-150];
%     end
%     
%     
%     % Subproblem 3 variables: x_3, w_3i, w_4i
%     % Subproblem 3 slack variables: z_3i, z_4i
%     % numbers = 4*scene_num + 2
%     for i = 1:scene_num
%     A3(i,1)=-yield(i,3);
%     A3(i,2*i:(2*i+1))=[1,1];
%     A3(i,2*scene_num+1+i)=1;
%     b3(i)=0;
%     end
%     for i = 1:scene_num
%     A3(i+scene_num,2*i)=1;
%     A3(i+scene_num,3*scene_num+1+i)=1;
%     b3(i+scene_num)=6000;
%     end
%     A3(2*scene_num+1,1)=1;
%     A3(2*scene_num+1,4*scene_num+2)=1;
%     b3(2*scene_num+1)=500;
%     
%     c3(1) = 260;
%     for i = 1:scene_num
%     c3(2*i:(2*i+1))=1/scene_num*[-36;-10];
%     end
%     
%     % master problem
%     L1(1) = 1;
%     L2(1) = 1;
%     L3(1) = 1;
%     b0 = 500;
%     compute_matrices=0;
% end
% 
% while (compute_initial_extreme_points)
% % Initial Extreme Points
% x1_init=zeros(var_num_s1,1);
% basic_var_1=[2,4,6,11];    % NEED TO CHANGE FOR MORE SCENES
% B1=A1(:,basic_var_1);
% x1_init_B=B1\b1;
% x1_init(basic_var_1)=x1_init_B;
% 
% x2_init=zeros(3*scene_num + 1,1);
% basic_var_2=[2,4,6,11];    % NEED TO CHANGE FOR MORE SCENES
% B2=A2(:,basic_var_1);
% x2_init_B=B2\b2;
% x2_init(basic_var_2)=x2_init_B;
% 
% 
% x3_init=zeros(4*scene_num + 1,1);
% basic_var_3=[3,5,7,11,12,13,14];    % need to change for scene_num
% B3=A3(:,basic_var_3);
% x3_init_B=B3\b3;
% x3_init(basic_var_3)=x3_init_B;
% compute_initial_extreme_points=0;
% end

%% Instance with Unbounded Subproblems
while (compute_matrices)
    scene_num = 3;
    epsilon = 1e-16;
    flag_optimal=0;
    flag_unbounded=0;
    
    % Initialization of Matrices and Vectors
    % Subproblem 1 variables: x_1, y_1i, w_1i
    % Subproblem 1 slack variables: z_1i
    % Subproblem 1 x_1 upper bound slack variable z1
    % numbers = 3*scene_num + 1
    var_num_s1=3*scene_num+1; 
    var_num_s2=3*scene_num+1;
    var_num_s3=4*scene_num+1; % without x_4
    
    A1 = zeros(scene_num,var_num_s1);
    A2 = zeros(scene_num,var_num_s2);
    A3 = zeros(2*scene_num,var_num_s3);
    
    b1 = zeros(scene_num,1);
    b2 = zeros(scene_num,1);
    b3 = zeros(2*scene_num,1);
    
    c1 = zeros(var_num_s1,1);
    c2 = zeros(var_num_s2,1);
    c3 = zeros(var_num_s3,1);
    
    L1 = zeros(1,var_num_s1);
    L2 = zeros(1,var_num_s2);
    L3 = zeros(1,var_num_s3);

    % Initiate Yields
    yield = zeros(scene_num,3);
    yield(1,:)=[3.0, 3.6, 24];
    yield(2,:)=[2.5, 3.0, 20];
    yield(3,:)=[2.0, 2.4, 16];
    
    for i = 1:scene_num
    A1(i,1)=yield(i,1);
    A1(i,2*i:(2*i+1))=[1,-1];
    A1(i,2*scene_num+1+i)=-1;
    end
    
    b1(1:scene_num) = 200;
    c1(1) = 150;
    for i = 1:scene_num
    c1(2*i:(2*i+1))=1/scene_num*[238;-170];
    end
    
    % Subproblem 2 variables: x_2, y_2i, w_2i 
    % Subproblem 2 slack variables: z_2i
    % numbers = 3*scene_num + 1
    for i = 1:scene_num
    A2(i,1)=yield(i,2);
    A2(i,2*i:(2*i+1))=[1,-1];
    A2(i,2*scene_num+1+i)=-1;
    end
    
    b2(1:scene_num) = 240;
    c2(1) = 230;
    for i = 1:scene_num
    c2(2*i:(2*i+1))=1/scene_num*[210;-150];
    end
    
    
    % Subproblem 3 variables: x_3, w_3i, w_4i
    % Subproblem 3 slack variables: z_3i, z_4i
    % numbers = 4*scene_num + 1
    for i = 1:scene_num
    A3(i,1)=-yield(i,3);
    A3(i,2*i:(2*i+1))=[1,1];
    A3(i,2*scene_num+1+i)=1;
    b3(i)=0;
    end
    for i = 1:scene_num
    A3(i+scene_num,2*i)=1;
    A3(i+scene_num,3*scene_num+1+i)=1;
    b3(i+scene_num)=6000;
    end
    
    c3(1) = 260;
    for i = 1:scene_num
    c3(2*i:(2*i+1))=1/scene_num*[-36;-10];
    end
    
    % master problem
    L1(1) = 1;
    L2(1) = 1;
    L3(1) = 1;
    b0 = 500;
    compute_matrices=0;
end

while (compute_initial_extreme_points)
% Initial Extreme Points
x1_init=zeros(var_num_s1,1);
basic_var_1=[2,4,6];    % NEED TO CHANGE FOR MORE SCENES
B1=A1(:,basic_var_1);
x1_init_B=B1\b1;
x1_init(basic_var_1)=x1_init_B;

x2_init=zeros(3*scene_num + 1,1);
basic_var_2=[2,4,6];    % NEED TO CHANGE FOR MORE SCENES
B2=A2(:,basic_var_1);
x2_init_B=B2\b2;
x2_init(basic_var_2)=x2_init_B;


x3_init=zeros(4*scene_num + 1,1);
basic_var_3=[3,5,7,11,12,13];    % need to change for scene_num
B3=A3(:,basic_var_3);
x3_init_B=B3\b3;
x3_init(basic_var_3)=x3_init_B;
compute_initial_extreme_points=0;
end



%% Construct Structure as the Input to Func
while (construct_structure)
    master.L{1}=L1;
    master.L{2}=L2;
    master.L{3}=L3;
    master.b=b0;
    sub.A{1}=A1;
    sub.A{2}=A2;
    sub.A{3}=A3;
    sub.b{1}=b1;
    sub.b{2}=b2;
    sub.b{3}=b3;
    sub.c{1}=c1;
    sub.c{2}=c2;
    sub.c{3}=c3;
    sub.v0{1}=x1_init;
    sub.v0{2}=x2_init;
    sub.v0{3}=x3_init;
    sub.basis{1} = basic_var_1;
    sub.basis{2} = basic_var_2;
    sub.basis{3} = basic_var_3;
    construct_structure=0;
end

%% Try to call function
[soln, fval, flag] = dw_dec(master, sub, 3)

%% DW-Decomposition Main Part
% Step 0: Generate initial matrices
% Generate inital basis B for Master problem
basic_var = [4,1,2,3];
B = eye(4);
b = [500;ones(scene_num,1)];
c = zeros(4,1);
sol.c{1} = [0];   % unclear
sol.x{1} = [500]; % unclear
sol.c{2} = c1;
sol.x{2} = x1_init;
sol.c{3} = c2;
sol.x{3} = x2_init;
sol.c{4} = c3;
sol.x{4} = x3_init;
sol.p{1} = 0;
sol.p{2} = 1;
sol.p{3} = 2;
sol.p{4} = 3;
f_B = zeros(scene_num+1,1);

while ((1-flag_optimal) && (1-flag_unbounded))
    unbnd_sub = 0; % -- Jason: a flag to indicate the existence of an unbounded sub_problem
    % Step 1: Calculate pie
    % second iteration 
    f_B(1) = sol.c{1}'*sol.x{1};
    f_B(2) = sol.c{2}'*sol.x{2};
    f_B(3) = sol.c{3}'*sol.x{3};
    f_B(4) = sol.c{4}'*sol.x{4};
    x_B = B\b;
    pie = B'\f_B;
    pie_1 = pie(1);    % Recall pie_1 and pie_2

    % Step 2: Optimality Check
    % Use my SimplexMethod instead of linprog, SimplexMethod works
    [x1_linp,fval1,exitflag1]=SimplexMethod(A1,b1,(c1'-pie_1'*L1)',basic_var_1');
    [x2_linp,fval2,exitflag2]=SimplexMethod(A2,b2,(c2'-pie_1'*L2)',basic_var_2');
    [x3_linp,fval3,exitflag3]=SimplexMethod(A3,b3,(c3'-pie_1'*L3)',basic_var_3');
    r_N = [0;fval1-pie(2);fval2-pie(3);fval3-pie(4)];
    sol.tempx{1}=x1_linp;
    sol.tempx{2}=x2_linp;
    sol.tempx{3}=x3_linp;
    entering_var = find(r_N==min(r_N)); % Very Important
    % -- Jason: the following is to evaluate the existence of unbounded sub_problem
    if exitflag1==-3
        enterring_var = 2;
        unbnd_sub = 1;
    end
    if exitflag2==-3
        enterring_var = 3;
        unbnd_sub = 1;
    end
    if exitflag3==-3
        enterring_var = 4;
        unbnd_sub = 1;
    end

    % use null(Ai) to find extreme rays -- Do I have to do that?
    if sum(r_N<-1e-3==1)==0 % if # of all r_N>>0, then stop
        flag_optimal=1;
        %break
    end

    % Step 3: calculate direction -- Jason
    a=zeros(scene_num+1,1);
    entering_prob = sol.p{entering_var};
    % The case of bounded subproblem
    if entering_prob>0
        if unbnd_sub == 0 % -- Jason
            a(entering_prob+1)=1;
        else
            a(entering_prob+1)=0;
        end
        a(1)=master.L{entering_prob}*sol.tempx{entering_prob};
    end
    % The case of unbounded subproblem
    if entering_prob==0 % Why this is true?
       a(1)=1;
    end

    d = B\(-a);

    if sum(d<-epsilon==1)==0 % if # of all d>>0, then stop
        flag_unbounded=1;
        %break
    end

    % Step 4: Minimum ratio test
    departing_candiates = find(d<-epsilon);
    ratio=-x_B(departing_candiates)./d(departing_candiates);
    alpha=min(ratio);
    departing_candidates = departing_candiates(find(ratio==alpha));
    departing_var = departing_candidates(1);
    %if (sum(departing_candidates==entering_var))
    %    departing_var=entering_var;
    %end


    % Step 5: Update all values
    B(:,departing_var)=a;

    if entering_prob>0
       f_B(departing_var) = sub.c{entering_prob}'*sol.tempx{entering_prob};
       sol.x{departing_var} = sol.tempx{entering_prob};
       sol.c{departing_var} = sub.c{entering_prob};
       sol.p{departing_var} = entering_prob;
    end

    if entering_prob==0
       f_B(departing_var)=B\b(departing_var);
       sol.c{departing_var} = [0];
       sol.p{departing_var} = 0;
    end



% if entering_var ==1
%     c_B(1) = 0;
% end
% if entering_var ==2
%     c_B(departing_var) = c1'*x1_linp;
%     x_B(entering_var)=1;
% end
% if entering_var ==3
%     c_B(departing_var) = c2'*x2_linp;
%     x_B(entering_var)=1;
% end
% if entering_var ==4
%     c_B(departing_var) = c3'*x3_linp;
%     x_B(entering_var)=1;
% end
end

% output results
val = f_B' * (B\b);
display(val)
% output varaibles
x{1}=zeros(var_num_s1,1);
x{2}=zeros(var_num_s2,1);
x{3}=zeros(var_num_s3,1);
lambda=B\b;
for i = 1:4
    prob = sol.p{i};
    x{prob} = x{prob} + lambda(i)*sol.x{i};
end    