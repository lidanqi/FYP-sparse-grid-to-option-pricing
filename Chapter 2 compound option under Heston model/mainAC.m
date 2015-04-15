% this is a trial of representing matrix form:
% result would be D(0): daughter option price at time zero
% adjusted for boundaries
%% parameters:
function [resultD,resultM,est,conds,totalTime] = mainAC(S0,NT,varargin)
tic
initialTime = cputime;
if (mod(NT,2)~=0)
     fprintf('timesteps must be even number');
     resultD=0;resultM=0;est=0;conds=0;
     return;
 end

 c1 = 4;
 c2 = 4;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 v0 = 0.04;
 r0 = 0.03;
 q = 0.05;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 max_itr  = 1000;
 condition_number = zeros(4,2);
 p = inputParser;
 p.addParamValue('level', [2 2]);
 p.addParamValue('w', 1.2);
 p.addParamValue('cor',0.5);
 % THIS TYPE DETERMINES EURO VS AME
 % EURO: 1;   AMERICAN:2
 p.addParamValue('type',1); 
 p.parse(varargin{:});
% ----------------------------------------------------------------------- 
% type
global OptionType
OptionType = p.Results.type;
%------------------------------------------------------------------------
% levels 
l = p.Results.level;
l1=l(1); l2=l(2); 
% define levels, number of segments:
N1= c1*2^l1; N2=c2*2^l2; 
% over-relaxation parameter
w = p.Results.w;
%------------------------------------------------------------------------
% strike price, time to maturity
TD=1; KD=1;
TM=0.6; KM=0.07;
% define limits 
Smax = 8; Smin = 0;
vmax = 0.5; vmin = 0;

% define variable: steps and vector
hS = (Smax-Smin)/N1; Svec = Smin:hS:Smax;
hv = (vmax-vmin)/N2; vvec = vmin:hv:vmax;
tau = TD/NT;
matrix_size = (N1+1)*(N2+1);
location_index = zeros(matrix_size,2);

% initialize location_index: map 3D to 1D
for iii=1:N1+1
    for jjj=1:N2+1
          loc = locate_D(iii,jjj);
    end
end




parameters = struct;
parameters.theta_v = 0.04; parameters.k_v=2; parameters.sig_v=0.2; parameters.mu_v=0; 
parameters.cor_Sv = p.Results.cor;
parameters.r =r0;
parameters.q =q;
% K matrix( Komogrov Operator): used for both daughter and mother
K = generate_sym_K(N1,N2,[Smax,vmax],parameters);
% save('K.mat', 'K');      
% loadingK = load('K.mat');
% K= loadingK.K;
R = r0 * ones(matrix_size, 1);
A = sparse(K-diag(R));
fprintf('Obtaining K,A completed, run time: %f s',toc);
start_time = cputime;

%% solve the daughter option problem
% ( I + theta * tau * A) * D(l+1) = ( I - (1-theta) * tau * A) * D(l) 
% terminal condition
global payoffVect;
payoffVect = computePayoff;
D_temp = payoffVect;
theta=1;
I = eye(matrix_size);
B   =  sparse(I - theta * tau * A);
RHS_temp = ( I + (1-theta) * tau * A) ;

DD = zeros(matrix_size,NT+1);
DD(:,NT+1) = D_temp;
D_new = zeros(matrix_size,1);

for l=1:NT
    if (l==4)
         theta=0.5;
         B   =  sparse(I - theta * tau * A);
         RHS_temp = ( I + (1-theta) * tau * A) ;
    end  
    % right hand side: known
    RHS = sparse(RHS_temp * D_temp);
    % RHS =  sparse( ( I + (1-theta) * tau * A) * D_temp);  
    % define for boundaries in S direction   
    [LHS,RHS] = define_boundaries_A(B,RHS,l,1);
    % ==========================================
    % Solve LHS * D_new = RHS
    if (l==1 || l==4)
     condition_number(l,:)=[l condest(LHS)];
    % save(['LHS_D_' int2str(l) '.mat'],'LHS');
    end
        lower_LHS = -tril(LHS,-1);
        upper_LHS = -triu(LHS,1);
        diag_LHS = LHS + lower_LHS + upper_LHS;
        [x,m]=solut_PSOR_D(diag_LHS,lower_LHS,upper_LHS,D_temp,RHS,w,max_itr);
    D_new = x;
    if m==max_itr 
       % if not converging, quit the program
       fprintf('\nstep %d, not converging, quit',l);
       resultD=0;resultM=0;est=0;conds=0;totalTime=0;
       return;
       break
    end
    if m>max_itr/2 || mod(l,floor(NT/5))==0
       % constantly printing iterations
       % so that I know it's still running :p
       fprintf('\nTime: %d ;Iterations: %d',l,m);
    end
    DD(:,NT-l+1)=D_new;
    D_temp = D_new;
end

conds = condition_number([1 4],:);

% daughter option price at maturity for mother option
nt=NT/5*3;
D_TM = DD(:,nt+1);

%% solve the mother option problem
M_temp = max(D_TM-KM,0);
theta=1;
B   =  sparse(I - theta * tau * A);
RHS_temp = ( I + (1-theta) * tau * A) ;
MM = zeros(matrix_size,nt+1);
MM(:,nt+1) = M_temp;
M_new = zeros(matrix_size,1);
fprintf('\n-------------------------');
for l=1:nt
    if (l==4)
         theta=0.5;
         B   =  sparse(I - theta * tau * A);
         RHS_temp = ( I + (1-theta) * tau * A) ;
    end
    % right hand side: known
    RHS = sparse(RHS_temp * M_temp);
    % define for boundaries in S direction   
    [LHS,RHS] = define_boundaries_A(B,RHS,l,2);
    % ==========================================
    % compute exercise value for mother options first: payoffVect_M
    payoffVect_M = max(0,DD(:,nt-l+1)-KM);
    % Solve LHS * D_new = RHS
    if (l==1 || l==4)
     condition_number(l+3,:)=[l*100 condest(LHS)];
   %  save(['LHS_M_' int2str(l) '.mat'],'LHS');
    end
        lower_LHS = -tril(LHS,-1);
        upper_LHS = -triu(LHS,1);
        diag_LHS = LHS + lower_LHS + upper_LHS;
        [x,m]=solut_PSOR_M(diag_LHS,lower_LHS,upper_LHS,M_temp,RHS,w,max_itr,payoffVect_M);
    M_new = x;
    if m==max_itr 
       % if not converging, quit the program
       fprintf('\nstep %d, not converging, quit',l);
       break
    end
    if m>max_itr/2 || mod(l,floor(NT/5))==0
       % constantly printing iterations
       % so that I know it's still running :p
       fprintf('\nTime: %d ;Iterations: %d',l,m);
    end
    MM(:,nt-l+1)=M_new;
    M_temp = M_new;
end

% interpolation to get results
est_D=interpolation(Svec',vvec',DD(:,1),S0,v0);
est_M=interpolation(Svec',vvec',MM(:,1),S0,v0);
est=[est_D,est_M];
resultD = DD;
resultM = MM;
% display information
fprintf('\nlevels: %2.0f %2.0f , grid size: %d * %d = %d\n',l1,l2,N1+1,N2+1,matrix_size);
fprintf('Variance: %2.2f Interest: %2.2f\n',v0,r0);
fprintf('Input Stock  Price: '); fprintf('%6g ',S0);fprintf('\n');
fprintf('Estid D Option Price: '); fprintf('%6g ',est_D); fprintf('\n');
fprintf('Estid M Option Price: '); fprintf('%6g ',est_M); fprintf('\n');
fprintf('run time: %f seconds\n\n',cputime-start_time); 

totalTime=cputime-initialTime;

%% sub functions

    % this is to define some of the diagonal elements on boundary points
    % in S direction
    function [LHS_new,RHS_new] = define_boundaries_A(LHS,RHS,l,mode)
        % mode = 1: daughter option
        % mode = 2: mother option
        % current time:[  l   * tau ] to maturity
        % next time:   [ (l+1)* tau ] to maturity
        if N1>1
           % i = 1,N1+1
           idx_Smin =  0*(N2+1)+1 : 1*(N2+1);
           idx_Smax = N1*(N2+1)+1 : matrix_size;
           LHS(idx_Smin,idx_Smin) = eye((N2+1));        
           LHS(idx_Smax,idx_Smax) = eye((N2+1));    
        if (mode==1)
                   % i = 1, S = 0 : daughter value is discounted strike
                   for idx = idx_Smax
                       interest = r0;
                       RHS(idx) = Smax-KD*exp(- (interest) * (l*tau));
                   end   
                   % i = 1, S = Smax: daughter option value is 0 
                   RHS(idx_Smin) = 0;        
        else
                   % i = 1, S = 0 : daughter value is discounted strike >km
                   %                mother value is 0
                   RHS(idx_Smin)=0;
                   % i = 1, S = Smax:  daughter option value is 0;
                   %                   mother option value is discounted strike KM
                   for idx = idx_Smax
                       interest = r0;
                       RHS(idx) = (DD(idx,l)-KM)*exp(- (interest) * (l*tau));
                   end          
        end

        end
        LHS_new = LHS; 
        RHS_new = RHS;  
    end


%   This function defines the payoff values of all points: (N1+1)*(N2+1)*(N3+1)
    function prices = computePayoff()
    % put option : max(0,KD-stock_price);
     prices = zeros(matrix_size,1);
        for index=1:N1+1
            index_i = (index-1)*(N2+1)+1 : (index)*(N2+1); 
            Stock = Smin+(index-1) * hS;
            prices(index_i) =  max(0,Stock-KD);
        end
    end

%% mapping
%   location mapping: 3D to 1D
    function location = locate_D(i,j)
     % mapping matric D(i,j,k) to D((N1+1)*(N2+1)*(N3+1))
     location = (i-1)*(N2+1)...
              + j;
     location_index(location,:)=[i j];     
    end

%   location mapping: 1D to 3D
    function [ii,jj] = unlocate(loc)
      ii=location_index(loc,1);
      jj=location_index(loc,2);
    end

end


%% PSOR daughter
    % This function computes LHS*x = RHS. 
    % left hand side: d,l,u
    % right hand side: b
    % relaxation factor: W
    function [x,m]=solut_PSOR_D(d,l,u,x0,b,W,max_itr)
        global payoffVect
        m=0; %number of steps taken in the SOR method
        iteration_condition = 1; 
        x=0;
        global OptionType
       if (OptionType==2)
       try
        while iteration_condition && m<max_itr % do the iteration
           x1=max((d-W*l)\(((1-W)*d+W*u)*x0+W*b),payoffVect);
           if norm(x1-x0,inf) < 1*10^(-6) %iteration condition
               iteration_condition = 0;
               x = x1;
               m = m+1; %count the iteration steps
           else
               x0 = x1;
               m = m+1 ;%count the iteration steps
           end
        end 
        if x==0
            x = x1;
        end
       catch exception
           % do nothing
       end
       else
       try
        while iteration_condition && m<max_itr % do the iteration
           x1=(d-W*l)\(((1-W)*d+W*u)*x0+W*b);;
           if norm(x1-x0,inf) < 1*10^(-6) %iteration condition
               iteration_condition = 0;
               x = x1;
               m = m+1; %count the iteration steps
           else
               x0 = x1;
               m = m+1 ;%count the iteration steps
           end
        end 
        if x==0
            x = x1;
        end
       catch exception
           % do nothing
       end    
       end
    end
    
    
    %% PSOR mother
    % This function computes LHS*x = RHS. 
    % left hand side: d,l,u
    % right hand side: b
    % relaxation factor: W
    function [x,m]=solut_PSOR_M(d,l,u,x0,b,W,max_itr,payoffVect_M)
        m=0; %number of steps taken in the SOR method
        iteration_condition = 1; 
        x=0;
        global OptionType
    if (OptionType==2)
       try
        while iteration_condition && m<max_itr % do the iteration
            x1=max((d-W*l)\(((1-W)*d+W*u)*x0+W*b),payoffVect_M);
           if norm(x1-x0,inf) < 1*10^(-6) %iteration condition
               iteration_condition = 0;
               x = x1;
               m = m+1; %count the iteration steps
           else
               x0 = x1;
               m = m+1 ;%count the iteration steps
           end
        end 
        if x==0
            x = x1;
        end
        catch exception
           % do nothing
       end
    else
        try
        while iteration_condition && m<max_itr % do the iteration
            x1=((d-W*l)\(((1-W)*d+W*u)*x0+W*b));
           if norm(x1-x0,inf) < 1*10^(-6) %iteration condition
               iteration_condition = 0;
               x = x1;
               m = m+1; %count the iteration steps
           else
               x0 = x1;
               m = m+1 ;%count the iteration steps
           end
        end 
        if x==0
            x = x1;
        end
        catch exception
           % do nothing
       end 
    end
    end





