% this is a trial of representing matrix form:
% result would be D(0): daughter option price at time zero
% adjusted for boundaries
%% parameters:
function [resultD,est,conds,totalTime] = mainLB(NT,varargin)
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
 S0=10;
 m0=[6, 7, 8, 9, 10];
 x0=log(m0/S0);
 v0 = 0.5;
 r0 = 0.05;
 q = 0.0;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 max_itr  = 1000;
 condition_number = zeros(4,2);
 p = inputParser;
 p.addParamValue('level', [2 2]);
 p.addParamValue('w', 1.2);
 p.addParamValue('cor',0);
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
TD=1; 
% define limits 
xmax = 0; xmin = -4;
vmax = 2; vmin = 0;

% define variable: steps and vector
hx = (xmax-xmin)/N1; xvec = xmin:hx:xmax;
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
parameters.theta_v = 0.5; parameters.k_v=0.1; parameters.sig_v=0.25; parameters.mu_v=0; 
parameters.cor_Sv = p.Results.cor;
parameters.r =r0;
parameters.q =q;
% K matrix( Komogrov Operator): used for both daughter and mother
K = generate_sym_K(N1,N2,[xmin,vmax],parameters);
% save('K.mat', 'K');      
% loadingK = load('K.mat');
% K= loadingK.K;
Q = q * ones(matrix_size, 1);
A = sparse(K-diag(Q));
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
% x0=0; % log(m/S0)
% interpolation to get results
est_D=interpolation(xvec',vvec',DD(:,1),x0,v0);
est=[est_D]*S0; %transform back
resultD = DD;
% display information
fprintf('\nlevels: %2.0f %2.0f , grid size: %d * %d = %d\n',l1,l2,N1+1,N2+1,matrix_size);
fprintf('Variance: %2.2f Interest: %2.2f\n',v0,r0);
fprintf('Input Stock  Price: '); fprintf('%6g ',S0);fprintf('\n');
fprintf('Estid Option Price: '); fprintf('%6g ',est); fprintf('\n');
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
%         if N1>1
%            % i = 1,N1+1
            idx_Smin =  0*(N2+1)+1 : 1*(N2+1);
%            idx_Smax = N1*(N2+1)+1 : matrix_size;
%            LHS(idx_Smin,idx_Smin) = eye((N2+1));        
%            LHS(idx_Smax,idx_Smax) = eye((N2+1));         
%                    % i = 1, x = -4? : daughter value is discounted strike
                   for idx = idx_Smin
                       LHS(idx,:)=0;
                       LHS(idx,idx)=1;
                       interest = r0;
                       RHS(idx) = 1- exp(xmin)*exp(- (interest) * (l*tau));
                   end   
%                    % i = 1, S = Smax: daughter option value is 0 
%                    RHS(idx_Smax) = 0;        
%         end

%         idx_Smax=N1*(N2+1)+1 : matrix_size;
%         idx_Smax_1=(N1-1)*(N2+1)+1 : N1*(N2+1);
%         idx_Smax_2=(N1-2)*(N2+1)+1 : (N1-1)*(N2+1);
%         idx_Smax_3=(N1-3)*(N2+1)+1 : (N1-2)*(N2+1);
%         idx_Smax_4=(N1-4)*(N2+1)+1 : (N1-3)*(N2+1);
%         for row_idx=1:N2+1
%             LHS(idx_Smax(row_idx),:)=0;
%             LHS(idx_Smax(row_idx),idx_Smax(row_idx))= -25/12;
%             LHS(idx_Smax(row_idx),idx_Smax_1(row_idx))= 4;
%             LHS(idx_Smax(row_idx),idx_Smax_2(row_idx))= - 3;
%             LHS(idx_Smax(row_idx),idx_Smax_3(row_idx))= 4/3;
%             LHS(idx_Smax(row_idx),idx_Smax_4(row_idx))= -1/4;
%             RHS(idx_Smax(row_idx))=0;
%         end

        LHS_new = LHS; 
        RHS_new = RHS;  
    end


%   This function defines the payoff values of all points: (N1+1)*(N2+1)*(N3+1)
    function prices = computePayoff()
    % put option : max(0,KD-stock_price);
     prices = zeros(matrix_size,1);
        for index=1:N1+1
            index_i = (index-1)*(N2+1)+1 : (index)*(N2+1); 
            x = xmin+(index-1) * hx;
            prices(index_i) =  max(0,1-exp(x));
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
    
    
   


