% This fuction generate the Kolmogrov Matrix
% inputs: N1,N2,N3,maxValues
% Other parameters are hard-coded inside the function
function output = generate_sym_K(N1,N2,maxValues,params)

Smax = maxValues(1);
vmax = maxValues(2);

 Smin = 0;
 vmin = 0;

% define parameters( globalize so can be used globally)
global theta_v k_v sig_v mu_v;
global cor_Sv r q;


% theta_v = 0.04; k_v=2; sig_v=0.2; mu_v=0; 
% cor_Sv = -0.5;
% r=0.03;
% q=0.05;

theta_v = params.theta_v; k_v=params.k_v; sig_v=params.sig_v; mu_v=params.mu_v; 
cor_Sv = params.cor_Sv;
r=params.r;
q=params.q;


% define variable: steps and vector
hS = (Smax-Smin)/N1;
hv = (vmax-vmin)/N2; 
matrix_size = (N1+1)*(N2+1);
% generate a sparse square matrix containing only zero
K = sparse(matrix_size,matrix_size);

% coefficients for center points
coeff_inner = get_coeff(2,2);
coeff_inner = subs(coeff_inner,{'hS', 'hv'},{ hS, hv});
for i=(2:N1)
    for j=(2:N2)     
            S = (i-1) * hS;
            v = (j-1) * hv;
             % evaluate the coefficients
            % temp_coe = subs(coeff_inner,{'S','v', 'r'},{S,v,r});
            temp_coe = eval(coeff_inner);
            list = get_neighbours(i,j);
            K(locate_D(i,j),list)=temp_coe'; 
    end 
end

% coefficients for boundary points
 coeff_Smin = get_coeff(1,2);
 coeff_Smin = subs(coeff_Smin,{'hS', 'hv'},{ hS, hv});
 coeff_Smax = get_coeff(3,2);
 coeff_Smax = subs(coeff_Smax,{'hS', 'hv'},{ hS, hv });
 for j=(2:N2)        
            i = 1; S=0;
            v = (j-1) * hv; 
%             temp_coe = subs(coeff_Smin,{'S','v', 'r'},{S,v,r});
            temp_coe = eval(coeff_Smin); 
            list = get_neighbours(i,j);
            K(locate_D(i,j),list)=temp_coe';         
            i = N1+1; S=Smax;
            v = (j-1) * hv;  
%             temp_coe = subs(coeff_Smax,{'S','v', 'r'},{S,v,r});
            temp_coe = eval(coeff_Smax);
            list = get_neighbours(i,j);
            K(locate_D(i,j),list)=temp_coe'; 
 end 

 coeff_vmin = get_coeff(2,1);
 coeff_vmin = subs(coeff_vmin,{'hS', 'hv'},{ hS, hv });
 coeff_vmax = get_coeff(2,3);
 coeff_vmax = subs(coeff_vmax,{'hS', 'hv'},{ hS, hv });
 for i=(2:N1)    
            j = 1; v=0;
            S = (i-1) * hS;
%             temp_coe = subs(coeff_vmin,{'S','v', 'r'},{S,v,r});
            temp_coe = eval(coeff_vmin);
            list = get_neighbours(i,j);
            K(locate_D(i,j),list)=temp_coe'; 
            j = N2+1; v=vmax;
            S = (i-1) * hS;  
%             temp_coe = subs(coeff_vmax,{'S','v', 'r'},{S,v,r});
            temp_coe = eval(coeff_vmax);
            list = get_neighbours(i,j);
            K(locate_D(i,j),list)=temp_coe'; 
 end 
 

%  % coefficients on edges
 loc_edges = [1 1 ; 1 3; 3 1; 3 3 ];
     boundary_S = [0 0 hS*N1]; index_S = [1 0 N1+1];
     boundary_v = [0 0 hv*N2]; index_v = [1 0 N2+1]; 
 for index=1:size(loc_edges,1)
     loc9 = loc_edges(index,:);
     coeff_edge = get_coeff(loc9(1),loc9(2));
     coeff_edge = subs(coeff_edge,{'hS', 'hv'},{ hS, hv });
                 S = boundary_S(loc9(1));
                 v = boundary_v(loc9(2)); 
%                  temp_coe = subs(coeff_edge,{'S','v','r'},{S,v,r});
                 temp_coe = eval(coeff_edge);
                 list = get_neighbours(index_S(loc9(1)),index_v(loc9(2)));
                 K(locate_D(index_S(loc9(1)),index_v(loc9(2))),list)=temp_coe'; 
 end
          
          
location_index = zeros(matrix_size,2);
    for iii=1:N1+1
        for jjj=1:N2+1
              loc = locate_D(iii,jjj);
        end
    end

% coefficients for whatever that's left    
% D = sym('D', [matrix_size 1]);
% [~, list_b] = get_inner_boundary();
% for idx = list_b
%             [i,j,k] = unlocate(idx);
%             S = (i-1) * hS;
%             v = (j-1) * hv;
%             r = (k-1) * hr;
%             % discretization
%             second_S = derivative_S2(i,j,k);
%             second_v = derivative_v2(i,j,k);
%             second_r = derivative_r2(i,j,k);
%             first_S = derivative_S1(i,j,k);
%             first_v = derivative_v1(i,j,k);
%             first_r = derivative_r1(i,j,k);
%             second_Sv = derivative_Sv(i,j,k);
% %           second_Sr = derivative_Sr(i,j,k);
% %           second_vr = derivative_vr(i,j,k);
%             % matrix tryout
%             matrix_K = (v*S^2/2) * (second_S) + (sig_v^2*v/2)*(second_v) + (sig_r^2*r/2)*(second_r)...
%                      + (cor_Sv * sig_v *v*S) *(second_Sv)...
%                      + (k_r*(theta_r-r)-mu_r*r) *(first_r) + (r-q)*S * (first_S) + (k_v*(theta_v-v)-mu_v*v) * (first_v);   
%             %        + (cor_Sr * sig_r *sqrt(r)*S) *(second_Sr) +  (cor_vr * sig_v*sig_r* sqrt(v*r)) *(second_vr)...                   
%             K(locate_D(i,j,k),:)=get_coeffs_boundary(matrix_K,i,j,k); 
% end

output = (K);


%% Subfunctions
    function [list_inner, list_left] = get_inner_boundary()
      listResult=[];
      for idxS=[1 N1+1]
          for idxv = [1 N2+1]
              for idxr =[1 N3+1]
                  listResult=[listResult locate_D(idxS,idxv,idxr)];
              end
          end
      end
      list_left = listResult;
      list_inner=[];
    end

    function coeffs_in_K = get_coeffs_boundary(matrix_K,i,j,k)
        location_list = get_neighbours(i,j,k);
        coeffs_in_K=zeros(1,matrix_size);
        listsize=size(location_list,2);
        for list_idx=(1:listsize)
         coeffs_temp = coeffs(matrix_K,[D(location_list(list_idx))]);
         if size(coeffs_temp,2)~=1
             coeffs_in_K(location_list(list_idx)) = coeffs_temp(2);
         end
        end
    end

    function list=get_neighbours(i,j)
        list=[];
        for ii=max(1,i-1):min(N1+1,i+1)
            for jj=max(1,j-1):min(N2+1,j+1)
                    list=[list locate_D(ii,jj)];
            end
        end
    end

%   location mapping: 3D to 1D
    function location = locate_D(i,j)
     % mapping matric D(i,j,k) to D((N1+1)*(N2+1)*(N3+1))
     location = (i-1)*(N2+1)...
              + j;
     location_index(location,:)=[i j ];     
    end

%   location mapping: 1D to 3D
    function [ii,jj] = unlocate(loc)
      ii=location_index(loc,1);
      jj=location_index(loc,2);
    end

%% Derivatives discretization
    function  second_S=derivative_S2(i,j,k)
     if (i==N1+1)||(i==1)
        second_S =0;
     else 
        second_S = (D(locate_D(i+1,j,k)) - 2*D(locate_D(i,j,k)) + D(locate_D(i-1,j,k)))...
                 / (hS^2);
     end
    end

    function second_v=derivative_v2(i,j,k)
     if (j==1)||(j==N2+1)
        second_v = 0;
     else
        second_v = (D(locate_D(i,j+1,k)) - 2*D(locate_D(i,j,k)) + D(locate_D(i,j-1,k)))...
                 / (hv^2);    
     end 
    end

    function second_r=derivative_r2(i,j,k)
     if (k==1)||(k==N3+1)
        second_r = 0;
     else
        second_r = (D(locate_D(i,j,k+1)) - 2*D(locate_D(i,j,k)) + D(locate_D(i,j,k-1)))...
                 / (hr^2);    
     end 
    end

    function first_S = derivative_S1(i,j,k)
        if (i==1)
            first_S = ( D(locate_D(2,j,k)) - D(locate_D(1,j,k)) )...
                    / hS;
            return;
        end
        if (i==N1+1)
             first_S = ( D(locate_D(N1+1,j,k)) - D(locate_D(N1,j,k)) )...
                     / hS ;           
        else
             first_S = ( D(locate_D(i+1,j,k)) - D(locate_D(i-1,j,k)) )...
                     / hS / 2;         
        end
    end

    function first_v = derivative_v1(i,j,k)
        if (j==1)
            first_v = ( D(locate_D(i,2,k)) - D(locate_D(i,1,k)) )...
                    / hv;
            return;
        end
        if (j==N2+1)
             first_v = ( D(locate_D(i,N2+1,k)) - D(locate_D(i,N2,k)) )...
                     / hv ;           
        else
             first_v = ( D(locate_D(i,j+1,k)) - D(locate_D(i,j-1,k)) )...
                     / hv / 2;         
        end
    end

    function first_r = derivative_r1(i,j,k)
        if (k==1)
            first_r = ( D(locate_D(i,j,2)) - D(locate_D(i,j,1)) )...
                    / hr;
            return;
        end
        if (k==N3+1)
             first_r = ( D(locate_D(i,j,N3+1)) - D(locate_D(i,j,N3)) )...
                     / hr ;           
        else
             first_r = ( D(locate_D(i,j,k+1)) - D(locate_D(i,j,k-1)) )...
                     / hr / 2;         
        end
    end

    function second_Sv = derivative_Sv(i,j,k)
        if (i==1 || j==1) || (i==N1+1 || j==N2+1)
            second_Sv=0;
        else
            if cor_Sv>0
            second_Sv = 1/2 *( ( D(locate_D(i+1,j+1,k)) - D(locate_D(i,j+1,k)) - ( D(locate_D(i+1,j,k))-D(locate_D(i,j,k)) ) )/hS/hv...
                      + (D(locate_D(i,j,k))-D(locate_D(i-1,j,k))-(D(locate_D(i,j-1,k))-D(locate_D(i-1,j-1,k))))/hS/hv );
            else
            second_Sv = 1/2 *( ( D(locate_D(i,j+1,k)) - D(locate_D(i-1,j+1,k)) - ( D(locate_D(i,j,k))-D(locate_D(i-1,j,k)) ) )/hS/hv...
                      + (D(locate_D(i+1,j,k))-D(locate_D(i,j,k))-(D(locate_D(i+1,j-1,k))-D(locate_D(i,j-1,k))))/hS/hv );      
            end
        end
    end

%     function second_Sr = derivative_Sr(i,j,k)
%         if (i==1 || k==1) || (i==N1+1 || k==N3+1)
%             second_Sr=0;
%         else
%             if cor_Sr>0
%             second_Sr = 1/2 *( ( D(locate_D(i+1,j,k+1)) - D(locate_D(i,j,k+1)) - ( D(locate_D(i+1,j,k))-D(locate_D(i,j,k)) ) )/hS/hr...
%                       + (D(locate_D(i,j,k))-D(locate_D(i-1,j,k))-(D(locate_D(i,j,k-1))-D(locate_D(i-1,j,k-1))))/hS/hr );
%             else
%             second_Sr = 1/2 *( ( D(locate_D(i,j,k+1)) - D(locate_D(i-1,j,k+1)) - ( D(locate_D(i,j,k))-D(locate_D(i-1,j,k)) ) )/hS/hr...
%                       + (D(locate_D(i+1,j,k))-D(locate_D(i,j,k))-(D(locate_D(i+1,j,k-1))-D(locate_D(i,j,k-1))))/hS/hr );      
%             end
%         end
%     end

%     function second_vr = derivative_vr(i,j,k)
%         if (j==1 || k==1) || (j==N2+1 || k==N3+1)
%             second_vr=0;
%         else
%             if cor_vr>0
%             second_vr = 1/2 *( ( D(locate_D(i,j+1,k+1)) - D(locate_D(i,j,k+1)) - ( D(locate_D(i,j+1,k))-D(locate_D(i,j,k)) ) )/hv/hr...
%                       + (D(locate_D(i,j,k))-D(locate_D(i,j-1,k))-(D(locate_D(i,j,k-1))-D(locate_D(i,j-1,k-1))))/hv/hr );
%             else
%             second_vr = 1/2 *( ( D(locate_D(i,j,k+1)) - D(locate_D(i,j-1,k+1)) - ( D(locate_D(i,j,k))-D(locate_D(i,j-1,k)) ) )/hv/hr...
%                       + (D(locate_D(i,j+1,k))-D(locate_D(i,j,k))-(D(locate_D(i,j+1,k-1))-D(locate_D(i,j,k-1))))/hv/hr );      
%             end
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%% get coefficient of "center" points 
    % (i,j,k) is an abstract way of representing the location in grid
    %  E.g: i=1 is lower boundary along S direction
    %       i=2 is centering points regarding S direction
    %       i=3 is upper boundary along S direction
    % (2,2,2) is for all points not on any boudary
    function output = get_coeff(i,j)
      
        n1=2; n2=2;
        global theta_v k_v sig_v mu_v;
        global cor_Sv r q;
        
        S=sym('S'); v = sym('v'); 
        hS = sym('hS'); hv = sym('hv');     
        location_index = zeros(9,2);
        D = sym('D', [9 1]);
        
        % initialize location_index.
            for iii=1:3
                for jjj=1:3
                        loc = locate_D(iii,jjj);
                end
            end   

         output = calculate_coefficients(i,j);
%        save('coeff27.mat','output');
      
%% sub functions        
    function result = calculate_coefficients(i_input,j_input)
        i=i_input; j=j_input; 
        % discretization
        second_S=derivative_S2;        second_v=derivative_v2;       
        first_S = derivative_S1;        first_v = derivative_v1;  
        second_Sv = derivative_Sv;%       second_Sr = derivative_Sr;%       second_vr = derivative_vr;
        % matrix tryout
        matrix_K = (v*S^2/2) * (second_S) + (sig_v^2*v/2)*(second_v) ...
                 + (cor_Sv * sig_v *v*S) *(second_Sv) ...
                 + (r-q)*S * (first_S) + (k_v*(theta_v-v)-mu_v*v) * (first_v);

        coeffs_in_K = sym('coeffs_in_K',[9 1]);           
        for idx=(1:9)
            coeffs_temp = coeffs(matrix_K,[D(idx)]);
            if size(coeffs_temp,2)~=1
                coeffs_in_K(idx) = coeffs_temp(2);
            end
        end
        % remove zeros:
        neibourlist = get_neighbours(i,j);
        coeffs_in_K = coeffs_in_K(neibourlist,:);
        list =cell(9,1);
        for i = 1:9
            list{i} = ['coeffs_in_K' num2str(i)];
        end
        result = subs(coeffs_in_K,list,{zeros(9,1)});
    end
    
    % this function returns all the neighbours of (i,j,k)
    function list=get_neighbours(i,j)
        list=[];
        for ii=max(1,i-1):min(3,i+1)
            for jj=max(1,j-1):min(3,j+1)
                    list=[list locate_D(ii,jj)];
            end
        end
    end

    % location mapping: 3D to 1D
    function location = locate_D(i,j)
     % mapping matric D(i,j,k) to D((n1+1)*(n2+1)*(N3+1))
     location = (i-1)*(3)...
              + j;
     location_index(location,:)=[i j];     
    end


%% Derivatives discretization
    function  second_S=derivative_S2()
     if (i==3)||(i==1)
        second_S = 0;%(D(locate_D(3,j)) - 2*D(locate_D(2,j)) + D(locate_D(1,j)))/ (hS^2);
     else 
        second_S = (D(locate_D(i+1,j)) - 2*D(locate_D(i,j)) + D(locate_D(i-1,j)))...
                 / (hS^2);
     end
    end

    function second_v=derivative_v2()
     if (j==1)||(j==3)
        second_v = 0;
     else
        second_v = (D(locate_D(i,j+1)) - 2*D(locate_D(i,j)) + D(locate_D(i,j-1)))...
                 / (hv^2); 
     end
    end

    function first_S = derivative_S1()
        if (i==1)
            first_S = ( D(locate_D(2,j)) - D(locate_D(1,j)) )/ hS;
            return;
        end
        if (i==3)
             first_S = ( D(locate_D(n1+1,j)) - D(locate_D(n1,j)) ) / hS ;           
        else
             first_S = ( D(locate_D(i+1,j)) - D(locate_D(i-1,j)) )...
                     / hS / 2;  
        end
    end

    function first_v = derivative_v1()
        if (j==1)
            first_v = ( D(locate_D(i,2)) - D(locate_D(i,1)) )...
                    / hv;
            return;
        end
        if (j==3)
             first_v = ( D(locate_D(i,n2+1)) - D(locate_D(i,n2)) )...
                     / hv ;           
        else
             first_v = ( D(locate_D(i,j+1)) - D(locate_D(i,j-1)) )...
                     / hv / 2;  
        end
    end

    function second_Sv = derivative_Sv()
        if (i==1 || j==1) || (i==3 || j==3)
            second_Sv=0;
        else
            if cor_Sv>0
            second_Sv = 1/2 *( ( D(locate_D(i+1,j+1)) - D(locate_D(i,j+1)) - ( D(locate_D(i+1,j))-D(locate_D(i,j)) ) )/hS/hv...
                      + (D(locate_D(i,j))-D(locate_D(i-1,j))-(D(locate_D(i,j-1))-D(locate_D(i-1,j-1))))/hS/hv );
            else
            second_Sv = 1/2 *( ( D(locate_D(i,j+1)) - D(locate_D(i-1,j+1)) - ( D(locate_D(i,j))-D(locate_D(i-1,j)) ) )/hS/hv...
                      + (D(locate_D(i+1,j))-D(locate_D(i,j))-(D(locate_D(i+1,j-1))-D(locate_D(i,j-1))))/hS/hv );      
            end
        end
    end


    end
 