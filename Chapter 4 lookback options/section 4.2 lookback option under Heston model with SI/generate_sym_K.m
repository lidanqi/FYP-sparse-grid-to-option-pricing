% This fuction generate the Kolmogrov Matrix
% inputs: N1,N2,N3,maxValues
% Other parameters are hard-coded inside the function
function output = generate_sym_K(N1,N2,N3,maxValues)

xmin = maxValues(1);
vmax = maxValues(2);
rmax = maxValues(3); 

 xmax = 0;
 vmin = 0;
 rmin = 0;
% define parameters( globalize so can be used globally)
global theta_v k_v sig_v mu_v;
global theta_r k_r sig_r mu_r;
global cor_Sv q;

theta_v = 0.02; k_v=1.5; sig_v=0.15; 
theta_r = 0.04; k_r=0.3; sig_r=0.1;  
cor_Sv = -0.5;
q=0;

% define variable: steps and vector
hx = (xmax-xmin)/N1;
hv = (vmax-vmin)/N2; 
hr = (rmax-rmin)/N3;
matrix_size = (N1+1)*(N2+1)*(N3+1);
% generate a sparse square matrix containing only zero
K = sparse(matrix_size,matrix_size);

% coefficients for center points
coeff_inner = get_coeff(2,2,2);
coeff_inner = subs(coeff_inner,{'hx', 'hv', 'hr'},{ hx, hv, hr });
for i=(2:N1)
    for j=(2:N2)
        for k=(2:N3)           
            x = (i-1) * hx;
            v = (j-1) * hv;
            r = (k-1) * hr;
            % evaluate the coefficients
            % temp_coe = subs(coeff_inner,{'x','v', 'r'},{x,v,r});
            temp_coe = eval(coeff_inner);
            list = get_neighbours(i,j,k);
            K(locate_D(i,j,k),list)=temp_coe'; 
        end
    end 
end

% coefficients for boundary points
 coeff_xmin = get_coeff(1,2,2);
 coeff_xmin = subs(coeff_xmin,{'hx', 'hv', 'hr'},{ hx, hv, hr });
 coeff_xmax = get_coeff(3,2,2);
 coeff_xmax = subs(coeff_xmax,{'hx', 'hv', 'hr'},{ hx, hv, hr });
 for j=(2:N2)
        for k=(2:N3)           
            i = 1; x=0;
            v = (j-1) * hv;  r = (k-1) * hr;
%             temp_coe = subs(coeff_Smin,{'x','v', 'r'},{x,v,r});
            temp_coe = eval(coeff_xmin); 
            list = get_neighbours(i,j,k);
            K(locate_D(i,j,k),list)=temp_coe';         
            i = N1+1; x=xmax;
            v = (j-1) * hv;  r = (k-1) * hr;
%             temp_coe = subs(coeff_Smax,{'x','v', 'r'},{x,v,r});
            temp_coe = eval(coeff_xmin);
            list = get_neighbours(i,j,k);
            K(locate_D(i,j,k),list)=temp_coe'; 
        end
 end 

 coeff_vmin = get_coeff(2,1,2);
 coeff_vmin = subs(coeff_vmin,{'hx', 'hv', 'hr'},{ hx, hv, hr });
 coeff_vmax = get_coeff(2,3,2);
 coeff_vmax = subs(coeff_vmax,{'hx', 'hv', 'hr'},{ hx, hv, hr });
 for i=(2:N1)
        for k=(2:N3)           
            j = 1; v=0;
            x = (i-1) * hx;  r = (k-1) * hr;
%             temp_coe = subs(coeff_vmin,{'x','v', 'r'},{x,v,r});
            temp_coe = eval(coeff_vmin);
            list = get_neighbours(i,j,k);
            K(locate_D(i,j,k),list)=temp_coe'; 
            j = N2+1; v=vmax;
            x = (i-1) * hx;  r = (k-1) * hr;
%             temp_coe = subs(coeff_vmax,{'x','v', 'r'},{x,v,r});
            temp_coe = eval(coeff_vmax);
            list = get_neighbours(i,j,k);
            K(locate_D(i,j,k),list)=temp_coe'; 
        end
 end 
 
 coeff_rmin = get_coeff(2,2,1 );
 coeff_rmin = subs(coeff_rmin,{'hx', 'hv', 'hr'},{ hx, hv, hr });
 coeff_rmax = get_coeff(2,2,3 );
 coeff_rmax = subs(coeff_rmax,{'hx', 'hv', 'hr'},{ hx, hv, hr });
 for i=(2:N1)
        for j=(2:N2)           
            % discrete value at point (i,j,k)
            k = 1; r=0;
            x = (i-1) * hx;  v = (j-1) * hv;
            % evaluate the coefficients
%             temp_coe = subs(coeff_rmin,{'x','v', 'r'},{x,v,r});
            temp_coe = eval(coeff_rmin);
            list = get_neighbours(i,j,k);
            K(locate_D(i,j,k),list)=temp_coe'; 
            
            k = N3+1; r=rmax;
            x = (i-1) * hx;  v = (j-1) * hv;
            % evaluate the coefficients
%             temp_coe = subs(coeff_rmax,{'x','v', 'r'},{x,v,r});
            temp_coe = eval(coeff_rmax);
            list = get_neighbours(i,j,k);
            K(locate_D(i,j,k),list)=temp_coe'; 
        end
 end 

%  % coefficients on edges
 loc_edges = [1 1 2; 1 3 2; 3 1 2; 3 3 2;...
              2 1 1; 2 1 3; 2 3 1; 2 3 3;...
              1 2 1; 1 2 3; 3 2 1; 3 2 3];
     boundary_x = [0 0 hx*N1]; index_S = [1 0 N1+1];
     boundary_v = [0 0 hv*N2]; index_v = [1 0 N2+1]; 
     boundary_r = [0 0 hr*N3]; index_r = [1 0 N3+1];
 for index=1:size(loc_edges,1)
     loc27 = loc_edges(index,:);
     coeff_edge = get_coeff(loc27(1),loc27(2),loc27(3));
     coeff_edge = subs(coeff_edge,{'hx', 'hv', 'hr'},{ hx, hv, hr });
     idx_center = find(loc27==2);
     switch (idx_center)
         case 1
             for subidx = 2:N1
                 x = (subidx-1) * hx;
                 v = boundary_v(loc27(2)); r=boundary_r(loc27(3));
%                  temp_coe = subs(coeff_edge,{'x','v','r'},{x,v,r});
                 temp_coe = eval(coeff_edge);
                 list = get_neighbours(subidx,index_v(loc27(2)),index_r(loc27(3)));
                 K(locate_D(subidx,index_v(loc27(2)),index_r(loc27(3))),list)=temp_coe'; 
             end
         case 2
             for subidx = 2:N2
                 v = (subidx-1) * hv;
                 x = boundary_x(loc27(1)); r=boundary_r(loc27(3));
%                  temp_coe = subs(coeff_edge,{'x','v','r'},{x,v,r});
                 temp_coe = eval(coeff_edge);
                 list = get_neighbours(index_S(loc27(1)),subidx,index_r(loc27(3)));
                 K(locate_D(index_S(loc27(1)),subidx,index_r(loc27(3))),list)=temp_coe'; 
             end
         case 3
             for subidx = 2:N3
                 r = (subidx-1) * hr;
                 x = boundary_x(loc27(1)); v=boundary_v(loc27(2));
%                  temp_coe = subs(coeff_edge,{'x','v','r'},{x,v,r});
                 temp_coe = eval(coeff_edge);
                 list = get_neighbours(index_S(loc27(1)),index_v(loc27(2)),subidx);
                 K(locate_D(index_S(loc27(1)),index_v(loc27(2)),subidx),list)=temp_coe'; 
             end
     end
 end
          
          
location_index = zeros(matrix_size,3);
    for iii=1:N1+1
        for jjj=1:N2+1
            for kkk=1:N3+1
              loc = locate_D(iii,jjj,kkk);
            end
        end
    end

% coefficients for whatever that's left    
D = sym('D', [matrix_size 1]);
[~, list_b] = get_inner_boundary();
for idx = list_b
            [i,j,k] = unlocate(idx);
            x = (i-1) * hx;
            v = (j-1) * hv;
            r = (k-1) * hr;
            % discretization
            second_x = derivative_x2(i,j,k);
            second_v = derivative_v2(i,j,k);
            second_r = derivative_r2(i,j,k);
            first_x = derivative_x1(i,j,k);
            first_v = derivative_v1(i,j,k);
            first_r = derivative_r1(i,j,k);
            second_xv = derivative_xv(i,j,k);
%           second_Sr = derivative_Sr(i,j,k);
%           second_vr = derivative_vr(i,j,k);
            % matrix tryout
            matrix_K = (v/2) * (second_x) + (sig_v^2*v/2)*(second_v) + (sig_r^2*r/2)*(second_r)...
                     - (cor_Sv * sig_v *v) *(second_xv)...
                     + (k_r*(theta_r-r))*(first_r) - ((r-q)+v/2)*(first_x) + (k_v*(theta_v-v) + cor_Sv*sig_v*v) * (first_v);   
            %        + (cor_Sr * sig_r *sqrt(r)*x) *(second_Sr) +  (cor_vr * sig_v*sig_r* sqrt(v*r)) *(second_vr)...                   
            K(locate_D(i,j,k),:)=get_coeffs_boundary(matrix_K,i,j,k); 
end

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

    function list=get_neighbours(i,j,k)
        list=[];
        for ii=max(1,i-1):min(N1+1,i+1)
            for jj=max(1,j-1):min(N2+1,j+1)
                for kk=max(1,k-1):min(N3+1,k+1)
                    list=[list locate_D(ii,jj,kk)];
                end
            end
        end
    end

%   location mapping: 3D to 1D
    function location = locate_D(i,j,k)
     % mapping matric D(i,j,k) to D((N1+1)*(N2+1)*(N3+1))
     location = (i-1)*(N2+1)*(N3+1)...
              + (j-1)*(N3+1)...
              + (k);
     location_index(location,:)=[i j k];     
    end

%   location mapping: 1D to 3D
    function [ii,jj,kk] = unlocate(loc)
      ii=location_index(loc,1);
      jj=location_index(loc,2);
      kk=location_index(loc,3);
    end

%% Derivatives discretization
    function  second_x=derivative_x2(i,j,k)
     if (i==N1+1)||(i==1)
         if (i==1)
             % error = O(h);
             second_x =0;
%              second_x=(D(locate_D(i+2,j,k)) - 2*D(locate_D(i+1,j,k)) + D(locate_D(i,j,k)))/ (hx^2);
         else
             second_x=(D(locate_D(i-2,j,k)) - 2*D(locate_D(i-1,j,k)) + D(locate_D(i,j,k)))/ (hx^2);
%              second_x =0;
         end
     else 
        second_x = (D(locate_D(i+1,j,k)) - 2*D(locate_D(i,j,k)) + D(locate_D(i-1,j,k)))...
                 / (hx^2);
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

    function first_x = derivative_x1(i,j,k)
        if (i==1)
            first_x = ( D(locate_D(2,j,k)) - D(locate_D(1,j,k)) )/ hx;
%             first_x=0;
            return;
        end
        if (i==N1+1)
%              first_x = ( D(locate_D(N1+1,j,k)) - D(locate_D(N1,j,k)) )...
%                      / hx ;           
            first_x=0;
        else
             first_x = ( D(locate_D(i+1,j,k)) - D(locate_D(i-1,j,k)) )...
                     / hx / 2;         
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

    function second_xv = derivative_xv(i,j,k)
        if (i==1 || j==1) || (i==N1+1 || j==N2+1)
            second_xv=0;
        else
            if cor_Sv>0
            second_xv = 1/2 *( ( D(locate_D(i+1,j+1,k)) - D(locate_D(i,j+1,k)) - ( D(locate_D(i+1,j,k))-D(locate_D(i,j,k)) ) )/hx/hv...
                      + (D(locate_D(i,j,k))-D(locate_D(i-1,j,k))-(D(locate_D(i,j-1,k))-D(locate_D(i-1,j-1,k))))/hx/hv );
            else
            second_xv = 1/2 *( ( D(locate_D(i,j+1,k)) - D(locate_D(i-1,j+1,k)) - ( D(locate_D(i,j,k))-D(locate_D(i-1,j,k)) ) )/hx/hv...
                      + (D(locate_D(i+1,j,k))-D(locate_D(i,j,k))-(D(locate_D(i+1,j-1,k))-D(locate_D(i,j-1,k))))/hx/hv );      
            end
        end
    end

%     function second_Sr = derivative_Sr(i,j,k)
%         if (i==1 || k==1) || (i==N1+1 || k==N3+1)
%             second_Sr=0;
%         else
%             if cor_Sr>0
%             second_Sr = 1/2 *( ( D(locate_D(i+1,j,k+1)) - D(locate_D(i,j,k+1)) - ( D(locate_D(i+1,j,k))-D(locate_D(i,j,k)) ) )/hx/hr...
%                       + (D(locate_D(i,j,k))-D(locate_D(i-1,j,k))-(D(locate_D(i,j,k-1))-D(locate_D(i-1,j,k-1))))/hx/hr );
%             else
%             second_Sr = 1/2 *( ( D(locate_D(i,j,k+1)) - D(locate_D(i-1,j,k+1)) - ( D(locate_D(i,j,k))-D(locate_D(i-1,j,k)) ) )/hx/hr...
%                       + (D(locate_D(i+1,j,k))-D(locate_D(i,j,k))-(D(locate_D(i+1,j,k-1))-D(locate_D(i,j,k-1))))/hx/hr );      
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
    %  E.g: i=1 is lower boundary along x direction
    %       i=2 is centering points regarding x direction
    %       i=3 is upper boundary along x direction
    % (2,2,2) is for all points not on any boudary
    function output = get_coeff(i,j,k)
      
        n1=2; n2=2; n3=2;
        global theta_v k_v sig_v mu_v;
        global theta_r k_r sig_r mu_r;
        global cor_Sv  q;
        
        x=sym('x'); v = sym('v'); r = sym('r');
        hx = sym('hx'); hv = sym('hv'); hr = sym('hr');    
        location_index = zeros(27,3);
        D = sym('D', [27 1]);
        
        % initialize location_index.
            for iii=1:3
                for jjj=1:3
                    for kkk=1:3
                      loc = locate_D(iii,jjj,kkk);
                    end
                end
            end   

         output = calculate_coefficients(i,j,k);
%        save('coeff27.mat','output');
      
%% sub functions        
    function result = calculate_coefficients(i_input,j_input,k_input)
        i=i_input; j=j_input; k=k_input;
        % discretization
        second_x=derivative_x2;        second_v=derivative_v2;        second_r=derivative_r2;
        first_x = derivative_x1;        first_v = derivative_v1;        first_r = derivative_r1;
        second_xv = derivative_xv;%       second_Sr = derivative_Sr;%       second_vr = derivative_vr;
        % matrix tryout
            matrix_K = (v/2) * (second_x) + (sig_v^2*v/2)*(second_v) + (sig_r^2*r/2)*(second_r)...
                     - (cor_Sv * sig_v *v) *(second_xv)...
                     + (k_r*(theta_r-r))*(first_r) - ((r-q)+v/2)*(first_x) + (k_v*(theta_v-v) + cor_Sv*sig_v*v) * (first_v);  
%          + (cor_Sr * sig_r *sqrt(r)*x) *(second_Sr) +  (cor_vr * sig_v*sig_r* sqrt(v*r)) *(second_vr)...      

        coeffs_in_K = sym('coeffs_in_K',[27 1]);           
        for idx=(1:27)
            coeffs_temp = coeffs(matrix_K,[D(idx)]);
            if size(coeffs_temp,2)~=1
                coeffs_in_K(idx) = coeffs_temp(2);
            end
        end
        % remove zeros:
        neibourlist = get_neighbours(i,j,k);
        coeffs_in_K = coeffs_in_K(neibourlist,:);
        list =cell(27,1);
        for i = 1:27
            list{i} = ['coeffs_in_K' num2str(i)];
        end
        result = subs(coeffs_in_K,list,{zeros(27,1)});
    end
    
    % this function returns all the neighbours of (i,j,k)
    function list=get_neighbours(i,j,k)
        list=[];
        for ii=max(1,i-1):min(3,i+1)
            for jj=max(1,j-1):min(3,j+1)
                for kk=max(1,k-1):min(3,k+1)
                    list=[list locate_D(ii,jj,kk)];
                end
            end
        end
    end

    % location mapping: 3D to 1D
    function location = locate_D(i,j,k)
     % mapping matric D(i,j,k) to D((n1+1)*(n2+1)*(N3+1))
     location = (i-1)*(3)*(3)...
              + (j-1)*(3)...
              + (k);
     location_index(location,:)=[i j k];     
    end


%% Derivatives discretization
    function  second_x=derivative_x2()
     if (i==3)||(i==1)
         if (i==1)
             % error = O(h);
             second_x =0;
%              second_x=(D(locate_D(i+2,j,k)) - 2*D(locate_D(i+1,j,k)) + D(locate_D(i,j,k)))/ (hx^2);
         else
             second_x=(D(locate_D(i-2,j,k)) - 2*D(locate_D(i-1,j,k)) + D(locate_D(i,j,k)))/ (hx^2);
%              second_x =0;
         end
     else 
        second_x = (D(locate_D(i+1,j,k)) - 2*D(locate_D(i,j,k)) + D(locate_D(i-1,j,k)))...
                 / (hx^2);
     end
    end

    function second_v=derivative_v2()
     if (j==1)||(j==3)
        second_v = 0;
     else
        second_v = (D(locate_D(i,j+1,k)) - 2*D(locate_D(i,j,k)) + D(locate_D(i,j-1,k)))...
                 / (hv^2); 
     end
    end

    function second_r=derivative_r2()
     if (k==1)||(k==3)
        second_r = 0;
     else
        second_r = (D(locate_D(i,j,k+1)) - 2*D(locate_D(i,j,k)) + D(locate_D(i,j,k-1)))...
                 / (hr^2);    
     end
    end

    function first_x = derivative_x1()
        if (i==1)
            first_x = ( D(locate_D(2,j,k)) - D(locate_D(1,j,k)) )/ hx;
%             first_x=0;
            return;
        end
        if (i==3)
%              first_x = ( D(locate_D(N1+1,j,k)) - D(locate_D(N1,j,k)) )...
%                      / hx ;           
            first_x=0;
        else
             first_x = ( D(locate_D(i+1,j,k)) - D(locate_D(i-1,j,k)) )...
                     / hx / 2;         
        end
    end

    function first_v = derivative_v1()
        if (j==1)
            first_v = ( D(locate_D(i,2,k)) - D(locate_D(i,1,k)) )...
                    / hv;
            return;
        end
        if (j==3)
             first_v = ( D(locate_D(i,n2+1,k)) - D(locate_D(i,n2,k)) )...
                     / hv ;           
        else
             first_v = ( D(locate_D(i,j+1,k)) - D(locate_D(i,j-1,k)) )...
                     / hv / 2;  
        end
    end

    function first_r = derivative_r1()
        if (k==1)
            first_r = ( D(locate_D(i,j,2)) - D(locate_D(i,j,1)) )...
                    / hr;
            return;
        end
        if (k==3)
             first_r = ( D(locate_D(i,j,n3+1)) - D(locate_D(i,j,n3)) )...
                     / hr ;           
        else
             first_r = ( D(locate_D(i,j,k+1)) - D(locate_D(i,j,k-1)) )...
                     / hr / 2;   
        end
    end

    function second_xv = derivative_xv()
        if (i==1 || j==1) || (i==3 || j==3)
            second_xv=0;
        else
            if cor_Sv>0
            second_xv = 1/2 *( ( D(locate_D(i+1,j+1,k)) - D(locate_D(i,j+1,k)) - ( D(locate_D(i+1,j,k))-D(locate_D(i,j,k)) ) )/hx/hv...
                      + (D(locate_D(i,j,k))-D(locate_D(i-1,j,k))-(D(locate_D(i,j-1,k))-D(locate_D(i-1,j-1,k))))/hx/hv );
            else
            second_xv = 1/2 *( ( D(locate_D(i,j+1,k)) - D(locate_D(i-1,j+1,k)) - ( D(locate_D(i,j,k))-D(locate_D(i-1,j,k)) ) )/hx/hv...
                      + (D(locate_D(i+1,j,k))-D(locate_D(i,j,k))-(D(locate_D(i+1,j-1,k))-D(locate_D(i,j-1,k))))/hx/hv );      
            end
        end
    end

%     function second_Sr = derivative_Sr()
%         if (i==1 || k==1) || (i==n1+1 || k==N3+1)
%             second_Sr=0;
%         else
%             if cor_Sr>0
%             second_Sr = 1/2 *( ( D(locate_D(i+1,j,k+1)) - D(locate_D(i,j,k+1)) - ( D(locate_D(i+1,j,k))-D(locate_D(i,j,k)) ) )/hx/hr...
%                       + (D(locate_D(i,j,k))-D(locate_D(i-1,j,k))-(D(locate_D(i,j,k-1))-D(locate_D(i-1,j,k-1))))/hx/hr );
%             else
%             second_Sr = 1/2 *( ( D(locate_D(i,j,k+1)) - D(locate_D(i-1,j,k+1)) - ( D(locate_D(i,j,k))-D(locate_D(i-1,j,k)) ) )/hx/hr...
%                       + (D(locate_D(i+1,j,k))-D(locate_D(i,j,k))-(D(locate_D(i+1,j,k-1))-D(locate_D(i,j,k-1))))/hx/hr );      
%             end
%         end
%     end
% 
%     function second_vr = derivative_vr()
%         if (j==1 || k==1) || (j==n2+1 || k==N3+1)
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
    end
 