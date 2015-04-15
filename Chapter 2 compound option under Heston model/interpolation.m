% this function is used for multilinear interpolation
function estimation=interpolation(Svec,vvec,resultVect,S,v)
sizeS = size(Svec,1);
sizev = size(vvec,1);
[Y,X] =meshgrid(vvec,Svec);

% reshape result vector into 3-D
[res] = mymesh(resultVect);

% estimation = interp3( X,Y,Z,resultVect ,S,v,r);
 F = griddedInterpolant(X,Y,res,'cubic');
% F = griddedInterpolant(X,Y,Z,res,'linear');
numS = max(size(S,1),size(S,2));
estimation = zeros(numS,1);
for i=1:numS
 estimation(i) = F(S(i),v);
end
% estimation = interp3(x,y,z, res, S,v,r); 


function [resultVect] = mymesh(V)
    N1 = sizeS-1;
    N2 = sizev-1;
    
    resultVect=zeros(sizeS,sizev);
%     x=resultVec;y=x; z=x; 
        for ii=1:N1+1
            for jj=1:N2+1
%                     x(ii,jj,kk)=s(ii);
%                     y(ii,jj,kk)=v(jj);
%                     z(ii,jj,kk)=r(kk);
                    resultVect(ii,jj) = V(locate_D(ii,jj));
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

end

