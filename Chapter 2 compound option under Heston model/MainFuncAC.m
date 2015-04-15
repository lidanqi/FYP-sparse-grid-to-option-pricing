% sparse grid combination step
% Compound Option
function [estimation, details,timespent] = MainFuncAC(requiredlevel,varargin)
 requiredlevel=requiredlevel+1;
 p = inputParser;
 % THIS TYPE DETERMINES EURO VS AME
 % EURO: 1;   AMERICAN:2
 p.addParamValue('type',1); 
 p.parse(varargin{:});
 
 OptionType=p.Results.type;
 
% set benchmark for comparison: MC + POSR
% alter by commenting
S = [0.8, 0.9, 1, 1.1, 1.2];
if (OptionType==1)
%     benchmark_AC = [0.0685, 0.6450, 2.7880, 7.2315, 13.7291]; % nega cor
    benchmark_AC = [0.3088, 1.1427, 3.2084, 7.1530, 13.1855]; % posi cor
    filename_last = '_p_E.mat';
else
%     benchmark_AC = [0.0785, 0.6898, 2.9567, 7.6920, 14.7089];
    benchmark_AC = [0.3128,1.1940, 3.4061, 7.6646, 14.2716];
    filename_last = '_p_A.mat';
end

mapping=zeros(200,2);
count=0;
StartEnd=zeros(2,2);
for level=requiredlevel-1:requiredlevel
    StartEnd(level-requiredlevel+2,1)=count+1;
    for i=0:level
             j=level-i;
            count=count+1;
            mapping(count,:)=[i j];
    end
    StartEnd(level-requiredlevel+2,2)=count;
end

mapping = mapping(1:count,:);
% start 
timer = clock;
sums=zeros(2,10);
details=zeros(count,13);
answers=zeros(count,10);
maxSubTime=0;

for level=requiredlevel-1:requiredlevel

level_s = level - requiredlevel + 2;
% list=zeros(nchoosek(level+dimension-1,dimension-1),4);

parfor idx=StartEnd(level_s,1):StartEnd(level_s,2)
        % number of points on this grid of combination (i,j,k)
        i = mapping(idx,1);
        j = mapping(idx,2);

        fprintf('levels: %d %d \n', i,j);
        % generate a file containing estimation 
        
        filename=[num2str(i) num2str(j) filename_last];
        IScomputed = exist(filename);
        if (IScomputed==0)
            [~,~,est,~,subTime] = mainAC(S,200,'level',[i j],'type',OptionType);
            ts=500;
            while ( (est(:,1)==0) )
                [~,~,est,~,subTime] = mainAC(S,ts,'level',[i j],'type',OptionType,'w',1.8);
                ts=ts+1000;
            end
            mysave(filename,est);
            maxSubTime=max(subTime,maxSubTime);
        else
            loadedEst = load(filename);
            est=loadedEst.est;
        end
        temp = est(:,1);
        tempC = est(:,2);
    %  temp = 0;
        details(idx,:)=[i j 0 temp' tempC'];
        if isnan(temp)
            temp=0;
        end
        answers(idx,:) = [temp' tempC'];
    %   list(index,:)=[i j k points];
end

sums(level_s,:)=sum(answers(StartEnd(level_s,1):StartEnd(level_s,2),:));

end

estimation = (sums(2,:) - sums(1,:))*100; 
estimation_A = estimation(1:5);
estimation_AC = estimation(6:10);

fprintf('=====================================================================');
fprintf('\nEstid Mother Option price: '); fprintf('%8.4f ',estimation_AC); 
fprintf('\n               Benchmark : '); fprintf('%8.4f ',benchmark_AC);
timespent = etime(clock,timer);
fprintf('\nMax sub time spent: %4d s ',maxSubTime);
fprintf('\nTotal time spent  : %4d s \n',timespent);
end



function status = mysave(filename, est)
    save(filename,'est');
end