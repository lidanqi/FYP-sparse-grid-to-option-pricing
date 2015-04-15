% sparse grid combination step
% Compound Option
function [estimation, details,timespent] = MainFuncAC(requiredlevel)
% set benchmark for comparison: MC + POSR
S = [0.8, 0.9, 1, 1.1, 1.2];
benchmark_A = [19.9987, 10.9820, 5.4899, 2.6295, 1.2388];
benchmark_AC = [0.1072, 0.6119, 1.5618, 2.5233, 3.1928];

mapping=zeros(200,3);
count=0;
StartEnd=zeros(3,2);
for level=requiredlevel:requiredlevel+2
    StartEnd(level-requiredlevel+1,1)=count+1;
    for i=0:level
        for j=0:(level-i)
            k=level-i-j;
            count=count+1;
            mapping(count,:)=[i j k];
        end
    end
    StartEnd(level-requiredlevel+1,2)=count;
end

mapping = mapping(1:count,:);
% start 
timer = clock;
sums=zeros(3,10);
details=zeros(count,13);
answers=zeros(count,10);

for level=requiredlevel:requiredlevel+2

level_s = level - requiredlevel + 1;
% list=zeros(nchoosek(level+dimension-1,dimension-1),4);

parfor idx=StartEnd(level_s,1):StartEnd(level_s,2)
        % number of points on this grid of combination (i,j,k)
        i = mapping(idx,1);
        j = mapping(idx,2);
        k = mapping(idx,3);

        fprintf('levels: %d %d %d \n', i,j,k);
        % generate a file containing estimation 
        filename=[num2str(i) num2str(j) num2str(k) '_new.mat'];
        IScomputed = exist(filename);
        if (IScomputed==0)
            [~,~,est,~] = mainAC(S,128,'level',[i j k]);
            ts=500;
            while (est(:,2)==0)
               [~,~,est,~] = mainAC(S,ts,'level',[i j k],'w',1.8);
               ts=ts+500;
            end
            mysave(filename,est);
        else
            loadedEst = load(filename);
            est=loadedEst.est;
        end
        temp = est(:,1);
        tempC = est(:,2);
    %  temp = 0;
        details(idx,:)=[i j k temp' tempC'];
        if isnan(temp)
            temp=0;
        end
        answers(idx,:) = [temp' tempC'];
    %   list(index,:)=[i j k points];
end

sums(level_s,:)=sum(answers(StartEnd(level_s,1):StartEnd(level_s,2),:));

end

sums';
estimation = (sums(1,:) - 2*sums(2,:) + sums(3,:))*100; 
estimation_A = estimation(1:5);
estimation_AC = estimation(6:10);

fprintf('=====================================================================');
fprintf('\nEstid Daughter Option price: '); fprintf('%6g ',estimation_A); 
fprintf('\n               Benchmark   : '); fprintf('%6g ',benchmark_A);
fprintf('\n---------------------------------------------------------------------');
fprintf('\nEstid Mother Option price: '); fprintf('%6g ',estimation_AC); 
fprintf('\n               Benchmark   : '); fprintf('%6g ',benchmark_AC);
timespent = etime(clock,timer);
fprintf('\nTotal time spent: %4d s \n',timespent);
end



function status = mysave(filename, est)
    save(filename,'est');
end