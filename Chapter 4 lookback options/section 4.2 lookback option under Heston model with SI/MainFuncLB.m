% sparse grid combination step
% Compound Option
function [estimation, details,timespent] = MainFuncLB(requiredlevel)
% set benchmark for comparison: MC + POSR

benchmark = 13.2449;


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
sums=zeros(3,5);
details=zeros(count,8);
answers=zeros(count,5);

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
        filename=[ num2str(i) num2str(j) num2str(k) '.mat'];
        path_filename=['data\' filename];
        IScomputed = exist(path_filename);
        if (IScomputed==0)
            [~,est,~] = mainLB(100,'level',[i j k]);
            ts=500;
            while (est==0)
                [~,est,~] = mainLB(ts,'level',[i j k],'w',1.8);
                ts=ts+500;
            end
            mysave(filename,est);
        else
            loadedEst = load(path_filename);
            est=loadedEst.est;
        end
        temp = est';
    %  temp = 0;
        details(idx,:)=[i j k temp];
        if isnan(temp)
            temp=0;
        end
        answers(idx,:) = [temp];
    %   list(index,:)=[i j k points];
end

sums(level_s,:)=sum(answers(StartEnd(level_s,1):StartEnd(level_s,2),:));

end

sums';
estimation = (sums(1,:) - 2*sums(2,:) + sums(3,:)); 
estimation_A = estimation;

fprintf('=====================================================================');
fprintf('\nEstid Daughter Option price: '); fprintf('%6g ',estimation_A); 
fprintf('\n               Benchmark   : '); fprintf('%6g ',benchmark);
timespent = etime(clock,timer);
fprintf('\nTotal time spent: %6g s \n',timespent);
end



function status = mysave(filename, est)
    save(['data\' filename],'est');
    status=1;
end