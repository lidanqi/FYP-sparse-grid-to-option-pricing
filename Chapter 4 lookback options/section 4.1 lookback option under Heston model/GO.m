function status=GO()
%   MC_benchmark=[5.360888540966589,5.002942737696325,4.764786807062194,4.631523508332920,4.588000169844665];
  
  [estimation_4, details_4,timespent_4] = MainFuncLB(4);
  [estimation_5, details_5,timespent_5] = MainFuncLB(5);
  [estimation_6, details_6,timespent_6] = MainFuncLB(6);
  [estimation_7, details_7,timespent_7] = MainFuncLB(7);
%   [estimation_8, details_8,timespent_8] = MainFuncLB(8);
  summary=[estimation_4;...
         estimation_5;estimation_6;estimation_7];%estimation_8];
  tt=[timespent_4,...
      timespent_5,timespent_6,timespent_7];%timespent_8];    
  tt=tt';
  a=clock;
  filename=[num2str(a(2)*1000000+a(3)*10000+a(4)*100+a(5)) '_order3.mat'];
  save(filename);
  load handel
  sound(y,Fs)
  status=1;
end