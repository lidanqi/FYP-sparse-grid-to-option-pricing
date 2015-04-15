function status=GO()
  
  [estimation_1, details_1,timespent_1] = MainFuncC(1);
  [estimation_2, details_2,timespent_2] = MainFuncC(2);
  [estimation_3, details_3,timespent_3] = MainFuncC(3);
  [estimation_4, details_4,timespent_4] = MainFuncC(4);
  [estimation_5, details_5,timespent_5] = MainFuncC(5);
  a=clock;
  filename=[num2str(a(2)*1000000+a(3)*10000+a(4)*100+a(5)) '.mat'];
  % save results
  save(filename);
  % play sound to remind 
  load handel
  sound(y,Fs)
  status=1;
end