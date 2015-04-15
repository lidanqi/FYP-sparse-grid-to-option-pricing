function status=GO()
  [estimation_4, details_4,timespent_4] = MainFuncAC(4);
  [estimation_5, details_5,timespent_5] = MainFuncAC(5);
  [estimation_6, details_6,timespent_6] = MainFuncAC(6);
  [estimation_7, details_7,timespent_7] = MainFuncAC(7);
  [estimation_8, details_8,timespent_8] = MainFuncAC(8);
  summary=[estimation_4;estimation_5;estimation_6;estimation_7;estimation_8];
  time_summary=[timespent_4; timespent_5;timespent_6;timespent_7;timespent_8];
  a=clock;
  filename=[num2str(a(2)*1000000+a(3)*10000+a(4)*100+a(5)) '.mat'];
  % save everything in a file 
  save(filename);
  % play a sound to remind user to check result
  load handel
  sound(y,Fs)
  status=1;
end