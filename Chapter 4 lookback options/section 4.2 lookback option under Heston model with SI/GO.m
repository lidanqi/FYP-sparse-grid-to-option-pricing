function status=GO()
  MC_benchmark=[43.394729035906046;34.397140548268084;26.488850648389604;20.747655526437434;18.509054076545922];
  
  [estimation_1, details_1,timespent_1] = MainFuncLB(1);
  [estimation_2, details_2,timespent_2] = MainFuncLB(2);
  [estimation_3, details_3,timespent_3] = MainFuncLB(3);
  [estimation_4, details_4,timespent_4] = MainFuncLB(4);
  [estimation_5, details_5,timespent_5] = MainFuncLB(5);
  summary=[estimation_1;estimation_2;estimation_3;estimation_4;estimation_5];
  a=clock;
  filename=[num2str(a(2)*1000000+a(3)*10000+a(4)*100+a(5)) '.mat'];
  save(filename);
  load handel
  sound(y,Fs)
  status=1;
end