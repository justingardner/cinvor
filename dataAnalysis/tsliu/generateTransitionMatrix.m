function transition_total = generateTransitionMatrix(retDir,cinDir,prf_scan_num)

origDirectory = pwd;
cd(retDir);
v_ret = newView;
v_ret = viewSet(v_ret,'curGroup','Averages');
transition_ret = viewGet(v_ret,'scan2mag',prf_scan_num,1);
deleteView(v_ret);
mrQuit;
cd(cinDir);
v_cinvor = newView;
v_cinvor = viewSet(v_cinvor,'curGroup','Raw');
transition_cinvor = viewGet(v_cinvor,'scan2mag',1,1);
deleteView(v_cinvor);
mrQuit;
cd(origDirectory);

transition_total = inv(transition_ret)*transition_cinvor;