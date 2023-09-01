function sectime=getimagesec(timestring)

%%
timehr=str2double(timestring(1:2));
timemin=str2double(timestring(4:5));
timesec=str2double(timestring(7:8));

%%
sectime=timehr*3600+timemin*60+timesec;