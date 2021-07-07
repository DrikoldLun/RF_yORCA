addpath('Tinycodes');
Lineloca = [-133.55 -8; %A
            -132.45 -4]; %B
[sta, stlo, stla]=textread('latlon_all.dat','%s %f %f %*f %*f %*f');
[dis_along,projlat,projlon] = staproject(stla,stlo,Lineloca);
fid = fopen('sta_project_all.dat', 'w+');
for i = 1:length(dis_along)
    fprintf(fid, '%s %f %f %f\n',sta{i},dis_along(i),projlon(i),projlat(i));
end
fclose(fid);

function [dis_along,projlat,projlon] = staproject(stla,stlo,Lineloca)
[dis, azi] = distance(Lineloca(1,2),Lineloca(1,1),Lineloca(2,2),Lineloca(2,1));
[dis_center, azi_center] = distance(Lineloca(1,2),Lineloca(1,1),stla,stlo);
dis_along = deg2km(atand(tand(dis_center)).*cosd(azi - azi_center)); %km
[projlat,projlon] =latlon_from(Lineloca(1,2),Lineloca(1,1),azi,dis_along);
end