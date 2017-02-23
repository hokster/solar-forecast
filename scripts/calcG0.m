function G0 = calcG0(year,jday,dt,lat,lon,lonstd)
Gsc = 1361; % Solar constant, W/m^2
% daylight savings time arrays
dstjsp = [2012 71; ...
          2013 69; ...
          2014 68; ...
          2015 67];
dstjfa = [2012 309; ...
          2013 307; ...
          2014 306; ...
          2015 305];
% correct decimal time to EST/EDT instead of UTC
dt(jday >= dstjsp(dstjsp(:,1) == year,2) & jday < dstjsp(dstjfa(:,1) == year,2)) = ...
    dt(jday >= dstjsp(dstjsp(:,1) == year,2) & jday < dstjsp(dstjfa(:,1) == year,2)) - 4;
dt(jday < dstjsp(dstjsp(:,1) == year,2) | jday >= dstjsp(dstjfa(:,1) == year,2)) = ...
    dt(jday < dstjsp(dstjsp(:,1) == year,2) | jday >= dstjsp(dstjfa(:,1) == year,2)) - 5;

% calculate if it's a leap year
leap = (mod(year,4) == 0);
% calculate hour angle
if leap == 1
    B = (jday - 1)*360/366;
else
    B = (jday - 1)*360/365;
end
% analemma correction
Et = 229.2*(0.000075+0.001868*cosd(B)-0.032077*sind(B)-0.014615*cosd(2*B)-0.04089*sind(2*B));
tl = 4*(lon-lonstd);
TC = tl + Et;
% if convert decimal time to minutes and add time correction
tsol = TC + 60*dt;
omega = 15*(tsol/60-12); % hour angle
if leap == 1
    delta = 23.45*sind(360*(284+jday)/366); % declination angle
    G0 = Gsc*(1+0.0033*cosd(360*jday/366)).*(sind(lat)*sind(delta) + cosd(lat)*(cosd(delta).*cosd(omega)));
else
    delta = 23.45*sind(360*(284+jday)/365); % declination angle
    G0 = Gsc*(1+0.0033*cosd(360*jday/365)).*(sind(lat)*sind(delta) + cosd(lat)*(cosd(delta).*cosd(omega)));
end