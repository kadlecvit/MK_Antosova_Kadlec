function [XC,YC] = drawContinent(file, R, uk, vk, u0, proj)
%Draw continent in selected projection
points = load(file);

%Calculate lat and lon
u = points(:,1) * pi / 180;
v = points(:,2) * pi / 180;

%Convert to oblique aspect
[s,d] = uvTosd(u, v, uk, vk);

%Remove points near equator (gnomonic projection diverges at s=0)
idx = find(s > 5*pi/180);
s = s(idx);
d = d(idx);

%Draw in selected projection 
[XC, YC] = proj(R, s, d, u0);

end