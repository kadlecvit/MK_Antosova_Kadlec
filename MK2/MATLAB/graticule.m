function [XM, YM, XP, YP] = graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, R, uk, vk, u0, proj)
%Create graticule: list of meridians and parallels

%Create list of meridians
XM = []; YM = [];
for v = vmin:Dv:vmax

    %Create meridian
    um = umin:du:umax;
    n = length(um);
    vm = ones(1, n)*v;
    
    %Convert to oblique aspect
    [sm, dm] = uvTosd(um, vm, uk, vk);

    %Compute xm, ym
    [xm, ym] = proj(R, sm, dm, u0);

    %Add meridian to the list
    XM = [XM; xm];
    YM = [YM; ym];

end

%Create list of parallels
XP = []; YP = [];
for u = umin:Du:umax
    
    %Create parallel
    vp = vmin:dv:vmax;
    n = length(vp);
    up = ones(1, n)*u;
    
    %Convert to oblique aspect
    [sp, dp] = uvTosd(up, vp, uk, vk);

    %Compute xp, yp
    [xp, yp] = proj(R, sp, dp, u0);

    %Add parallel to the list
    XP = [XP; xp];
    YP = [YP; yp];

end

end