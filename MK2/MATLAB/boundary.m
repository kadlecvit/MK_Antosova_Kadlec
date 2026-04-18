function[XB,YB] = boundary(R, uk, vk, s0, proj, ub, vb)
%Draw boundary lines (cutting edges)

%Transform to oblique aspect
% OPRAVA: uv_sd -> uvTosd (nazev vasi funkce)
[sb, db] = uvTosd(ub, vb, uk, vk);

%Project points
[XB, YB] = proj(R, sb, db, s0);
end