% Code used to plot Figure 17
tic
clear all
% Machine's characteristics
tp   = 60e-3;        % pole pitch (m)
tm   = 55e-3;        % PM length in x direction (m)
hm   = 10e-3;        % PM height in y direction (m)
e    = 1e-3;         % Air-gap thickness (m)
hst  = 30e-3;        % Stator total height (m)
hs   = 20e-3;        % Slot height (m)
hmbi = 10e-3;        % Moving armature height (moving back iron height)
ws   = 10e-3;        % Slot opening (m)
ts   = 2*ws;         % Slot pitch (m)
la   = 1;            % Active length (m)
Br   = 1.2;          % PM remanent induction (residual induction) (T)
mu0  = 4*pi*1e-7;    % Permeability of vacuum (H/m)
mur1 = 1;            % Relative permeability of air
mur2 = 7500;         % Relative permeability of statot iron
mur3 = 7500;         % Relative permeability of moving armature iron

% Units
% m   : meter
% T   : Tesla
% H/m : Henry per meter

% Definition of the reluctance elements
h  = 0.5e-3;    % Height of the reluctance element (m)
l  = 0.5e-3;    % Length of the reluctance element (m)

% Reluctances and permeances of reluctace elements located in air and slots regions
Rv1 = (1/mur1/mu0)*l/(h*la);
Rw1 = (1/mur1/mu0)*h/(l*la);
Pv1 = mur1*mu0*(h*la/l);
Pw1 = mur1*mu0*(l*la/h);

% Reluctances and permeances of reluctace elements located in stator iron core
Rv20 = (1/mur2/mu0)*l/(h*la);
Rw20 = (1/mur2/mu0)*h/(l*la);
Pv20 = mur2*mu0*(h*la/l);
Pw20 = mur2*mu0*(l*la/h);

% Reluctances and permeances of reluctace elements located in movinr armature iron core
Rv21 = (1/mur3/mu0)*l/(h*la);
Rw21 = (1/mur3/mu0)*h/(l*la);
Pv21 = mur3*mu0*(h*la/l);
Pw21 = mur3*mu0*(l*la/h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Constitution of matrix [P] (equation (8) of the contribution) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of elements in the stator armature
m0s = round((ts-ws)/2/l);  % Number of elements in half a tooth in x direction
p0s = round((hst-hs)/h);   % Number of elements in the stator back iron in y direction
m  = 12*m0s;               % Total number of elements of the stator in x direction
p  = 3*p0s;                % Total number of element of stator in y direction

% Number of elements in the moving armature (the air-gap is supposed to be part of the moving armature)
m0m = round((tp-tm)/2/l);  % Number of elements in half the air-gap between two adjacent PM in x direction
p0m = round(e/h);          % Number of elements in the air-gap in y direction
p0  = round(hmbi/h);       % Number of elements in the moving armature iron in y direction
p1  = round((hm+e)/h);     % Number of elements in the magnetic air-gap (hm + e) in y direction

nn = (p+p0+p1)*m;          % Total number of elements or nodes

% Mobile armature iron 
% First line
r01 = 1;
s01 = 1;
nnum(r01,s01)=(r01-1)*m+s01;
MP(nnum(r01,s01),nnum(r01,s01)+m-1) = -1/Rv21;
MP(nnum(r01,s01),nnum(r01,s01)+1)   =  1/Rv21;
MP(nnum(r01,s01),nnum(r01,s01)+m)   =  1/Rw21;
M(nnum(r01,s01),nnum(r01,s01))      = -MP(nnum(r01,s01),nnum(r01,s01)+m-1)+MP(nnum(r01,s01),nnum(r01,s01)+1)+MP(nnum(r01,s01),nnum(r01,s01)+m);
M(nnum(r01,s01),nnum(r01,s01)+m-1)  = -MP(nnum(r01,s01),nnum(r01,s01)+m-1);
M(nnum(r01,s01),nnum(r01,s01)+1)    = -MP(nnum(r01,s01),nnum(r01,s01)+1);
M(nnum(r01,s01),nnum(r01,s01)+m)    = -MP(nnum(r01,s01),nnum(r01,s01)+m);
for s01 = 2:m-1
    nnum(r01,s01)=(r01-1)*m+s01;
    MP(nnum(r01,s01),nnum(r01,s01)-1) =  1/Rv21;
    MP(nnum(r01,s01),nnum(r01,s01)+1) =  1/Rv21;
    MP(nnum(r01,s01),nnum(r01,s01)+m) =  1/Rw21;
    M(nnum(r01,s01),nnum(r01,s01))    =  MP(nnum(r01,s01),nnum(r01,s01)-1)+MP(nnum(r01,s01),nnum(r01,s01)+1)+MP(nnum(r01,s01),nnum(r01,s01)+m);
    M(nnum(r01,s01),nnum(r01,s01)-1)  = -MP(nnum(r01,s01),nnum(r01,s01)-1);
    M(nnum(r01,s01),nnum(r01,s01)+1)  = -MP(nnum(r01,s01),nnum(r01,s01)+1);
    M(nnum(r01,s01),nnum(r01,s01)+m)  = -MP(nnum(r01,s01),nnum(r01,s01)+m);
end
s01 = m;
nnum(r01,s01)=(r01-1)*m+s01;
MP(nnum(r01,s01),nnum(r01,s01)-1)   =  1/Rv21;
MP(nnum(r01,s01),nnum(r01,s01)+m)   =  1/Rw21;
MP(nnum(r01,s01),nnum(r01,s01)-m+1) = -1/Rv21;
M(nnum(r01,s01),nnum(r01,s01))      =  MP(nnum(r01,s01),nnum(r01,s01)-1)+MP(nnum(r01,s01),nnum(r01,s01)+m)-MP(nnum(r01,s01),nnum(r01,s01)-m+1);
M(nnum(r01,s01),nnum(r01,s01)-1)    = -MP(nnum(r01,s01),nnum(r01,s01)-1);
M(nnum(r01,s01),nnum(r01,s01)+m)    = -MP(nnum(r01,s01),nnum(r01,s01)+m);
M(nnum(r01,s01),nnum(r01,s01)-m+1)  = -MP(nnum(r01,s01),nnum(r01,s01)-m+1);

% First column
s01 = 1;
for r01 = 2:p0-1
    nnum(r01,s01)=(r01-1)*m+s01;
    MP(nnum(r01,s01),nnum(r01,s01)+m-1) = -1/Rv21;
    MP(nnum(r01,s01),nnum(r01,s01)+1)   =  1/Rv21;
    MP(nnum(r01,s01),nnum(r01,s01)+m)   =  1/Rw21;
    MP(nnum(r01,s01),nnum(r01,s01)-m)   =  1/Rw21;
    M(nnum(r01,s01),nnum(r01,s01))      = -MP(nnum(r01,s01),nnum(r01,s01)+m-1)+MP(nnum(r01,s01),nnum(r01,s01)+1)+MP(nnum(r01,s01),nnum(r01,s01)+m)+MP(nnum(r01,s01),nnum(r01,s01)-m);
    M(nnum(r01,s01),nnum(r01,s01)+m-1)  = -MP(nnum(r01,s01),nnum(r01,s01)+m-1);
    M(nnum(r01,s01),nnum(r01,s01)+1)    = -MP(nnum(r01,s01),nnum(r01,s01)+1);
    M(nnum(r01,s01),nnum(r01,s01)+m)    = -MP(nnum(r01,s01),nnum(r01,s01)+m);
    M(nnum(r01,s01),nnum(r01,s01)-m)    = -MP(nnum(r01,s01),nnum(r01,s01)-m);
end

% Last column
s01 = m;
for r01 = 2:p0-1
    nnum(r01,s01)=(r01-1)*m+s01;
    MP(nnum(r01,s01),nnum(r01,s01)-m+1) = -1/Rv21;
    MP(nnum(r01,s01),nnum(r01,s01)-1)   =  1/Rv21;
    MP(nnum(r01,s01),nnum(r01,s01)+m)   =  1/Rw21;
    MP(nnum(r01,s01),nnum(r01,s01)-m)   =  1/Rw21;
    M(nnum(r01,s01),nnum(r01,s01))      = -MP(nnum(r01,s01),nnum(r01,s01)-m+1)+MP(nnum(r01,s01),nnum(r01,s01)-1)+MP(nnum(r01,s01),nnum(r01,s01)+m)+MP(nnum(r01,s01),nnum(r01,s01)-m);
    M(nnum(r01,s01),nnum(r01,s01)-m+1)  = -MP(nnum(r01,s01),nnum(r01,s01)-m+1);
    M(nnum(r01,s01),nnum(r01,s01)-1)    = -MP(nnum(r01,s01),nnum(r01,s01)-1);
    M(nnum(r01,s01),nnum(r01,s01)+m)    = -MP(nnum(r01,s01),nnum(r01,s01)+m);
    M(nnum(r01,s01),nnum(r01,s01)-m)    = -MP(nnum(r01,s01),nnum(r01,s01)-m);
end

% The rest of elements
for r01 = 2:p0-1
    for s01 = 2:m-1
        nnum(r01,s01)=(r01-1)*m+s01;
        MP(nnum(r01,s01),nnum(r01,s01)-1) =  1/Rv21;
        MP(nnum(r01,s01),nnum(r01,s01)+1) =  1/Rv21;
        MP(nnum(r01,s01),nnum(r01,s01)+m) =  1/Rw21;
        MP(nnum(r01,s01),nnum(r01,s01)-m) =  1/Rw21;
        M(nnum(r01,s01),nnum(r01,s01))    =  MP(nnum(r01,s01),nnum(r01,s01)-1)+MP(nnum(r01,s01),nnum(r01,s01)+1)+MP(nnum(r01,s01),nnum(r01,s01)+m)+MP(nnum(r01,s01),nnum(r01,s01)-m);
        M(nnum(r01,s01),nnum(r01,s01)-1)  = -MP(nnum(r01,s01),nnum(r01,s01)-1);
        M(nnum(r01,s01),nnum(r01,s01)+1)  = -MP(nnum(r01,s01),nnum(r01,s01)+1);
        M(nnum(r01,s01),nnum(r01,s01)+m)  = -MP(nnum(r01,s01),nnum(r01,s01)+m);
        M(nnum(r01,s01),nnum(r01,s01)-m)  = -MP(nnum(r01,s01),nnum(r01,s01)-m);
    end
end


% Mobile iron / Magnetic air-gap boundary
% Last line of mobile iron having a link with the magnetic air-gap
r01 = p0;
s01 = 1;
nnum(r01,s01)=(r01-1)*m+s01;
MP(nnum(r01,s01),nnum(r01,s01)-m)   =  1/Rw21;
MP(nnum(r01,s01),nnum(r01,s01)+m-1) = -1/Rv21;
MP(nnum(r01,s01),nnum(r01,s01)+1)   =  1/Rv21;
MP(nnum(r01,s01),nnum(r01,s01)+m)   =  1/(Rw21/2+Rw1/2);
M(nnum(r01,s01),nnum(r01,s01))      =  MP(nnum(r01,s01),nnum(r01,s01)+1)-MP(nnum(r01,s01),nnum(r01,s01)+m-1)+MP(nnum(r01,s01),nnum(r01,s01)+m)+MP(nnum(r01,s01),nnum(r01,s01)-m);
M(nnum(r01,s01),nnum(r01,s01)-m)    = -MP(nnum(r01,s01),nnum(r01,s01)-m);
M(nnum(r01,s01),nnum(r01,s01)+m-1)  = -MP(nnum(r01,s01),nnum(r01,s01)+m-1);
M(nnum(r01,s01),nnum(r01,s01)+1)    = -MP(nnum(r01,s01),nnum(r01,s01)+1);
M(nnum(r01,s01),nnum(r01,s01)+m)    = -MP(nnum(r01,s01),nnum(r01,s01)+m);
for s01 = 2:m-1
    nnum(r01,s01)=(r01-1)*m+s01;
    MP(nnum(r01,s01),nnum(r01,s01)-m) =  1/Rw21;
    MP(nnum(r01,s01),nnum(r01,s01)-1) =  1/Rv21;
    MP(nnum(r01,s01),nnum(r01,s01)+1) =  1/Rv21;
    MP(nnum(r01,s01),nnum(r01,s01)+m) =  1/(Rw21/2+Rw1/2);
    M(nnum(r01,s01),nnum(r01,s01))    =  MP(nnum(r01,s01),nnum(r01,s01)+1)+MP(nnum(r01,s01),nnum(r01,s01)-1)+MP(nnum(r01,s01),nnum(r01,s01)+m)+MP(nnum(r01,s01),nnum(r01,s01)-m);
    M(nnum(r01,s01),nnum(r01,s01)-1)  = -MP(nnum(r01,s01),nnum(r01,s01)-1);
    M(nnum(r01,s01),nnum(r01,s01)+1)  = -MP(nnum(r01,s01),nnum(r01,s01)+1);
    M(nnum(r01,s01),nnum(r01,s01)+m)  = -MP(nnum(r01,s01),nnum(r01,s01)+m);
    M(nnum(r01,s01),nnum(r01,s01)-m)  = -MP(nnum(r01,s01),nnum(r01,s01)-m);
end
s01 = m;
nnum(r01,s01)=(r01-1)*m+s01;
MP(nnum(r01,s01),nnum(r01,s01)-m)   =  1/Rw21;
MP(nnum(r01,s01),nnum(r01,s01)-1)   =  1/Rv21;
MP(nnum(r01,s01),nnum(r01,s01)+m)   =  1/(Rw21/2+Rw1/2);
MP(nnum(r01,s01),nnum(r01,s01)-m+1) = -1/Rv21;
M(nnum(r01,s01),nnum(r01,s01))      =  MP(nnum(r01,s01),nnum(r01,s01)-1)-MP(nnum(r01,s01),nnum(r01,s01)-m+1)+MP(nnum(r01,s01),nnum(r01,s01)+m)+MP(nnum(r01,s01),nnum(r01,s01)-m);
M(nnum(r01,s01),nnum(r01,s01)-1)    = -MP(nnum(r01,s01),nnum(r01,s01)-1);
M(nnum(r01,s01),nnum(r01,s01)-m+1)  = -MP(nnum(r01,s01),nnum(r01,s01)-m+1);
M(nnum(r01,s01),nnum(r01,s01)+m)    = -MP(nnum(r01,s01),nnum(r01,s01)+m);
M(nnum(r01,s01),nnum(r01,s01)-m)    = -MP(nnum(r01,s01),nnum(r01,s01)-m);

% Magnetic air-gap
% First line of the Magnetic air-gap having a link with the mobile iron
r11 = p0+1;
s11 = 1;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw21/2+Rw1/2);
MP(nnum(r11,s11),nnum(r11,s11)+m-1) = -1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw1;
M(nnum(r11,s11),nnum(r11,s11))      =  MP(nnum(r11,s11),nnum(r11,s11)-m)-MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
for s11 = 2:m-1
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m) =  1/(Rw21/2+Rw1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-1) =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) =  1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) =  1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    =  MP(nnum(r11,s11),nnum(r11,s11)-m)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
s11 = m;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw21/2+Rw1/2);
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw1;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) = -1/Rv1;
M(nnum(r11,s11),nnum(r11,s11))      =  MP(nnum(r11,s11),nnum(r11,s11)-m)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)-MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% First colum
s11 = 1;
for r11 = p0+2:p0+p1-1
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+m-1) =  -1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =   1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =   1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =   1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = -MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% Last colum
s11 = m;
for r11 = p0+2:p0+p1-1
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m+1) =  -1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =   1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =   1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =   1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))      = -MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% The rest of elements in magnetic air-gap region
for r11 = p0+2:p0+p1-1
    for s11 = 2:m-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

% Magnetic air-gap / Stator boundary
% Last line of magnetic air-gap havin a ling with the first line of the stator armature
r11 = p0+p1;
s11 = 1;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) = -1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw1/2+Rw20/2);
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)-MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m)    =  -MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  =  -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    =  -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    =  -MP(nnum(r11,s11),nnum(r11,s11)+m);
for s11 = 2:m0s
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
for dent = 1:2
    for s11 = 3*m0s+1+(dent-1)*4*m0s:5*m0s+(dent-1)*4*m0s
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end
for enc = 1:3
    for s11 = m0s+1+(enc-1)*4*m0s:3*m0s+(enc-1)*4*m0s
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end
for s11 = 11*m0s+1:m-1
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
s11 = m;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw1;
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv1;
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/(Rw1/2+Rw20/2);
MP(nnum(r11,s11),nnum(r11,s11)-m+1) = -1/Rv1;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)-MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% Stator armature
% First column
s11 = 1;
r11 = p0+p1+1;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) = -1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw20/2+Rw1/2);
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)-MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
for r11 = p0+p1+2:p0+p1+p-1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+m-1) = -1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)-MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
r11 = p0+p1+p;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)+m-1) = -1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)+1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)+1)-MP(nnum(r11,s11),nnum(r11,s11)+m-1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)+m-1)  = -MP(nnum(r11,s11),nnum(r11,s11)+m-1);
M(nnum(r11,s11),nnum(r11,s11)+1)    = -MP(nnum(r11,s11),nnum(r11,s11)+1);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% Last column
s11 = m;
r11 = p0+p1+1;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) = -1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/(Rw20/2+Rw1/2);
MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)-MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
for r11 = p0+p1+2:p0+p1+p-1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-m+1) = -1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
    MP(nnum(r11,s11),nnum(r11,s11)+m)   =  1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)-MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
    M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+m)    = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end
r11 = p0+p1+p;
nnum(r11,s11)=(r11-1)*m+s11;
MP(nnum(r11,s11),nnum(r11,s11)-m+1) = -1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-1)   =  1/Rv20;
MP(nnum(r11,s11),nnum(r11,s11)-m)   =  1/Rw20;
M(nnum(r11,s11),nnum(r11,s11))      = MP(nnum(r11,s11),nnum(r11,s11)-1)-MP(nnum(r11,s11),nnum(r11,s11)-m+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
M(nnum(r11,s11),nnum(r11,s11)-m+1)  = -MP(nnum(r11,s11),nnum(r11,s11)-m+1);
M(nnum(r11,s11),nnum(r11,s11)-1)    = -MP(nnum(r11,s11),nnum(r11,s11)-1);
M(nnum(r11,s11),nnum(r11,s11)-m)    = -MP(nnum(r11,s11),nnum(r11,s11)-m);

% Columns s11 = m0s, 5*m0s, 9*m0s
for dent = 1:3
    s11 = m0s+(dent-1)*4*m0s;
    r11 = p0+p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rw1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p0+p1+2:p0+p1+2*p0s
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    for r11 = p0+p1+2*p0s+1:p0+p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p0+p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% Columns s11 = 3*m0s+1, 7*m0s+1, 11*m0s+1
for dent = 1:3
    s11 = 3*m0s+1+(dent-1)*4*m0s;
    r11 = p0+p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rw1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p0+p1+2:p0+p1+2*p0s
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    for r11 = p0+p1+2*p0s+1:p0+p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p0+p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% First and last half teeth excluding s11 = 1 and 1+11*m0s and s11 = m0s and 12*m0s
for dent = 1:2
    for s11 = 2+(dent-1)*11*m0s:m0s-1+(dent-1)*11*m0s
        r11 = p0+p1+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rw1/2);
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        for r11=p0+p1+2:p0+p1+p-1
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
        r11 = p0+p1+p;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

% Other teeth excluding s11 = 3*m0s+1 and 3*m0s+1+4*m0s and s11 = 5*m0s and 9*m0s
for dent = 1:2
    for s11 = 3*m0s+2+(dent-1)*4*m0s:5*m0s-1+(dent-1)*4*m0s
        r11 = p0+p1+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rw1/2);
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        for r11=p0+p1+2:p0+p1+p-1
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
        r11 = p0+p1+p;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

% Slot/tooth boundaries
for dent = 1:3
    s11 = m0s+1+(dent-1)*4*m0s;
    r11 = p0+p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p0+p1+2:p0+p1+2*p0s-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p0+p1+2*p0s;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rw1/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p0+p1+2*p0s+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rw1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p0+p1+2*p0s+2:p0+p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p0+p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

for dent = 1:3
    s11 = 3*m0s+(dent-1)*4*m0s;
    r11 = p0+p1+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p0+p1+2:p0+p1+2*p0s-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p0+p1+2*p0s;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/(Rv20/2+Rv1/2);
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rw1/2);
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    r11 = p0+p1+2*p0s+1;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rw1/2);
    MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    for r11 = p0+p1+2*p0s+2:p0+p1+p-1
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
    r11 = p0+p1+p;
    nnum(r11,s11)=(r11-1)*m+s11;
    MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
    MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
    M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
    M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
    M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
    M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
end

% Slots region and back-iron regions on top of the slots
for enc = 1:3
    for s11 = m0s+2+(enc-1)*4*m0s:3*m0s-1+(enc-1)*4*m0s
        r11 = p0+p1+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        for r11 = p0+p1+2:p0+p1+2*p0s-1
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw1;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
        r11 = p0+p1+2*p0s;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv1;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw1;
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/(Rw20/2+Rw1/2);
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        r11 = p0+p1+2*p0s+1;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/(Rw20/2+Rw1/2);
        MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        for r11 = p0+p1+2*p0s+2:p0+p1+p-1
            nnum(r11,s11)=(r11-1)*m+s11;
            MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
            MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
            MP(nnum(r11,s11),nnum(r11,s11)+m) = 1/Rw20;
            M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)+m)+MP(nnum(r11,s11),nnum(r11,s11)-m);
            M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
            M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
            M(nnum(r11,s11),nnum(r11,s11)+m)  = -MP(nnum(r11,s11),nnum(r11,s11)+m);
            M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
        end
        r11 = p0+p1+p;
        nnum(r11,s11)=(r11-1)*m+s11;
        MP(nnum(r11,s11),nnum(r11,s11)+1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-1) = 1/Rv20;
        MP(nnum(r11,s11),nnum(r11,s11)-m) = 1/Rw20;
        M(nnum(r11,s11),nnum(r11,s11))    = MP(nnum(r11,s11),nnum(r11,s11)-1)+MP(nnum(r11,s11),nnum(r11,s11)+1)+MP(nnum(r11,s11),nnum(r11,s11)-m);
        M(nnum(r11,s11),nnum(r11,s11)-1)  = -MP(nnum(r11,s11),nnum(r11,s11)-1);
        M(nnum(r11,s11),nnum(r11,s11)+1)  = -MP(nnum(r11,s11),nnum(r11,s11)+1);
        M(nnum(r11,s11),nnum(r11,s11)-m)  = -MP(nnum(r11,s11),nnum(r11,s11)-m);
    end
end

toc
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Excitation vector [Phi] %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FMMPM = (Br*h/mu0);
for i=m0m+1:(tm/l)+m0m
    E1(i)     = -FMMPM/2/(Rw21/2+Rw1/2);
    E1(i+120) =  FMMPM/2/(Rw21/2+Rw1/2);
    E2(i)     = -FMMPM*Rw21/Rw1/(Rw21+Rw1);
    E2(i+120) =  FMMPM*Rw21/Rw1/(Rw21+Rw1);
    E3(i)     =  Br*la*l/2;
    E3(i+120) = -Br*la*l/2;    
end
E1(240) = 0;
E2(240) = 0;
E3(240) = 0;

% Displacement vector
xd   = [0:l:tp/2];

for u=1:length(xd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Excitation vector [Phi] adaptation with xd %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:120
    r1      = p0;
    k1      = (r1-1)*m+s;
    E(k1)   = E1(s+u-1);
    E(k1+m) = E2(s+u-1);
    r2      = 2*p0;
    k2      = (r2-1)*m+s;
    E(k2)   = E3(s+u-1);
    E(k2+m) = E3(s+u-1);    
end
E(nn) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Solution vector  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

C=M\E';

% Computation of the flux in a phase coil
for i=1:(p0+p+p1)
    for j=1:m
        nnum2(i,j)=(i-1)*m+j;
        Ups(i,j) = C(nnum2(i,j));        
    end
end

for i=1:6*m0s
    nnum01=(p0+p1+20)*m+i;
    nnum02=(p0+p1+20+1)*m+i;
    nnum11=(p0+p1+20)*m+i+6*m0s;
    nnum12=(p0+p1+20+1)*m+i+6*m0s;
    Phip(i)=[Ups((p0+p1+20),i)-Ups((p0+p1+20+1),i)]*MP(nnum01,nnum02);
    Phim(i)=[Ups((p0+p1+20),i+6*m0s)-Ups((p0+p1+20+1),i+6*m0s)]*MP(nnum11,nnum12);
end
FluxC(u)=sum(Phip)-sum(Phim);

clear C E

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Flux variation storage %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save FichFlux tp l xd FluxC

figure
plot(xd*1000,FluxC)

toc