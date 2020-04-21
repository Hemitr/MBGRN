% Code used to plot Figure 15
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
    E1(i)       = -FMMPM/2/(Rw21/2+Rw1/2);
    E1(i+120)   = -FMMPM*Rw21/Rw1/(Rw21+Rw1);
    E2(i)       =  Br*la*l/2;
    E2(i+120)   =  Br*la*l/2;    
end
E1(240) = 0;
E2(240) = 0;

for s = 1:120
    r1      = p0;
    k1      = (r1-1)*m+s;
    E(k1)   = E1(s);
    E(k1+m) = E1(s+m);
    r2      = 2*p0;
    k2      = (r2-1)*m+s;
    E(k2)   = E2(s);
    E(k2+m) = E2(s+m);    
end
E(nn) = 0;
toc
tic

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Solution vector  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

C=M\E';
toc
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computation of the spatial distributions of B components %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y = hmbi + hm +*e/4

% Bx component
for i = 2:m-1
    nnum3 = (2*p0+1)*m+i;
    Bx(i) = mu0*(C(nnum3-1)-C(nnum3+1))/2/l;
end
nnum3   = (2*p0+1)*m+1;
Bx(1) = mu0*(-C(nnum3-1+m)-C(nnum3+1))/2/l;
Bx(m) = mu0*(C(nnum3-1+m-1)+C(nnum3))/2/l;

% By component
for i = 1:m
    nnum30 = (2*p0+1)*m+i;
    nnum31 = 2*p0*m+i;
    nnum32 = (2*p0+2)*m+i;
    By1(i) = (C(nnum31)-C(nnum30))*MP(nnum30,nnum31)/l/la;
    By2(i) = (C(nnum30)-C(nnum32))*MP(nnum30,nnum32)/l/la;
    By(i)  = [By1(i)+By2(i)]/2;
end

x = [l/2:l:tp-l/2];
figure
plot(x,Bx)
hold on
plot(x,By)

%%%%%%%%%%%%%%%%%%%%%%
%%% Flux2D results %%%
%%%%%%%%%%%%%%%%%%%%%%

Flux=[0.00000000E+00      0.00000000E+00      2.07500000E-02      0.00000000E+00     -3.19575521E-02     -1.00201052E-04
      2.50000000E-04      2.50000000E-04      2.07500000E-02      0.00000000E+00     -3.25334996E-02      3.19395652E-02
      5.00000000E-04      5.00000000E-04      2.07500000E-02      0.00000000E+00     -3.44108316E-02      6.39793315E-02
      7.50000000E-04      7.50000000E-04      2.07500000E-02      0.00000000E+00     -3.76626875E-02      1.00754928E-01
      1.00000000E-03      1.00000000E-03      2.07500000E-02      0.00000000E+00     -4.27410965E-02      1.37797868E-01
      1.25000000E-03      1.25000000E-03      2.07500000E-02      0.00000000E+00     -5.02096361E-02      1.86623002E-01
      1.50000000E-03      1.50000000E-03      2.07500000E-02      0.00000000E+00     -6.09038401E-02      2.35879435E-01
      1.75000000E-03      1.75000000E-03      2.07500000E-02      0.00000000E+00     -7.39715072E-02      3.09569954E-01
      2.00000000E-03      2.00000000E-03      2.07500000E-02      0.00000000E+00     -8.68213551E-02      3.83382156E-01
      2.25000000E-03      2.25000000E-03      2.07500000E-02      0.00000000E+00     -1.13205286E-01      4.92513472E-01
      2.50000000E-03      2.50000000E-03      2.07500000E-02      0.00000000E+00     -1.09113595E-01      6.12611472E-01
      2.75000000E-03      2.75000000E-03      2.07500000E-02      0.00000000E+00     -1.17228293E-01      7.19346709E-01
      3.00000000E-03      3.00000000E-03      2.07500000E-02      0.00000000E+00     -9.53613170E-02      8.43349865E-01
      3.25000000E-03      3.25000000E-03      2.07500000E-02      0.00000000E+00     -8.75891611E-02      9.20131250E-01
      3.50000000E-03      3.50000000E-03      2.07500000E-02      0.00000000E+00     -8.14943214E-02      1.00710958E+00
      3.75000000E-03      3.75000000E-03      2.07500000E-02      0.00000000E+00     -8.24972630E-02      1.08154748E+00
      4.00000000E-03      4.00000000E-03      2.07500000E-02      0.00000000E+00     -9.50114911E-02      1.16113667E+00
      4.25000000E-03      4.25000000E-03      2.07500000E-02      0.00000000E+00     -1.09990842E-01      1.28478919E+00
      4.50000000E-03      4.50000000E-03      2.07500000E-02      0.00000000E+00     -1.60782520E-01      1.38618025E+00
      4.75000000E-03      4.75000000E-03      2.07500000E-02      0.00000000E+00     -3.35330345E-01      1.50818566E+00
      5.00000000E-03      5.00000000E-03      2.07500000E-02      0.00000000E+00     -8.86276094E-01      1.35794675E+00
      5.25000000E-03      5.25000000E-03      2.07500000E-02      0.00000000E+00     -8.72267109E-01      1.12693130E+00
      5.50000000E-03      5.50000000E-03      2.07500000E-02      0.00000000E+00     -7.92789337E-01      8.95915860E-01
      5.75000000E-03      5.75000000E-03      2.07500000E-02      0.00000000E+00     -6.84370840E-01      8.04425867E-01
      6.00000000E-03      6.00000000E-03      2.07500000E-02      0.00000000E+00     -5.96876096E-01      7.32857153E-01
      6.25000000E-03      6.25000000E-03      2.07500000E-02      0.00000000E+00     -5.30285547E-01      7.00560086E-01
      6.50000000E-03      6.50000000E-03      2.07500000E-02      0.00000000E+00     -4.73366410E-01      6.69767412E-01
      6.75000000E-03      6.75000000E-03      2.07500000E-02      0.00000000E+00     -4.22199731E-01      6.55428637E-01
      7.00000000E-03      7.00000000E-03      2.07500000E-02      0.00000000E+00     -3.76934145E-01      6.39568569E-01
      7.25000000E-03      7.25000000E-03      2.07500000E-02      0.00000000E+00     -3.36135238E-01      6.30258955E-01
      7.50000000E-03      7.50000000E-03      2.07500000E-02      0.00000000E+00     -2.98798537E-01      6.21083458E-01
      7.75000000E-03      7.75000000E-03      2.07500000E-02      0.00000000E+00     -2.64053640E-01      6.15689778E-01
      8.00000000E-03      8.00000000E-03      2.07500000E-02      0.00000000E+00     -2.31508037E-01      6.10142995E-01
      8.25000000E-03      8.25000000E-03      2.07500000E-02      0.00000000E+00     -2.00676257E-01      6.06837811E-01
      8.50000000E-03      8.50000000E-03      2.07500000E-02      0.00000000E+00     -1.71255434E-01      6.03507345E-01
      8.75000000E-03      8.75000000E-03      2.07500000E-02      0.00000000E+00     -1.42925869E-01      6.01725817E-01
      9.00000000E-03      9.00000000E-03      2.07500000E-02      0.00000000E+00     -1.15474644E-01      5.99918940E-01
      9.25000000E-03      9.25000000E-03      2.07500000E-02      0.00000000E+00     -8.86742472E-02      5.99292959E-01
      9.50000000E-03      9.50000000E-03      2.07500000E-02      0.00000000E+00     -6.23496811E-02      5.98656855E-01
      9.75000000E-03      9.75000000E-03      2.07500000E-02      0.00000000E+00     -3.63200470E-02      5.99037624E-01
      1.00000000E-02      1.00000000E-02      2.07500000E-02      0.00000000E+00     -1.04261694E-02      5.99414354E-01
      1.02500000E-02      1.02500000E-02      2.07500000E-02      0.00000000E+00      1.54927562E-02      6.00778176E-01
      1.05000000E-02      1.05000000E-02      2.07500000E-02      0.00000000E+00      4.16004248E-02      6.02143713E-01
      1.07500000E-02      1.07500000E-02      2.07500000E-02      0.00000000E+00      6.80578792E-02      6.04586339E-01
      1.10000000E-02      1.10000000E-02      2.07500000E-02      0.00000000E+00      9.50543654E-02      6.07036631E-01
      1.12500000E-02      1.12500000E-02      2.07500000E-02      0.00000000E+00      1.22771977E-01      6.10795459E-01
      1.15000000E-02      1.15000000E-02      2.07500000E-02      0.00000000E+00      1.51454116E-01      6.14568576E-01
      1.17500000E-02      1.17500000E-02      2.07500000E-02      0.00000000E+00      1.81327577E-01      6.20085669E-01
      1.20000000E-02      1.20000000E-02      2.07500000E-02      0.00000000E+00      2.12740183E-01      6.25633711E-01
      1.22500000E-02      1.22500000E-02      2.07500000E-02      0.00000000E+00      2.46020275E-01      6.33753904E-01
      1.25000000E-02      1.25000000E-02      2.07500000E-02      0.00000000E+00      2.81702667E-01      6.41906137E-01
      1.27500000E-02      1.27500000E-02      2.07500000E-02      0.00000000E+00      3.20219959E-01      6.54117596E-01
      1.30000000E-02      1.30000000E-02      2.07500000E-02      0.00000000E+00      3.62537804E-01      6.66497746E-01
      1.32500000E-02      1.32500000E-02      2.07500000E-02      0.00000000E+00      4.09759851E-01      6.86662132E-01
      1.35000000E-02      1.35000000E-02      2.07500000E-02      0.00000000E+00      4.63498624E-01      7.06698991E-01
      1.37500000E-02      1.37500000E-02      2.07500000E-02      0.00000000E+00      5.23656458E-01      7.41852689E-01
      1.40000000E-02      1.40000000E-02      2.07500000E-02      0.00000000E+00      5.94628579E-01      7.78606259E-01
      1.42500000E-02      1.42500000E-02      2.07500000E-02      0.00000000E+00      6.88600190E-01      8.61571474E-01
      1.45000000E-02      1.45000000E-02      2.07500000E-02      0.00000000E+00      8.06262654E-01      9.43209327E-01
      1.47500000E-02      1.47500000E-02      2.07500000E-02      0.00000000E+00      8.94677893E-01      1.22072445E+00
      1.50000000E-02      1.50000000E-02      2.07500000E-02      0.00000000E+00      9.18266565E-01      1.47508724E+00
      1.52500000E-02      1.52500000E-02      2.07500000E-02      0.00000000E+00      3.33250460E-01      1.66210814E+00
      1.55000000E-02      1.55000000E-02      2.07500000E-02      0.00000000E+00      1.45451712E-01      1.42185225E+00
      1.57500000E-02      1.57500000E-02      2.07500000E-02      0.00000000E+00      8.75555218E-02      1.47556515E+00
      1.60000000E-02      1.60000000E-02      2.07500000E-02      0.00000000E+00      6.57801954E-02      1.39664192E+00
      1.62500000E-02      1.62500000E-02      2.07500000E-02      0.00000000E+00      4.42558336E-02      1.33058749E+00
      1.65000000E-02      1.65000000E-02      2.07500000E-02      0.00000000E+00      3.17335469E-02      1.28858121E+00
      1.67500000E-02      1.67500000E-02      2.07500000E-02      0.00000000E+00      2.44702942E-02      1.27015863E+00
      1.70000000E-02      1.70000000E-02      2.07500000E-02      0.00000000E+00      1.95194749E-02      1.24634254E+00
      1.72500000E-02      1.72500000E-02      2.07500000E-02      0.00000000E+00      1.55002471E-02      1.23037097E+00
      1.75000000E-02      1.75000000E-02      2.07500000E-02      0.00000000E+00      1.24914753E-02      1.21519551E+00
      1.77500000E-02      1.77500000E-02      2.07500000E-02      0.00000000E+00      1.01527848E-02      1.20553152E+00
      1.80000000E-02      1.80000000E-02      2.07500000E-02      0.00000000E+00      8.28255261E-03      1.19555637E+00
      1.82500000E-02      1.82500000E-02      2.07500000E-02      0.00000000E+00      6.71334195E-03      1.18896529E+00
      1.85000000E-02      1.85000000E-02      2.07500000E-02      0.00000000E+00      5.39294966E-03      1.18236057E+00
      1.87500000E-02      1.87500000E-02      2.07500000E-02      0.00000000E+00      4.24677631E-03      1.17822054E+00
      1.90000000E-02      1.90000000E-02      2.07500000E-02      0.00000000E+00      3.23460089E-03      1.17403965E+00
      1.92500000E-02      1.92500000E-02      2.07500000E-02      0.00000000E+00      2.31384040E-03      1.17177655E+00
      1.95000000E-02      1.95000000E-02      2.07500000E-02      0.00000000E+00      1.45997394E-03      1.16949814E+00
      1.97500000E-02      1.97500000E-02      2.07500000E-02      0.00000000E+00      6.46064575E-04      1.16887019E+00
      2.00000000E-02      2.00000000E-02      2.07500000E-02      0.00000000E+00     -1.48309381E-04      1.16823432E+00
      2.02500000E-02      2.02500000E-02      2.07500000E-02      0.00000000E+00     -9.44011376E-04      1.16916452E+00
      2.05000000E-02      2.05000000E-02      2.07500000E-02      0.00000000E+00     -1.76199216E-03      1.17009498E+00
      2.07500000E-02      2.07500000E-02      2.07500000E-02      0.00000000E+00     -2.62273019E-03      1.17267039E+00
      2.10000000E-02      2.10000000E-02      2.07500000E-02      0.00000000E+00     -3.55341956E-03      1.17525426E+00
      2.12500000E-02      2.12500000E-02      2.07500000E-02      0.00000000E+00     -4.57879147E-03      1.17974764E+00
      2.15000000E-02      2.15000000E-02      2.07500000E-02      0.00000000E+00     -5.74195578E-03      1.18425695E+00
      2.17500000E-02      2.17500000E-02      2.07500000E-02      0.00000000E+00     -7.08361191E-03      1.19118624E+00
      2.20000000E-02      2.20000000E-02      2.07500000E-02      0.00000000E+00     -8.67935870E-03      1.19815730E+00
      2.22500000E-02      2.22500000E-02      2.07500000E-02      0.00000000E+00     -1.05821276E-02      1.20854282E+00
      2.25000000E-02      2.25000000E-02      2.07500000E-02      0.00000000E+00     -1.29615584E-02      1.21894269E+00
      2.27500000E-02      2.27500000E-02      2.07500000E-02      0.00000000E+00     -1.60217756E-02      1.23431822E+00
      2.30000000E-02      2.30000000E-02      2.07500000E-02      0.00000000E+00     -2.01076375E-02      1.25000892E+00
      2.32500000E-02      2.32500000E-02      2.07500000E-02      0.00000000E+00     -2.51401861E-02      1.27528829E+00
      2.35000000E-02      2.35000000E-02      2.07500000E-02      0.00000000E+00     -3.25168950E-02      1.29976407E+00
      2.37500000E-02      2.37500000E-02      2.07500000E-02      0.00000000E+00     -4.52171274E-02      1.33725089E+00
      2.40000000E-02      2.40000000E-02      2.07500000E-02      0.00000000E+00     -6.70291776E-02      1.38018993E+00
      2.42500000E-02      2.42500000E-02      2.07500000E-02      0.00000000E+00     -8.91137956E-02      1.48486970E+00
      2.45000000E-02      2.45000000E-02      2.07500000E-02      0.00000000E+00     -1.47742369E-01      1.56524927E+00
      2.47500000E-02      2.47500000E-02      2.07500000E-02      0.00000000E+00     -3.37701707E-01      1.67488363E+00
      2.50000000E-02      2.50000000E-02      2.07500000E-02      0.00000000E+00     -9.29262562E-01      1.91856388E+00
      2.52500000E-02      2.52500000E-02      2.07500000E-02      0.00000000E+00     -9.06622597E-01      1.23045515E+00
      2.55000000E-02      2.55000000E-02      2.07500000E-02      0.00000000E+00     -8.18084605E-01      9.73637272E-01
      2.57500000E-02      2.57500000E-02      2.07500000E-02      0.00000000E+00     -6.99875717E-01      8.68100336E-01
      2.60000000E-02      2.60000000E-02      2.07500000E-02      0.00000000E+00     -6.05488962E-01      7.85778823E-01
      2.62500000E-02      2.62500000E-02      2.07500000E-02      0.00000000E+00     -5.34278859E-01      7.47477922E-01
      2.65000000E-02      2.65000000E-02      2.07500000E-02      0.00000000E+00     -4.73952278E-01      7.10540023E-01
      2.67500000E-02      2.67500000E-02      2.07500000E-02      0.00000000E+00     -4.20085713E-01      6.92165476E-01
      2.70000000E-02      2.70000000E-02      2.07500000E-02      0.00000000E+00     -3.72782241E-01      6.72171630E-01
      2.72500000E-02      2.72500000E-02      2.07500000E-02      0.00000000E+00     -3.30419099E-01      6.59844590E-01
      2.75000000E-02      2.75000000E-02      2.07500000E-02      0.00000000E+00     -2.91883652E-01      6.47647197E-01
      2.77500000E-02      2.77500000E-02      2.07500000E-02      0.00000000E+00     -2.56203064E-01      6.39960096E-01
      2.80000000E-02      2.80000000E-02      2.07500000E-02      0.00000000E+00     -2.22940894E-01      6.32102923E-01
      2.82500000E-02      2.82500000E-02      2.07500000E-02      0.00000000E+00     -1.91557751E-01      6.26987727E-01
      2.85000000E-02      2.85000000E-02      2.07500000E-02      0.00000000E+00     -1.61721430E-01      6.21840772E-01
      2.87500000E-02      2.87500000E-02      2.07500000E-02      0.00000000E+00     -1.33080264E-01      6.18610592E-01
      2.90000000E-02      2.90000000E-02      2.07500000E-02      0.00000000E+00     -1.05403565E-01      6.15349920E-01
      2.92500000E-02      2.92500000E-02      2.07500000E-02      0.00000000E+00     -7.84441007E-02      6.13550289E-01
      2.95000000E-02      2.95000000E-02      2.07500000E-02      0.00000000E+00     -5.20150841E-02      6.11737232E-01
      2.97500000E-02      2.97500000E-02      2.07500000E-02      0.00000000E+00     -2.59227930E-02      6.11161276E-01
      3.00000000E-02      3.00000000E-02      2.07500000E-02      0.00000000E+00     -3.13717976E-15      6.10578890E-01
      3.02500000E-02      3.02500000E-02      2.07500000E-02      0.00000000E+00      2.59227930E-02      6.11161276E-01
      3.05000000E-02      3.05000000E-02      2.07500000E-02      0.00000000E+00      5.20150841E-02      6.11743662E-01
      3.07500000E-02      3.07500000E-02      2.07500000E-02      0.00000000E+00      7.84441007E-02      6.13550289E-01
      3.10000000E-02      3.10000000E-02      2.07500000E-02      0.00000000E+00      1.05403565E-01      6.15363346E-01
      3.12500000E-02      3.12500000E-02      2.07500000E-02      0.00000000E+00      1.33080264E-01      6.18610592E-01
      3.15000000E-02      3.15000000E-02      2.07500000E-02      0.00000000E+00      1.61721430E-01      6.21871264E-01
      3.17500000E-02      3.17500000E-02      2.07500000E-02      0.00000000E+00      1.91557751E-01      6.26987727E-01
      3.20000000E-02      3.20000000E-02      2.07500000E-02      0.00000000E+00      2.22940894E-01      6.32134682E-01
      3.22500000E-02      3.22500000E-02      2.07500000E-02      0.00000000E+00      2.56203064E-01      6.39960096E-01
      3.25000000E-02      3.25000000E-02      2.07500000E-02      0.00000000E+00      2.91883652E-01      6.47817270E-01
      3.27500000E-02      3.27500000E-02      2.07500000E-02      0.00000000E+00      3.30419099E-01      6.59844590E-01
      3.30000000E-02      3.30000000E-02      2.07500000E-02      0.00000000E+00      3.72782241E-01      6.72041984E-01
      3.32500000E-02      3.32500000E-02      2.07500000E-02      0.00000000E+00      4.20085713E-01      6.92165476E-01
      3.35000000E-02      3.35000000E-02      2.07500000E-02      0.00000000E+00      4.73952278E-01      7.12159322E-01
      3.37500000E-02      3.37500000E-02      2.07500000E-02      0.00000000E+00      5.34278859E-01      7.47477922E-01
      3.40000000E-02      3.40000000E-02      2.07500000E-02      0.00000000E+00      6.05488962E-01      7.84415822E-01
      3.42500000E-02      3.42500000E-02      2.07500000E-02      0.00000000E+00      6.99875717E-01      8.68100336E-01
      3.45000000E-02      3.45000000E-02      2.07500000E-02      0.00000000E+00      8.18084605E-01      9.50421848E-01
      3.47500000E-02      3.47500000E-02      2.07500000E-02      0.00000000E+00      9.06622597E-01      1.23045515E+00
      3.50000000E-02      3.50000000E-02      2.07500000E-02      0.00000000E+00      9.29262562E-01      1.48727303E+00
      3.52500000E-02      3.52500000E-02      2.07500000E-02      0.00000000E+00      3.37701707E-01      1.67488363E+00
      3.55000000E-02      3.55000000E-02      2.07500000E-02      0.00000000E+00      1.47742369E-01      1.43120337E+00
      3.57500000E-02      3.57500000E-02      2.07500000E-02      0.00000000E+00      8.91137956E-02      1.48486970E+00
      3.60000000E-02      3.60000000E-02      2.07500000E-02      0.00000000E+00      6.70291776E-02      1.40449014E+00
      3.62500000E-02      3.62500000E-02      2.07500000E-02      0.00000000E+00      4.52171274E-02      1.33725089E+00
      3.65000000E-02      3.65000000E-02      2.07500000E-02      0.00000000E+00      3.25168950E-02      1.29431185E+00
      3.67500000E-02      3.67500000E-02      2.07500000E-02      0.00000000E+00      2.51401861E-02      1.27528829E+00
      3.70000000E-02      3.70000000E-02      2.07500000E-02      0.00000000E+00      2.01076375E-02      1.25081251E+00
      3.72500000E-02      3.72500000E-02      2.07500000E-02      0.00000000E+00      1.60217756E-02      1.23431822E+00
      3.75000000E-02      3.75000000E-02      2.07500000E-02      0.00000000E+00      1.29615584E-02      1.21862752E+00
      3.77500000E-02      3.77500000E-02      2.07500000E-02      0.00000000E+00      1.05821276E-02      1.20854282E+00
      3.80000000E-02      3.80000000E-02      2.07500000E-02      0.00000000E+00      8.67935870E-03      1.19814295E+00
      3.82500000E-02      3.82500000E-02      2.07500000E-02      0.00000000E+00      7.08361191E-03      1.19118624E+00
      3.85000000E-02      3.85000000E-02      2.07500000E-02      0.00000000E+00      5.74195578E-03      1.18421519E+00
      3.87500000E-02      3.87500000E-02      2.07500000E-02      0.00000000E+00      4.57879147E-03      1.17974764E+00
      3.90000000E-02      3.90000000E-02      2.07500000E-02      0.00000000E+00      3.55341956E-03      1.17523834E+00
      3.92500000E-02      3.92500000E-02      2.07500000E-02      0.00000000E+00      2.62273019E-03      1.17267039E+00
      3.95000000E-02      3.95000000E-02      2.07500000E-02      0.00000000E+00      1.76199216E-03      1.17008652E+00
      3.97500000E-02      3.97500000E-02      2.07500000E-02      0.00000000E+00      9.44011376E-04      1.16916452E+00
      4.00000000E-02      4.00000000E-02      2.07500000E-02      0.00000000E+00      1.48309381E-04      1.16823406E+00
      4.02500000E-02      4.02500000E-02      2.07500000E-02      0.00000000E+00     -6.46064575E-04      1.16887019E+00
      4.05000000E-02      4.05000000E-02      2.07500000E-02      0.00000000E+00     -1.45997394E-03      1.16950605E+00
      4.07500000E-02      4.07500000E-02      2.07500000E-02      0.00000000E+00     -2.31384040E-03      1.17177655E+00
      4.10000000E-02      4.10000000E-02      2.07500000E-02      0.00000000E+00     -3.23460089E-03      1.17405496E+00
      4.12500000E-02      4.12500000E-02      2.07500000E-02      0.00000000E+00     -4.24677631E-03      1.17822054E+00
      4.15000000E-02      4.15000000E-02      2.07500000E-02      0.00000000E+00     -5.39294966E-03      1.18240142E+00
      4.17500000E-02      4.17500000E-02      2.07500000E-02      0.00000000E+00     -6.71334195E-03      1.18896529E+00
      4.20000000E-02      4.20000000E-02      2.07500000E-02      0.00000000E+00     -8.28255261E-03      1.19557002E+00
      4.22500000E-02      4.22500000E-02      2.07500000E-02      0.00000000E+00     -1.01527848E-02      1.20553152E+00
      4.25000000E-02      4.25000000E-02      2.07500000E-02      0.00000000E+00     -1.24914753E-02      1.21550668E+00
      4.27500000E-02      4.27500000E-02      2.07500000E-02      0.00000000E+00     -1.55002471E-02      1.23037097E+00
      4.30000000E-02      4.30000000E-02      2.07500000E-02      0.00000000E+00     -1.95194749E-02      1.24554642E+00
      4.32500000E-02      4.32500000E-02      2.07500000E-02      0.00000000E+00     -2.44702942E-02      1.27015863E+00
      4.35000000E-02      4.35000000E-02      2.07500000E-02      0.00000000E+00     -3.17335469E-02      1.29397472E+00
      4.37500000E-02      4.37500000E-02      2.07500000E-02      0.00000000E+00     -4.42558336E-02      1.33058749E+00
      4.40000000E-02      4.40000000E-02      2.07500000E-02      0.00000000E+00     -6.57801954E-02      1.37259377E+00
      4.42500000E-02      4.42500000E-02      2.07500000E-02      0.00000000E+00     -8.75555218E-02      1.47556515E+00
      4.45000000E-02      4.45000000E-02      2.07500000E-02      0.00000000E+00     -1.45451712E-01      1.55448837E+00
      4.47500000E-02      4.47500000E-02      2.07500000E-02      0.00000000E+00     -3.33250460E-01      1.66210814E+00
      4.50000000E-02      4.50000000E-02      2.07500000E-02      0.00000000E+00     -9.18266565E-01      1.90236403E+00
      4.52500000E-02      4.52500000E-02      2.07500000E-02      0.00000000E+00     -8.94677893E-01      1.22072445E+00
      4.55000000E-02      4.55000000E-02      2.07500000E-02      0.00000000E+00     -8.06262654E-01      9.66361658E-01
      4.57500000E-02      4.57500000E-02      2.07500000E-02      0.00000000E+00     -6.88600190E-01      8.61571474E-01
      4.60000000E-02      4.60000000E-02      2.07500000E-02      0.00000000E+00     -5.94628579E-01      7.79933621E-01
      4.62500000E-02      4.62500000E-02      2.07500000E-02      0.00000000E+00     -5.23656458E-01      7.41852689E-01
      4.65000000E-02      4.65000000E-02      2.07500000E-02      0.00000000E+00     -4.63498624E-01      7.05099118E-01
      4.67500000E-02      4.67500000E-02      2.07500000E-02      0.00000000E+00     -4.09759851E-01      6.86662132E-01
      4.70000000E-02      4.70000000E-02      2.07500000E-02      0.00000000E+00     -3.62537804E-01      6.66625273E-01
      4.72500000E-02      4.72500000E-02      2.07500000E-02      0.00000000E+00     -3.20219959E-01      6.54117596E-01
      4.75000000E-02      4.75000000E-02      2.07500000E-02      0.00000000E+00     -2.81702667E-01      6.41737446E-01
      4.77500000E-02      4.77500000E-02      2.07500000E-02      0.00000000E+00     -2.46020275E-01      6.33753904E-01
      4.80000000E-02      4.80000000E-02      2.07500000E-02      0.00000000E+00     -2.12740183E-01      6.25601672E-01
      4.82500000E-02      4.82500000E-02      2.07500000E-02      0.00000000E+00     -1.81327577E-01      6.20085669E-01
      4.85000000E-02      4.85000000E-02      2.07500000E-02      0.00000000E+00     -1.51454116E-01      6.14537627E-01
      4.87500000E-02      4.87500000E-02      2.07500000E-02      0.00000000E+00     -1.22771977E-01      6.10795459E-01
      4.90000000E-02      4.90000000E-02      2.07500000E-02      0.00000000E+00     -9.50543654E-02      6.07022342E-01
      4.92500000E-02      4.92500000E-02      2.07500000E-02      0.00000000E+00     -6.80578792E-02      6.04586339E-01
      4.95000000E-02      4.95000000E-02      2.07500000E-02      0.00000000E+00     -4.16004248E-02      6.02136047E-01
      4.97500000E-02      4.97500000E-02      2.07500000E-02      0.00000000E+00     -1.54927562E-02      6.00778176E-01
      5.00000000E-02      5.00000000E-02      2.07500000E-02      0.00000000E+00      1.04261694E-02      5.99412638E-01
      5.02500000E-02      5.02500000E-02      2.07500000E-02      0.00000000E+00      3.63200470E-02      5.99037624E-01
      5.05000000E-02      5.05000000E-02      2.07500000E-02      0.00000000E+00      6.23496811E-02      5.98660894E-01
      5.07500000E-02      5.07500000E-02      2.07500000E-02      0.00000000E+00      8.86742472E-02      5.99292959E-01
      5.10000000E-02      5.10000000E-02      2.07500000E-02      0.00000000E+00      1.15474644E-01      5.99929062E-01
      5.12500000E-02      5.12500000E-02      2.07500000E-02      0.00000000E+00      1.42925869E-01      6.01725817E-01
      5.15000000E-02      5.15000000E-02      2.07500000E-02      0.00000000E+00      1.71255434E-01      6.03532693E-01
      5.17500000E-02      5.17500000E-02      2.07500000E-02      0.00000000E+00      2.00676257E-01      6.06837811E-01
      5.20000000E-02      5.20000000E-02      2.07500000E-02      0.00000000E+00      2.31508037E-01      6.10168276E-01
      5.22500000E-02      5.22500000E-02      2.07500000E-02      0.00000000E+00      2.64053640E-01      6.15689778E-01
      5.25000000E-02      5.25000000E-02      2.07500000E-02      0.00000000E+00      2.98798537E-01      6.21236560E-01
      5.27500000E-02      5.27500000E-02      2.07500000E-02      0.00000000E+00      3.36135238E-01      6.30258955E-01
      5.30000000E-02      5.30000000E-02      2.07500000E-02      0.00000000E+00      3.76934145E-01      6.39434453E-01
      5.32500000E-02      5.32500000E-02      2.07500000E-02      0.00000000E+00      4.22199731E-01      6.55428637E-01
      5.35000000E-02      5.35000000E-02      2.07500000E-02      0.00000000E+00      4.73366410E-01      6.71288706E-01
      5.37500000E-02      5.37500000E-02      2.07500000E-02      0.00000000E+00      5.30285547E-01      7.00560086E-01
      5.40000000E-02      5.40000000E-02      2.07500000E-02      0.00000000E+00      5.96876096E-01      7.31352759E-01
      5.42500000E-02      5.42500000E-02      2.07500000E-02      0.00000000E+00      6.84370840E-01      8.04425867E-01
      5.45000000E-02      5.45000000E-02      2.07500000E-02      0.00000000E+00      7.92789337E-01      8.75994581E-01
      5.47500000E-02      5.47500000E-02      2.07500000E-02      0.00000000E+00      8.72267109E-01      1.12693130E+00
      5.50000000E-02      5.50000000E-02      2.07500000E-02      0.00000000E+00      8.86276094E-01      1.35794675E+00
      5.52500000E-02      5.52500000E-02      2.07500000E-02      0.00000000E+00      3.35330345E-01      1.50818566E+00
      5.55000000E-02      5.55000000E-02      2.07500000E-02      0.00000000E+00      1.60782520E-01      1.26099412E+00
      5.57500000E-02      5.57500000E-02      2.07500000E-02      0.00000000E+00      1.09990842E-01      1.28478919E+00
      5.60000000E-02      5.60000000E-02      2.07500000E-02      0.00000000E+00      9.50114911E-02      1.18339813E+00
      5.62500000E-02      5.62500000E-02      2.07500000E-02      0.00000000E+00      8.24972630E-02      1.08154748E+00
      5.65000000E-02      5.65000000E-02      2.07500000E-02      0.00000000E+00      8.14943214E-02      1.00195829E+00
      5.67500000E-02      5.67500000E-02      2.07500000E-02      0.00000000E+00      8.75891611E-02      9.20131250E-01
      5.70000000E-02      5.70000000E-02      2.07500000E-02      0.00000000E+00      9.53613170E-02      8.33152921E-01
      5.72500000E-02      5.72500000E-02      2.07500000E-02      0.00000000E+00      1.17228293E-01      7.19346709E-01
      5.75000000E-02      5.75000000E-02      2.07500000E-02      0.00000000E+00      1.09113595E-01      6.12611472E-01
      5.77500000E-02      5.77500000E-02      2.07500000E-02      0.00000000E+00      1.13205286E-01      4.92513472E-01
      5.80000000E-02      5.80000000E-02      2.07500000E-02      0.00000000E+00      8.68213551E-02      3.72415471E-01
      5.82500000E-02      5.82500000E-02      2.07500000E-02      0.00000000E+00      7.39715072E-02      3.09569954E-01
      5.85000000E-02      5.85000000E-02      2.07500000E-02      0.00000000E+00      6.09038401E-02      2.35757751E-01
      5.87500000E-02      5.87500000E-02      2.07500000E-02      0.00000000E+00      5.02096361E-02      1.86623002E-01
      5.90000000E-02      5.90000000E-02      2.07500000E-02      0.00000000E+00      4.27410965E-02      1.37366570E-01
      5.92500000E-02      5.92500000E-02      2.07500000E-02      0.00000000E+00      3.76626875E-02      1.00754928E-01
      5.95000000E-02      5.95000000E-02      2.07500000E-02      0.00000000E+00      3.44108316E-02      6.37119886E-02
      5.97500000E-02      5.97500000E-02      2.07500000E-02      0.00000000E+00      3.25334996E-02      3.19395652E-02
      6.00000000E-02      6.00000000E-02      2.07500000E-02      0.00000000E+00      3.19575521E-02     -1.00201052E-04];
figure
plot(Flux(:,1),Flux(:,5))
hold on
plot(Flux(:,1),Flux(:,6))

%%%%%%%%%%%%%%%%%%
%%% Comparison %%%
%%%%%%%%%%%%%%%%%%

figure
plot(Flux(:,1),Flux(:,5),'d')
hold on
plot(Flux(:,1),Flux(:,6),'o')
plot(x,Bx,'--')
plot(x,By)

% For clarity reasons the number of points in the FEM results is reduced
for i=1:81
    Flux1(i)=Flux(1+(i-1)*3,1);
    Flux2(i)=Flux(1+(i-1)*3,5);
    Flux3(i)=Flux(1+(i-1)*3,6);
end
figure
plot(Flux1*1000,Flux3,'o')
hold on
plot(x*1000,By)
plot(Flux1*1000,Flux2,'d')
plot(x*1000,Bx,'--')

toc