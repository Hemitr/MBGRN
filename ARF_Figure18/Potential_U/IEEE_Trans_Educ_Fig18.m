% Code used to plot Figure 18
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
Br   = 0;            % PM remanent induction (residual induction) (T)
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

% Terms coming from currents distribution in the slots

J   = 5;         % Maximum current denity (A/mm�)
NI  = 20*10*J;
wt = 0;
NIA = NI*cos(wt-2*pi/3);
NIB = NI*cos(wt+2*pi/3);
NIC = NI*cos(wt);
xenc = [l/2:l:ws-l/2];
for i=1:length(xenc)
    NIA0(i) =  NIA/2 - NIA*xenc(i)/ws;
    NIB0(i) =  NIB/2 - NIB*xenc(i)/ws;
    NIC0(i) = -NIC/2 + NIC*xenc(i)/ws;
end
for i = 1:m
    % Phase A
    if i < 9*m0s+1
        NIA1(i) = NIA/2;
    else
        if i > 11*m0s
            NIA1(i) = -NIA/2;
        else
            NIA1(i) = NIA0(i-9*m0s);
        end
    end
    % Phase B
    if i < m0s+1
        NIB1(i) = NIB/2;
    else
        if i > 3*m0s
            NIB1(i) = -NIB/2;
        else
            NIB1(i) = NIB0(i-m0s);
        end
    end
    % Phase C
    if i < 5*m0s+1
        NIC1(i) = -NIC/2;
    else
        if i > 7*m0s
            NIC1(i) = NIC/2;
        else
            NIC1(i) = NIC0(i-5*m0s);
        end
    end
    NIT(i) = NIA1(i)+NIB1(i)+NIC1(i);
    FMMT(i) = NIT(i)/4/p0s;
    r = p0+p1;
    nnum=(r-1)*m+i;
    E(nnum) = -FMMT(i)*MP(nnum,nnum+m);
    r = p0+p1+1;
    nnum=(r-1)*m+i;
    E(nnum) = -2*FMMT(i)*MP(nnum,nnum+m)+FMMT(i)*MP(nnum,nnum-m);
    r = p0+p1+2*p0s;
    nnum=(r-1)*m+i;
    E(nnum) = 2*FMMT(i)*MP(nnum,nnum-m)-FMMT(i)*MP(nnum,nnum+m);
    r = p0+p1+2*p0s+1;
    nnum=(r-1)*m+i;
    E(nnum) = FMMT(i)*MP(nnum,nnum-m);
end

% Terms coming from the PM

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

for s = 1:120
    r1      = p0;
    k1      = (r1-1)*m+s;
    E(k1)   = E1(s);
    E(k1+m) = E2(s);
    r2      = 2*p0;
    k2      = (r2-1)*m+s;
    E(k2)   = E3(s);
    E(k2+m) = E3(s);    
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
    By2(i) = (C(nnum30)-C(nnum32)+FMMT(i))*MP(nnum30,nnum32)/l/la;
    By(i)  = [By1(i)+By2(i)]/2;
end

x = [l/2:l:tp-l/2];
figure
plot(x,Bx)
hold on
plot(x,By)
toc
tic

%%%%%%%%%%%%%%%%%%%%%%
%%% Flux2D results %%%
%%%%%%%%%%%%%%%%%%%%%%

Flux=[0.00000000E+00      0.00000000E+00      2.07500000E-02      0.00000000E+00      5.40499075E-17     -1.31795861E-01
      2.50000000E-04      2.50000000E-04      2.07500000E-02      0.00000000E+00      1.59471291E-04     -1.31953202E-01
      5.00000000E-04      5.00000000E-04      2.07500000E-02      0.00000000E+00      3.22597660E-04     -1.32110543E-01
      7.50000000E-04      7.50000000E-04      2.07500000E-02      0.00000000E+00      4.92957702E-04     -1.32595439E-01
      1.00000000E-03      1.00000000E-03      2.07500000E-02      0.00000000E+00      6.75213719E-04     -1.33081762E-01
      1.25000000E-03      1.25000000E-03      2.07500000E-02      0.00000000E+00      8.73638577E-04     -1.33940646E-01
      1.50000000E-03      1.50000000E-03      2.07500000E-02      0.00000000E+00      1.09544707E-03     -1.34802306E-01
      1.75000000E-03      1.75000000E-03      2.07500000E-02      0.00000000E+00      1.34743570E-03     -1.36123580E-01
      2.00000000E-03      2.00000000E-03      2.07500000E-02      0.00000000E+00      1.64202685E-03     -1.37451738E-01
      2.25000000E-03      2.25000000E-03      2.07500000E-02      0.00000000E+00      1.98820651E-03     -1.39405796E-01
      2.50000000E-03      2.50000000E-03      2.07500000E-02      0.00000000E+00      2.41307209E-03     -1.41363245E-01
      2.75000000E-03      2.75000000E-03      2.07500000E-02      0.00000000E+00      2.94843610E-03     -1.44209244E-01
      3.00000000E-03      3.00000000E-03      2.07500000E-02      0.00000000E+00      3.64912251E-03     -1.47103256E-01
      3.25000000E-03      3.25000000E-03      2.07500000E-02      0.00000000E+00      4.50758550E-03     -1.51617618E-01
      3.50000000E-03      3.50000000E-03      2.07500000E-02      0.00000000E+00      5.73766997E-03     -1.56018535E-01
      3.75000000E-03      3.75000000E-03      2.07500000E-02      0.00000000E+00      7.78104945E-03     -1.62648830E-01
      4.00000000E-03      4.00000000E-03      2.07500000E-02      0.00000000E+00      1.12184587E-02     -1.70079366E-01
      4.25000000E-03      4.25000000E-03      2.07500000E-02      0.00000000E+00      1.47874531E-02     -1.87068544E-01
      4.50000000E-03      4.50000000E-03      2.07500000E-02      0.00000000E+00      2.38804704E-02     -2.00533778E-01
      4.75000000E-03      4.75000000E-03      2.07500000E-02      0.00000000E+00      5.23532837E-02     -2.19453963E-01
      5.00000000E-03      5.00000000E-03      2.07500000E-02      0.00000000E+00      1.40197981E-01     -1.97918111E-01
      5.25000000E-03      5.25000000E-03      2.07500000E-02      0.00000000E+00      1.43221185E-01     -1.59844031E-01
      5.50000000E-03      5.50000000E-03      2.07500000E-02      0.00000000E+00      1.34504997E-01     -1.21769950E-01
      5.75000000E-03      5.75000000E-03      2.07500000E-02      0.00000000E+00      1.20803498E-01     -1.06214075E-01
      6.00000000E-03      6.00000000E-03      2.07500000E-02      0.00000000E+00      1.09876741E-01     -9.30438652E-02
      6.25000000E-03      6.25000000E-03      2.07500000E-02      0.00000000E+00      1.01927601E-01     -8.60984953E-02
      6.50000000E-03      6.50000000E-03      2.07500000E-02      0.00000000E+00      9.52883226E-02     -7.94672792E-02
      6.75000000E-03      6.75000000E-03      2.07500000E-02      0.00000000E+00      8.94208995E-02     -7.55926893E-02
      7.00000000E-03      7.00000000E-03      2.07500000E-02      0.00000000E+00      8.43705579E-02     -7.14667000E-02
      7.25000000E-03      7.25000000E-03      2.07500000E-02      0.00000000E+00      7.99446839E-02     -6.84900412E-02
      7.50000000E-03      7.50000000E-03      2.07500000E-02      0.00000000E+00      7.60036031E-02     -6.55336790E-02
      7.75000000E-03      7.75000000E-03      2.07500000E-02      0.00000000E+00      7.24295550E-02     -6.32567766E-02
      8.00000000E-03      8.00000000E-03      2.07500000E-02      0.00000000E+00      6.91725956E-02     -6.09539714E-02
      8.25000000E-03      8.25000000E-03      2.07500000E-02      0.00000000E+00      6.61691853E-02     -5.90570663E-02
      8.50000000E-03      8.50000000E-03      2.07500000E-02      0.00000000E+00      6.33815241E-02     -5.71548044E-02
      8.75000000E-03      8.75000000E-03      2.07500000E-02      0.00000000E+00      6.07694655E-02     -5.55222235E-02
      9.00000000E-03      9.00000000E-03      2.07500000E-02      0.00000000E+00      5.83086016E-02     -5.38843761E-02
      9.25000000E-03      9.25000000E-03      2.07500000E-02      0.00000000E+00      5.59722971E-02     -5.24324671E-02
      9.50000000E-03      9.50000000E-03      2.07500000E-02      0.00000000E+00      5.37424809E-02     -5.09777010E-02
      9.75000000E-03      9.75000000E-03      2.07500000E-02      0.00000000E+00      5.15999995E-02     -4.96562255E-02
      1.00000000E-02      1.00000000E-02      2.07500000E-02      0.00000000E+00      4.95306383E-02     -4.83327328E-02
      1.02500000E-02      1.02500000E-02      2.07500000E-02      0.00000000E+00      4.75195536E-02     -4.71070547E-02
      1.05000000E-02      1.05000000E-02      2.07500000E-02      0.00000000E+00      4.55544261E-02     -4.58799979E-02
      1.07500000E-02      1.07500000E-02      2.07500000E-02      0.00000000E+00      4.36226718E-02     -4.47261085E-02
      1.10000000E-02      1.10000000E-02      2.07500000E-02      0.00000000E+00      4.17123963E-02     -4.35712349E-02
      1.12500000E-02      1.12500000E-02      2.07500000E-02      0.00000000E+00      3.98117354E-02     -4.24722608E-02
      1.15000000E-02      1.15000000E-02      2.07500000E-02      0.00000000E+00      3.79076974E-02     -4.13725862E-02
      1.17500000E-02      1.17500000E-02      2.07500000E-02      0.00000000E+00      3.59877128E-02     -4.03174911E-02
      1.20000000E-02      1.20000000E-02      2.07500000E-02      0.00000000E+00      3.40355663E-02     -3.92620136E-02
      1.22500000E-02      1.22500000E-02      2.07500000E-02      0.00000000E+00      3.20358211E-02     -3.82467718E-02
      1.25000000E-02      1.25000000E-02      2.07500000E-02      0.00000000E+00      2.99655774E-02     -3.72312733E-02
      1.27500000E-02      1.27500000E-02      2.07500000E-02      0.00000000E+00      2.78047880E-02     -3.62620698E-02
      1.30000000E-02      1.30000000E-02      2.07500000E-02      0.00000000E+00      2.55154041E-02     -3.52940908E-02
      1.32500000E-02      1.32500000E-02      2.07500000E-02      0.00000000E+00      2.30591024E-02     -3.44095826E-02
      1.35000000E-02      1.35000000E-02      2.07500000E-02      0.00000000E+00      2.03684903E-02     -3.35248731E-02
      1.37500000E-02      1.37500000E-02      2.07500000E-02      0.00000000E+00      1.74135326E-02     -3.28405986E-02
      1.40000000E-02      1.40000000E-02      2.07500000E-02      0.00000000E+00      1.40243040E-02     -3.21624516E-02
      1.42500000E-02      1.42500000E-02      2.07500000E-02      0.00000000E+00      9.90727797E-03     -3.20913104E-02
      1.45000000E-02      1.45000000E-02      2.07500000E-02      0.00000000E+00      4.74745652E-03     -3.21207194E-02
      1.47500000E-02      1.47500000E-02      2.07500000E-02      0.00000000E+00     -6.98463099E-04     -3.64060929E-02
      1.50000000E-02      1.50000000E-02      2.07500000E-02      0.00000000E+00     -7.45890295E-03     -3.93536564E-02
      1.52500000E-02      1.52500000E-02      2.07500000E-02      0.00000000E+00     -1.13320015E-03     -4.80394138E-02
      1.55000000E-02      1.55000000E-02      2.07500000E-02      0.00000000E+00      6.43694458E-04     -4.79686813E-02
      1.57500000E-02      1.57500000E-02      2.07500000E-02      0.00000000E+00      8.97254630E-04     -5.07412459E-02
      1.60000000E-02      1.60000000E-02      2.07500000E-02      0.00000000E+00      8.62453757E-04     -5.17034156E-02
      1.62500000E-02      1.62500000E-02      2.07500000E-02      0.00000000E+00      9.19597373E-04     -5.22852654E-02
      1.65000000E-02      1.65000000E-02      2.07500000E-02      0.00000000E+00      9.23470688E-04     -5.32128220E-02
      1.67500000E-02      1.67500000E-02      2.07500000E-02      0.00000000E+00      8.96096470E-04     -5.41697913E-02
      1.70000000E-02      1.70000000E-02      2.07500000E-02      0.00000000E+00      8.70294805E-04     -5.50594968E-02
      1.72500000E-02      1.72500000E-02      2.07500000E-02      0.00000000E+00      8.54782289E-04     -5.58905007E-02
      1.75000000E-02      1.75000000E-02      2.07500000E-02      0.00000000E+00      8.44821365E-04     -5.67366371E-02
      1.77500000E-02      1.77500000E-02      2.07500000E-02      0.00000000E+00      8.40751311E-04     -5.75686081E-02
      1.80000000E-02      1.80000000E-02      2.07500000E-02      0.00000000E+00      8.43697022E-04     -5.83991579E-02
      1.82500000E-02      1.82500000E-02      2.07500000E-02      0.00000000E+00      8.53975435E-04     -5.92396630E-02
      1.85000000E-02      1.85000000E-02      2.07500000E-02      0.00000000E+00      8.71640860E-04     -6.00821859E-02
      1.87500000E-02      1.87500000E-02      2.07500000E-02      0.00000000E+00      8.96879490E-04     -6.09649132E-02
      1.90000000E-02      1.90000000E-02      2.07500000E-02      0.00000000E+00      9.30229446E-04     -6.18490177E-02
      1.92500000E-02      1.92500000E-02      2.07500000E-02      0.00000000E+00      9.72136104E-04     -6.28051948E-02
      1.95000000E-02      1.95000000E-02      2.07500000E-02      0.00000000E+00      1.02347676E-03     -6.37630765E-02
      1.97500000E-02      1.97500000E-02      2.07500000E-02      0.00000000E+00      1.08504352E-03     -6.48301339E-02
      2.00000000E-02      2.00000000E-02      2.07500000E-02      0.00000000E+00      1.15828738E-03     -6.58991874E-02
      2.02500000E-02      2.02500000E-02      2.07500000E-02      0.00000000E+00      1.24451482E-03     -6.71230676E-02
      2.05000000E-02      2.05000000E-02      2.07500000E-02      0.00000000E+00      1.34607442E-03     -6.83494620E-02
      2.07500000E-02      2.07500000E-02      2.07500000E-02      0.00000000E+00      1.46509381E-03     -6.97902443E-02
      2.10000000E-02      2.10000000E-02      2.07500000E-02      0.00000000E+00      1.60544316E-03     -7.12344494E-02
      2.12500000E-02      2.12500000E-02      2.07500000E-02      0.00000000E+00      1.77051807E-03     -7.29757328E-02
      2.15000000E-02      2.15000000E-02      2.07500000E-02      0.00000000E+00      1.96708793E-03     -7.47214974E-02
      2.17500000E-02      2.17500000E-02      2.07500000E-02      0.00000000E+00      2.20141114E-03     -7.68839170E-02
      2.20000000E-02      2.20000000E-02      2.07500000E-02      0.00000000E+00      2.48572387E-03     -7.90545981E-02
      2.22500000E-02      2.22500000E-02      2.07500000E-02      0.00000000E+00      2.82895782E-03     -8.18371876E-02
      2.25000000E-02      2.25000000E-02      2.07500000E-02      0.00000000E+00      3.25789345E-03     -8.46251869E-02
      2.27500000E-02      2.27500000E-02      2.07500000E-02      0.00000000E+00      3.80321839E-03     -8.83187431E-02
      2.30000000E-02      2.30000000E-02      2.07500000E-02      0.00000000E+00      4.51941731E-03     -9.20588916E-02
      2.32500000E-02      2.32500000E-02      2.07500000E-02      0.00000000E+00      5.40368197E-03     -9.74478265E-02
      2.35000000E-02      2.35000000E-02      2.07500000E-02      0.00000000E+00      6.66114066E-03     -1.02738449E-01
      2.37500000E-02      2.37500000E-02      2.07500000E-02      0.00000000E+00      8.70064683E-03     -1.10363565E-01
      2.40000000E-02      2.40000000E-02      2.07500000E-02      0.00000000E+00      1.20809125E-02     -1.18721657E-01
      2.42500000E-02      2.42500000E-02      2.07500000E-02      0.00000000E+00      1.56847077E-02     -1.36327298E-01
      2.45000000E-02      2.45000000E-02      2.07500000E-02      0.00000000E+00      2.45241649E-02     -1.50754702E-01
      2.47500000E-02      2.47500000E-02      2.07500000E-02      0.00000000E+00      5.12200836E-02     -1.71414549E-01
      2.50000000E-02      2.50000000E-02      2.07500000E-02      0.00000000E+00      1.32739078E-01     -2.09795418E-01
      2.52500000E-02      2.52500000E-02      2.07500000E-02      0.00000000E+00      1.42522722E-01     -1.23437938E-01
      2.55000000E-02      2.55000000E-02      2.07500000E-02      0.00000000E+00      1.39252454E-01     -8.83114204E-02
      2.57500000E-02      2.57500000E-02      2.07500000E-02      0.00000000E+00      1.30710776E-01     -7.41227641E-02
      2.60000000E-02      2.60000000E-02      2.07500000E-02      0.00000000E+00      1.23901045E-01     -6.09819637E-02
      2.62500000E-02      2.62500000E-02      2.07500000E-02      0.00000000E+00      1.19341133E-01     -5.32578967E-02
      2.65000000E-02      2.65000000E-02      2.07500000E-02      0.00000000E+00      1.15656813E-01     -4.59485336E-02
      2.67500000E-02      2.67500000E-02      2.07500000E-02      0.00000000E+00      1.12480002E-01     -4.11831067E-02
      2.70000000E-02      2.70000000E-02      2.07500000E-02      0.00000000E+00      1.09885962E-01     -3.61724079E-02
      2.72500000E-02      2.72500000E-02      2.07500000E-02      0.00000000E+00      1.07749472E-01     -3.22279713E-02
      2.75000000E-02      2.75000000E-02      2.07500000E-02      0.00000000E+00      1.05969180E-01     -2.83036302E-02
      2.77500000E-02      2.77500000E-02      2.07500000E-02      0.00000000E+00      1.04465376E-01     -2.50100048E-02
      2.80000000E-02      2.80000000E-02      2.07500000E-02      0.00000000E+00      1.03208162E-01     -2.16917011E-02
      2.82500000E-02      2.82500000E-02      2.07500000E-02      0.00000000E+00      1.02156898E-01     -1.87395751E-02
      2.85000000E-02      2.85000000E-02      2.07500000E-02      0.00000000E+00      1.01289222E-01     -1.57818357E-02
      2.87500000E-02      2.87500000E-02      2.07500000E-02      0.00000000E+00      1.00581201E-01     -1.30499627E-02
      2.90000000E-02      2.90000000E-02      2.07500000E-02      0.00000000E+00      1.00020998E-01     -1.03124408E-02
      2.92500000E-02      2.92500000E-02      2.07500000E-02      0.00000000E+00      9.95949689E-02     -7.70635857E-03
      2.95000000E-02      2.95000000E-02      2.07500000E-02      0.00000000E+00      9.92969070E-02     -5.09671899E-03
      2.97500000E-02      2.97500000E-02      2.07500000E-02      0.00000000E+00      9.91195530E-02     -2.54917083E-03
      3.00000000E-02      3.00000000E-02      2.07500000E-02      0.00000000E+00      9.90612766E-02      1.37861722E-06
      3.02500000E-02      3.02500000E-02      2.07500000E-02      0.00000000E+00      9.91195530E-02      2.54917083E-03
      3.05000000E-02      3.05000000E-02      2.07500000E-02      0.00000000E+00      9.92969070E-02      5.09972027E-03
      3.07500000E-02      3.07500000E-02      2.07500000E-02      0.00000000E+00      9.95949689E-02      7.70635857E-03
      3.10000000E-02      3.10000000E-02      2.07500000E-02      0.00000000E+00      1.00020998E-01      1.03159981E-02
      3.12500000E-02      3.12500000E-02      2.07500000E-02      0.00000000E+00      1.00581201E-01      1.30499627E-02
      3.15000000E-02      3.15000000E-02      2.07500000E-02      0.00000000E+00      1.01289222E-01      1.57874847E-02
      3.17500000E-02      3.17500000E-02      2.07500000E-02      0.00000000E+00      1.02156898E-01      1.87395751E-02
      3.20000000E-02      3.20000000E-02      2.07500000E-02      0.00000000E+00      1.03208162E-01      2.16973146E-02
      3.22500000E-02      3.22500000E-02      2.07500000E-02      0.00000000E+00      1.04465376E-01      2.50100048E-02
      3.25000000E-02      3.25000000E-02      2.07500000E-02      0.00000000E+00      1.05969180E-01      2.83283086E-02
      3.27500000E-02      3.27500000E-02      2.07500000E-02      0.00000000E+00      1.07749472E-01      3.22279713E-02
      3.30000000E-02      3.30000000E-02      2.07500000E-02      0.00000000E+00      1.09885962E-01      3.61523125E-02
      3.32500000E-02      3.32500000E-02      2.07500000E-02      0.00000000E+00      1.12480002E-01      4.11831067E-02
      3.35000000E-02      3.35000000E-02      2.07500000E-02      0.00000000E+00      1.15656813E-01      4.61938055E-02
      3.37500000E-02      3.37500000E-02      2.07500000E-02      0.00000000E+00      1.19341133E-01      5.32578967E-02
      3.40000000E-02      3.40000000E-02      2.07500000E-02      0.00000000E+00      1.23901045E-01      6.05672598E-02
      3.42500000E-02      3.42500000E-02      2.07500000E-02      0.00000000E+00      1.30710776E-01      7.41227641E-02
      3.45000000E-02      3.45000000E-02      2.07500000E-02      0.00000000E+00      1.39252454E-01      8.72635645E-02
      3.47500000E-02      3.47500000E-02      2.07500000E-02      0.00000000E+00      1.42522722E-01      1.23437938E-01
      3.50000000E-02      3.50000000E-02      2.07500000E-02      0.00000000E+00      1.32739078E-01      1.58564455E-01
      3.52500000E-02      3.52500000E-02      2.07500000E-02      0.00000000E+00      5.12200836E-02      1.71414549E-01
      3.55000000E-02      3.55000000E-02      2.07500000E-02      0.00000000E+00      2.45241649E-02      1.33033680E-01
      3.57500000E-02      3.57500000E-02      2.07500000E-02      0.00000000E+00      1.56847077E-02      1.36327298E-01
      3.60000000E-02      3.60000000E-02      2.07500000E-02      0.00000000E+00      1.20809125E-02      1.21899894E-01
      3.62500000E-02      3.62500000E-02      2.07500000E-02      0.00000000E+00      8.70064683E-03      1.10363565E-01
      3.65000000E-02      3.65000000E-02      2.07500000E-02      0.00000000E+00      6.66114066E-03      1.02005473E-01
      3.67500000E-02      3.67500000E-02      2.07500000E-02      0.00000000E+00      5.40368197E-03      9.74478265E-02
      3.70000000E-02      3.70000000E-02      2.07500000E-02      0.00000000E+00      4.51941731E-03      9.21572041E-02
      3.72500000E-02      3.72500000E-02      2.07500000E-02      0.00000000E+00      3.80321839E-03      8.83187431E-02
      3.75000000E-02      3.75000000E-02      2.07500000E-02      0.00000000E+00      3.25789345E-03      8.45785946E-02
      3.77500000E-02      3.77500000E-02      2.07500000E-02      0.00000000E+00      2.82895782E-03      8.18371876E-02
      3.80000000E-02      3.80000000E-02      2.07500000E-02      0.00000000E+00      2.48572387E-03      7.90491884E-02
      3.82500000E-02      3.82500000E-02      2.07500000E-02      0.00000000E+00      2.20141114E-03      7.68839170E-02
      3.85000000E-02      3.85000000E-02      2.07500000E-02      0.00000000E+00      1.96708793E-03      7.47132359E-02
      3.87500000E-02      3.87500000E-02      2.07500000E-02      0.00000000E+00      1.77051807E-03      7.29757328E-02
      3.90000000E-02      3.90000000E-02      2.07500000E-02      0.00000000E+00      1.60544316E-03      7.12299681E-02
      3.92500000E-02      3.92500000E-02      2.07500000E-02      0.00000000E+00      1.46509381E-03      6.97902443E-02
      3.95000000E-02      3.95000000E-02      2.07500000E-02      0.00000000E+00      1.34607442E-03      6.83460392E-02
      3.97500000E-02      3.97500000E-02      2.07500000E-02      0.00000000E+00      1.24451482E-03      6.71230676E-02
      4.00000000E-02      4.00000000E-02      2.07500000E-02      0.00000000E+00      1.15828738E-03      6.58966731E-02
      4.02500000E-02      4.02500000E-02      2.07500000E-02      0.00000000E+00      1.08504352E-03      6.48301339E-02
      4.05000000E-02      4.05000000E-02      2.07500000E-02      0.00000000E+00      1.02347676E-03      6.37610805E-02
      4.07500000E-02      4.07500000E-02      2.07500000E-02      0.00000000E+00      9.72136104E-04      6.28051948E-02
      4.10000000E-02      4.10000000E-02      2.07500000E-02      0.00000000E+00      9.30229446E-04      6.18473131E-02
      4.12500000E-02      4.12500000E-02      2.07500000E-02      0.00000000E+00      8.96879490E-04      6.09649132E-02
      4.15000000E-02      4.15000000E-02      2.07500000E-02      0.00000000E+00      8.71640860E-04      6.00808086E-02
      4.17500000E-02      4.17500000E-02      2.07500000E-02      0.00000000E+00      8.53975435E-04      5.92396630E-02
      4.20000000E-02      4.20000000E-02      2.07500000E-02      0.00000000E+00      8.43697022E-04      5.83971401E-02
      4.22500000E-02      4.22500000E-02      2.07500000E-02      0.00000000E+00      8.40751311E-04      5.75686081E-02
      4.25000000E-02      4.25000000E-02      2.07500000E-02      0.00000000E+00      8.44821365E-04      5.67380583E-02
      4.27500000E-02      4.27500000E-02      2.07500000E-02      0.00000000E+00      8.54782289E-04      5.58905007E-02
      4.30000000E-02      4.30000000E-02      2.07500000E-02      0.00000000E+00      8.70294805E-04      5.50443643E-02
      4.32500000E-02      4.32500000E-02      2.07500000E-02      0.00000000E+00      8.96096470E-04      5.41697913E-02
      4.35000000E-02      4.35000000E-02      2.07500000E-02      0.00000000E+00      9.23470688E-04      5.32800858E-02
      4.37500000E-02      4.37500000E-02      2.07500000E-02      0.00000000E+00      9.19597373E-04      5.22852654E-02
      4.40000000E-02      4.40000000E-02      2.07500000E-02      0.00000000E+00      8.62453757E-04      5.13577088E-02
      4.42500000E-02      4.42500000E-02      2.07500000E-02      0.00000000E+00      8.97254630E-04      5.07412459E-02
      4.45000000E-02      4.45000000E-02      2.07500000E-02      0.00000000E+00      6.43694458E-04      4.97790763E-02
      4.47500000E-02      4.47500000E-02      2.07500000E-02      0.00000000E+00     -1.13320015E-03      4.80394138E-02
      4.50000000E-02      4.50000000E-02      2.07500000E-02      0.00000000E+00     -7.45890295E-03      4.81101463E-02
      4.52500000E-02      4.52500000E-02      2.07500000E-02      0.00000000E+00     -6.98463099E-04      3.64060929E-02
      4.55000000E-02      4.55000000E-02      2.07500000E-02      0.00000000E+00      4.74745652E-03      3.34585295E-02
      4.57500000E-02      4.57500000E-02      2.07500000E-02      0.00000000E+00      9.90727797E-03      3.20913104E-02
      4.60000000E-02      4.60000000E-02      2.07500000E-02      0.00000000E+00      1.40243040E-02      3.20619015E-02
      4.62500000E-02      4.62500000E-02      2.07500000E-02      0.00000000E+00      1.74135326E-02      3.28405986E-02
      4.65000000E-02      4.65000000E-02      2.07500000E-02      0.00000000E+00      2.03684903E-02      3.35187456E-02
      4.67500000E-02      4.67500000E-02      2.07500000E-02      0.00000000E+00      2.30591024E-02      3.44095826E-02
      4.70000000E-02      4.70000000E-02      2.07500000E-02      0.00000000E+00      2.55154041E-02      3.52942922E-02
      4.72500000E-02      4.72500000E-02      2.07500000E-02      0.00000000E+00      2.78047880E-02      3.62620698E-02
      4.75000000E-02      4.75000000E-02      2.07500000E-02      0.00000000E+00      2.99655774E-02      3.72300489E-02
      4.77500000E-02      4.77500000E-02      2.07500000E-02      0.00000000E+00      3.20358211E-02      3.82467718E-02
      4.80000000E-02      4.80000000E-02      2.07500000E-02      0.00000000E+00      3.40355663E-02      3.92622703E-02
      4.82500000E-02      4.82500000E-02      2.07500000E-02      0.00000000E+00      3.59877128E-02      4.03174911E-02
      4.85000000E-02      4.85000000E-02      2.07500000E-02      0.00000000E+00      3.79076974E-02      4.13729687E-02
      4.87500000E-02      4.87500000E-02      2.07500000E-02      0.00000000E+00      3.98117354E-02      4.24722608E-02
      4.90000000E-02      4.90000000E-02      2.07500000E-02      0.00000000E+00      4.17123963E-02      4.35719353E-02
      4.92500000E-02      4.92500000E-02      2.07500000E-02      0.00000000E+00      4.36226718E-02      4.47261085E-02
      4.95000000E-02      4.95000000E-02      2.07500000E-02      0.00000000E+00      4.55544261E-02      4.58809820E-02
      4.97500000E-02      4.97500000E-02      2.07500000E-02      0.00000000E+00      4.75195536E-02      4.71070547E-02
      5.00000000E-02      5.00000000E-02      2.07500000E-02      0.00000000E+00      4.95306383E-02      4.83341115E-02
      5.02500000E-02      5.02500000E-02      2.07500000E-02      0.00000000E+00      5.15999995E-02      4.96562255E-02
      5.05000000E-02      5.05000000E-02      2.07500000E-02      0.00000000E+00      5.37424809E-02      5.09797181E-02
      5.07500000E-02      5.07500000E-02      2.07500000E-02      0.00000000E+00      5.59722971E-02      5.24324671E-02
      5.10000000E-02      5.10000000E-02      2.07500000E-02      0.00000000E+00      5.83086016E-02      5.38872331E-02
      5.12500000E-02      5.12500000E-02      2.07500000E-02      0.00000000E+00      6.07694655E-02      5.55222235E-02
      5.15000000E-02      5.15000000E-02      2.07500000E-02      0.00000000E+00      6.33815241E-02      5.71600709E-02
      5.17500000E-02      5.17500000E-02      2.07500000E-02      0.00000000E+00      6.61691853E-02      5.90570663E-02
      5.20000000E-02      5.20000000E-02      2.07500000E-02      0.00000000E+00      6.91725956E-02      6.09593282E-02
      5.22500000E-02      5.22500000E-02      2.07500000E-02      0.00000000E+00      7.24295550E-02      6.32567766E-02
      5.25000000E-02      5.25000000E-02      2.07500000E-02      0.00000000E+00      7.60036031E-02      6.55595819E-02
      5.27500000E-02      5.27500000E-02      2.07500000E-02      0.00000000E+00      7.99446839E-02      6.84900412E-02
      5.30000000E-02      5.30000000E-02      2.07500000E-02      0.00000000E+00      8.43705579E-02      7.14464033E-02
      5.32500000E-02      5.32500000E-02      2.07500000E-02      0.00000000E+00      8.94208995E-02      7.55926893E-02
      5.35000000E-02      5.35000000E-02      2.07500000E-02      0.00000000E+00      9.52883226E-02      7.97186786E-02
      5.37500000E-02      5.37500000E-02      2.07500000E-02      0.00000000E+00      1.01927601E-01      8.60984953E-02
      5.40000000E-02      5.40000000E-02      2.07500000E-02      0.00000000E+00      1.09876741E-01      9.27297114E-02
      5.42500000E-02      5.42500000E-02      2.07500000E-02      0.00000000E+00      1.20803498E-01      1.06214075E-01
      5.45000000E-02      5.45000000E-02      2.07500000E-02      0.00000000E+00      1.34504997E-01      1.19384284E-01
      5.47500000E-02      5.47500000E-02      2.07500000E-02      0.00000000E+00      1.43221185E-01      1.59844031E-01
      5.50000000E-02      5.50000000E-02      2.07500000E-02      0.00000000E+00      1.40197981E-01      1.97918111E-01
      5.52500000E-02      5.52500000E-02      2.07500000E-02      0.00000000E+00      5.23532837E-02      2.19453963E-01
      5.55000000E-02      5.55000000E-02      2.07500000E-02      0.00000000E+00      2.38804704E-02      1.81002361E-01
      5.57500000E-02      5.57500000E-02      2.07500000E-02      0.00000000E+00      1.47874531E-02      1.87068544E-01
      5.60000000E-02      5.60000000E-02      2.07500000E-02      0.00000000E+00      1.12184587E-02      1.73603310E-01
      5.62500000E-02      5.62500000E-02      2.07500000E-02      0.00000000E+00      7.78104945E-03      1.62648830E-01
      5.65000000E-02      5.65000000E-02      2.07500000E-02      0.00000000E+00      5.73766997E-03      1.55218295E-01
      5.67500000E-02      5.67500000E-02      2.07500000E-02      0.00000000E+00      4.50758550E-03      1.51617618E-01
      5.70000000E-02      5.70000000E-02      2.07500000E-02      0.00000000E+00      3.64912251E-03      1.47216701E-01
      5.72500000E-02      5.72500000E-02      2.07500000E-02      0.00000000E+00      2.94843610E-03      1.44209244E-01
      5.75000000E-02      5.75000000E-02      2.07500000E-02      0.00000000E+00      2.41307209E-03      1.41363245E-01
      5.77500000E-02      5.77500000E-02      2.07500000E-02      0.00000000E+00      1.98820651E-03      1.39405796E-01
      5.80000000E-02      5.80000000E-02      2.07500000E-02      0.00000000E+00      1.64202685E-03      1.37448346E-01
      5.82500000E-02      5.82500000E-02      2.07500000E-02      0.00000000E+00      1.34743570E-03      1.36123580E-01
      5.85000000E-02      5.85000000E-02      2.07500000E-02      0.00000000E+00      1.09544707E-03      1.34795422E-01
      5.87500000E-02      5.87500000E-02      2.07500000E-02      0.00000000E+00      8.73638577E-04      1.33940646E-01
      5.90000000E-02      5.90000000E-02      2.07500000E-02      0.00000000E+00      6.75213719E-04      1.33078986E-01
      5.92500000E-02      5.92500000E-02      2.07500000E-02      0.00000000E+00      4.92957702E-04      1.32595439E-01
      5.95000000E-02      5.95000000E-02      2.07500000E-02      0.00000000E+00      3.22597660E-04      1.32109116E-01
      5.97500000E-02      5.97500000E-02      2.07500000E-02      0.00000000E+00      1.59471291E-04      1.31953202E-01
      6.00000000E-02      6.00000000E-02      2.07500000E-02      0.00000000E+00     -5.40499075E-17      1.31795861E-01];
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