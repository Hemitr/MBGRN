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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Constitution of matrix [R] (equation (10) of the contribution) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
M(nnum(r01,s01),nnum(r01,s01)+m-1) =  Rw21;
M(nnum(r01,s01),nnum(r01,s01)+1)   = -Rw21;
M(nnum(r01,s01),nnum(r01,s01)+m)   = -Rv21;
M(nnum(r01,s01),nnum(r01,s01))     = 2*(Rw21+Rv21);
for s01 = 2:m-1
    nnum(r01,s01)=(r01-1)*m+s01;
    M(nnum(r01,s01),nnum(r01,s01)-1) = -Rw21;
    M(nnum(r01,s01),nnum(r01,s01)+1) = -Rw21;
    M(nnum(r01,s01),nnum(r01,s01)+m) = -Rv21;
    M(nnum(r01,s01),nnum(r01,s01))   = 2*(Rw21+Rv21);
end
s01 = m;
nnum(r01,s01)=(r01-1)*m+s01;
M(nnum(r01,s01),nnum(r01,s01)-1)   = -Rw21;
M(nnum(r01,s01),nnum(r01,s01)+m)   = -Rv21;
M(nnum(r01,s01),nnum(r01,s01)-m+1) =  Rw21;
M(nnum(r01,s01),nnum(r01,s01))     = 2*(Rw21+Rv21);

% First column
s01 = 1;
for r01 = 2:p0-1
    nnum(r01,s01)=(r01-1)*m+s01;
    M(nnum(r01,s01),nnum(r01,s01)+m-1) =  Rw21;
    M(nnum(r01,s01),nnum(r01,s01)+1)   = -Rw21;
    M(nnum(r01,s01),nnum(r01,s01)+m)   = -Rv21;
    M(nnum(r01,s01),nnum(r01,s01)-m)   = -Rv21;
    M(nnum(r01,s01),nnum(r01,s01))     = 2*(Rw21+Rv21);
end

% Last column
s01 = m;
for r01 = 2:p0-1
    nnum(r01,s01)=(r01-1)*m+s01;
    M(nnum(r01,s01),nnum(r01,s01)-m+1) =  Rw21;
    M(nnum(r01,s01),nnum(r01,s01)-1)   = -Rw21;
    M(nnum(r01,s01),nnum(r01,s01)+m)   = -Rv21;
    M(nnum(r01,s01),nnum(r01,s01)-m)   = -Rv21;
    M(nnum(r01,s01),nnum(r01,s01))     = 2*(Rw21+Rv21);
end

% The rest of elements
for r01 = 2:p0-1
    for s01 = 2:m-1
        nnum(r01,s01)=(r01-1)*m+s01;
        M(nnum(r01,s01),nnum(r01,s01)-1) = -Rw21;
        M(nnum(r01,s01),nnum(r01,s01)+1) = -Rw21;
        M(nnum(r01,s01),nnum(r01,s01)+m) = -Rv21;
        M(nnum(r01,s01),nnum(r01,s01)-m) = -Rv21;
        M(nnum(r01,s01),nnum(r01,s01))   = 2*(Rw21+Rv21);
    end
end

% Mobile iron / Magnetic air-gap boundary
% Last line of mobile iron having a link with the Magnetic air-gap
r01 = p0;
s01 = 1;
nnum(r01,s01)=(r01-1)*m+s01;
M(nnum(r01,s01),nnum(r01,s01)-m)   = -Rv21;
M(nnum(r01,s01),nnum(r01,s01)+m-1) =  (Rw21/2+Rw1/2);
M(nnum(r01,s01),nnum(r01,s01)+1)   = -(Rw21/2+Rw1/2);
M(nnum(r01,s01),nnum(r01,s01)+m)   = -Rv1;
M(nnum(r01,s01),nnum(r01,s01))     =  Rv21+Rv1+Rw21+Rw1;
for s01 = 2:m-1
    nnum(r01,s01)=(r01-1)*m+s01;
    M(nnum(r01,s01),nnum(r01,s01)-m) = -Rv21;
    M(nnum(r01,s01),nnum(r01,s01)-1) = -(Rw21/2+Rw1/2);
    M(nnum(r01,s01),nnum(r01,s01)+1) = -(Rw21/2+Rw1/2);
    M(nnum(r01,s01),nnum(r01,s01)+m) = -Rv1;
    M(nnum(r01,s01),nnum(r01,s01))   =  Rv21+Rv1+Rw21+Rw1;
end
s01 = m;
nnum(r01,s01)=(r01-1)*m+s01;
M(nnum(r01,s01),nnum(r01,s01)-m)   = -Rv21;
M(nnum(r01,s01),nnum(r01,s01)-1)   = -(Rw21/2+Rw1/2);
M(nnum(r01,s01),nnum(r01,s01)+m)   = -Rv1;
M(nnum(r01,s01),nnum(r01,s01)-m+1) =  (Rw21/2+Rw1/2);
M(nnum(r01,s01),nnum(r01,s01))     =  Rv21+Rv1+Rw21+Rw1;

% Magnetic air-gap
% First line of the Magnetic air-gap having a link with the mobile iron
r11 = p0+1;
s11 = 1;
nnum(r11,s11)=(r11-1)*m+s11;
M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv1;
M(nnum(r11,s11),nnum(r11,s11)+m-1) =  Rw1;
M(nnum(r11,s11),nnum(r11,s11)+1)   = -Rw1;
M(nnum(r11,s11),nnum(r11,s11)+m)   = -Rv1;
M(nnum(r11,s11),nnum(r11,s11))     = 2*(Rw1+Rv1);
for s11 = 2:m-1
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
    M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw1;
    M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw1;
    M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv1;
    M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rw1+Rv1);
end
s11 = m;
nnum(r11,s11)=(r11-1)*m+s11;
M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv1;
M(nnum(r11,s11),nnum(r11,s11)-1)   = -Rw1;
M(nnum(r11,s11),nnum(r11,s11)+m)   = -Rv1;
M(nnum(r11,s11),nnum(r11,s11)-m+1) =  Rw1;
M(nnum(r11,s11),nnum(r11,s11))      = 2*(Rw1+Rv1);

% First colum
s11 = 1;
for r11 = p0+2:p0+p1-1
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)+m-1) =  Rw1;
    M(nnum(r11,s11),nnum(r11,s11)+1)   = -Rw1;
    M(nnum(r11,s11),nnum(r11,s11)+m)   = -Rv1;
    M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv1;
    M(nnum(r11,s11),nnum(r11,s11))     = 2*(Rw1+Rv1);
end

% Last colum
s11 = m;
for r11 = p0+2:p0+p1-1
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-m+1) =  Rw1;
    M(nnum(r11,s11),nnum(r11,s11)-1)   = -Rw1;
    M(nnum(r11,s11),nnum(r11,s11)+m)   = -Rv1;
    M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv1;
    M(nnum(r11,s11),nnum(r11,s11))     = 2*(Rw1+Rv1);
end

% The rest of elements in magnetic air-gap region
for r11 = p0+2:p0+p1-1
    for s11 = 2:m-1
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw1;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw1;
        M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv1;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
        M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rw1+Rv1);
    end
end

% Magnetic air-gap / Stator boundary
% Last line of magnetic air-gap havin a ling with the first line of the stator armature
r11 = p0+p1;
s11 = 1;
nnum(r11,s11)=(r11-1)*m+s11;
M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv1;
M(nnum(r11,s11),nnum(r11,s11)+m-1) =  (Rw1/2+Rw20/2);
M(nnum(r11,s11),nnum(r11,s11)+1)   = -(Rw1/2+Rw20/2);
M(nnum(r11,s11),nnum(r11,s11)+m)   = -Rv20;
M(nnum(r11,s11),nnum(r11,s11))      = Rv1+Rv20+Rw1+Rw20;
for s11 = 2:m0s-1
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
    M(nnum(r11,s11),nnum(r11,s11)-1) = -(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11)+1) = -(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))    = Rv1+Rv20+Rw1+Rw20;
end
for k = 1:3
    for s11 = m0s+(k-1)*4*m0s
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw1;
        M(nnum(r11,s11),nnum(r11,s11)+m) = -(Rv1/2+Rv20/2);
        M(nnum(r11,s11),nnum(r11,s11))   = 3*Rv1/2+3*Rw1/2+Rw20/2+Rv20/2;
    end
    for s11 = 3*m0s+(k-1)*4*m0s
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw1;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -(Rv1/2+Rv20/2);
        M(nnum(r11,s11),nnum(r11,s11))    = 3*Rv1/2+3*Rw1/2+Rw20/2+Rv20/2;    
    end
end
for dent = 1:2
    for s11 = 3*m0s+1+(dent-1)*4*m0s:5*m0s-1+(dent-1)*4*m0s
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11)+1) = -(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11))    = Rv1+Rv20+Rw1+Rw20;
    end
end
for enc = 1:3
    for s11 = m0s+1+(enc-1)*4*m0s:3*m0s-1+(enc-1)*4*m0s
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw1;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw1;
        M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv1;
        M(nnum(r11,s11),nnum(r11,s11))    = 2*(Rv1+Rw1);
    end
end
for s11 = 11*m0s+1:m-1
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
    M(nnum(r11,s11),nnum(r11,s11)-1) = -(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11)+1) = -(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))   =  Rv1+Rv20+Rw1+Rw20;
end
s11 = m;
nnum(r11,s11)=(r11-1)*m+s11;
M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv1;
M(nnum(r11,s11),nnum(r11,s11)-1)   = -(Rw1/2+Rw20/2);
M(nnum(r11,s11),nnum(r11,s11)+m)   = -Rv20;
M(nnum(r11,s11),nnum(r11,s11)-m+1) =  (Rw1/2+Rw20/2);
M(nnum(r11,s11),nnum(r11,s11))     =  Rv1+Rv20+Rw1+Rw20;

% Stator armature
% First column
s11 = 1;
for r11 = p0+p1+1:p0+p1+p-2;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)+m-1) =  Rw20;
    M(nnum(r11,s11),nnum(r11,s11)+1)   = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11)+m)   = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))     = 2*(Rv20+Rw20);
end
r11 = p0+p1+p-1;
nnum(r11,s11)=(r11-1)*m+s11;
M(nnum(r11,s11),nnum(r11,s11)+1)   = -Rw20;
M(nnum(r11,s11),nnum(r11,s11)+m-1) =  Rw20;
M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv20;
M(nnum(r11,s11),nnum(r11,s11))    = 2*(Rv20+Rw20);

% Last column
s11 = m;
for r11 = p0+p1+1:p0+p1+p-2;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-1)   = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)-m+1) =  Rw20;
    M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11)+m)   = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))     = 2*(Rv20+Rw20);
end
r11 = p0+p1+p-1;
nnum(r11,s11)=(r11-1)*m+s11;
M(nnum(r11,s11),nnum(r11,s11)-1)   = -Rw20;
M(nnum(r11,s11),nnum(r11,s11)-m+1) =  Rw20;
M(nnum(r11,s11),nnum(r11,s11)-m)   = -Rv20;
M(nnum(r11,s11),nnum(r11,s11))    = 2*(Rv20+Rw20);

% Columns s11 = m0s, 5*m0s, 9*m0s
for dent = 1:3
    s11 = m0s+(dent-1)*4*m0s;
    for r11 = p0+p1+1:p0+p1+2*p0s-1
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw1;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -(Rv1/2+Rv20/2);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -(Rv1/2+Rv20/2);
        M(nnum(r11,s11),nnum(r11,s11))   =  Rw20+Rw1+Rv1+Rv20;
    end
    r11 = p0+p1+2*p0s;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)+1) = -(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11)-m) = -(Rv1/2+Rv20/2);
    M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))   =  3*Rw20/2+3*Rv20/2+Rv1/2+Rw1/2;
    for r11 = p0+p1+2*p0s+1:p0+p1+p-2
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
    end
    r11 = p0+p1+p-1;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
end

% Columns s11 = 3*m0s, 7*m0s, 11*m0s
for dent = 1:3
    s11 = 3*m0s+(dent-1)*4*m0s;
    for r11 = p0+p1+1:p0+p1+2*p0s-1
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw1;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -(Rv1/2+Rv20/2);
        M(nnum(r11,s11),nnum(r11,s11)+m) = -(Rv1/2+Rv20/2);
        M(nnum(r11,s11),nnum(r11,s11))   =  Rw1+Rw20+Rv1+Rv20;
    end
    r11 = p0+p1+2*p0s;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-1) = -(Rw1/2+Rw20/2);
    M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)-m) = -(Rv1/2+Rv20/2);
    M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))   =  3*Rw20/2+3*Rv20/2+Rv1/2+Rw1/2;
    for r11 = p0+p1+2*p0s+1:p0+p1+p-2
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
    end
    r11 = p0+p1+p-1;
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
end

% First and last half teeth excluding s11 = 1, m0s and 12*m0s
for dent = 1:2
    for s11 = 2+(dent-1)*11*m0s:m0s-1+(dent-1)*11*m0s
        for r11=p0+p1+1:p0+p1+p-2
            nnum(r11,s11)=(r11-1)*m+s11;
            M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
            M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
            M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
            M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
            M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
        end
        r11 = p0+p1+p-1;
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
    end
end
s11 = 11*m0s+1;
for r11=p0+p1+1:p0+p1+p-2
    nnum(r11,s11)=(r11-1)*m+s11;
    M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
    M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
    M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
end
r11 = p0+p1+p-1;
nnum(r11,s11)=(r11-1)*m+s11;
M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);

% Other teeth excluding s11 = 5*m0s and 9*m0s
for dent = 1:2
    for s11 = 3*m0s+1+(dent-1)*4*m0s:5*m0s-1+(dent-1)*4*m0s
        for r11=p0+p1+1:p0+p1+p-2
            nnum(r11,s11)=(r11-1)*m+s11;
            M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
            M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
            M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
            M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
            M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
        end
        r11 = p0+p1+p-1;
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
    end
end

% Slots region and back-iron regions on top of the slots
for enc = 1:3
    for s11 = m0s+1+(enc-1)*4*m0s:3*m0s-1+(enc-1)*4*m0s
        for r11 = p0+p1+1:p0+p1+2*p0s-1
            nnum(r11,s11)=(r11-1)*m+s11;
            M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw1;
            M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw1;
            M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
            M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv1;
            M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv1+Rw1);
        end
        r11 = p0+p1+2*p0s;
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11)+1) = -(Rw1/2+Rw20/2);
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv1;
        M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11))    = Rv1+Rv20+Rw1+Rw20;
        for r11 = p0+p1+2*p0s+1:p0+p1+p-2
            nnum(r11,s11)=(r11-1)*m+s11;
            M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
            M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
            M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
            M(nnum(r11,s11),nnum(r11,s11)+m) = -Rv20;
            M(nnum(r11,s11),nnum(r11,s11))    = 2*(Rv20+Rw20);
        end
        r11 = p0+p1+p-1;
        nnum(r11,s11)=(r11-1)*m+s11;
        M(nnum(r11,s11),nnum(r11,s11)-1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)+1) = -Rw20;
        M(nnum(r11,s11),nnum(r11,s11)-m) = -Rv20;
        M(nnum(r11,s11),nnum(r11,s11))   = 2*(Rv20+Rw20);
    end
end

toc
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Excitation vector [F] %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FMMPM = (Br*h/mu0);
s = m0m;
E0(1,s) = FMMPM/2;
s = m-m0m;
E0(1,s) = -FMMPM/2;
s = m0m+m;
E0(1,s) = -FMMPM/2;
s = 2*m-m0m;
E0(1,s) = FMMPM/2;
E0(1,240) = 0;
for i=2:p0
    s = m0m;
    E0(i,s) = FMMPM;
    s = m-m0m;
    E0(i,s) = -FMMPM;
    s = m0m+m;
    E0(i,s) = -FMMPM;
    s = 2*m-m0m;
    E0(i,s) = FMMPM;
    E0(i,240) = 0;
end
s = m0m;
E0(p0+1,s) = FMMPM/2;
s = m-m0m;
E0(p0+1,s) = -FMMPM/2;
s = m0m+m;
E0(p0+1,s) = -FMMPM/2;
s = 2*m-m0m;
E0(p0+1,s) = FMMPM/2;
E0(p0+1,240) = 0;

toc
tic

% Displacement vector
xd   = [0:l:tp/2];

for u=1:length(xd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Excitation vector [F] adaptation with xd %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:m
    for r = p0:2*p0
        k1    = (r-1)*m+s;
        E(k1) = E0(r-p0+1,s+u-1);
    end    
end
E(nn-m) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Solution vector  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

C=M\E';

% Computation of the flux in a phase coil
nnum01=(p0+p1+20)*m+1;
Phip(1) = -(C(nnum01)+C(nnum01+m-1));
for i=2:6*m0s
    nnum01=(p0+p1+20)*m+i;
    Phip(i)=-(C(nnum01)-C(nnum01-1));
end
for i=1:6*m0s
    nnum11=(p0+p1+20)*m+i+6*m0s;
    Phim(i)=-(C(nnum11)-C(nnum11-1));
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
