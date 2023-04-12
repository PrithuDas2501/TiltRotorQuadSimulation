% Breaking any general NACA Airfoil into discrete panels and plotting it
NACA = 0015;

T = rem(NACA,100);
P = rem((NACA-T)/100,10);
M = ((NACA-T)/100 - P)/10;

P = 0.1*P;
M = 0.01*M;
T = 0.01*T;


N = 11; % (Number of Panels = 2N-2)
x = linspace(0, 1, N);
yu = [];
yl = [];
yc = [];
yt = [];
for i = 1:N
    y = func(x(i),P,M,T);
    yu(i) = y(1);
    yl(i) = y(2);
    yc(i) = y(3);
    yt(i) = y(4);
end
yl(N) = [];
yl(N) = yu(N);
x2 = flip(x);
x2(N) = [];
x2(1) = [];
y2 = flip(yu);
y2(N) = [];
y2(1) = [];
panel_end_x = cat(2,x,x2);  % --------------------
panel_end_y = cat(2,yl,y2); % --------------------


% Panel Inclinations
upanel_inclination = [];
lpanel_inclination = [];
for i = 1:N-1
    upanel_inclination(i) = atan((yu(i+1)-yu(i))/(x(i+1) - x(i)));
    lpanel_inclination(i) = atan((yl(i+1)-yl(i))/(x(i+1) - x(i)));
end
rev_upanel_inclination = flip(upanel_inclination);
panel_inclination = cat(2,lpanel_inclination,rev_upanel_inclination);


% Panel Control Point Coordinates
CPCoord = zeros(2,2*N-2);
X = [];
Y1 = [];
Y2 = [];
for i = 1:N-1
    X(i)  = (x(i) + x(i+1))/2;
    Y1(i) = (yl(i) + yl(i+1))/2;
    Y2(i) = (yu(N+1-i) + yu(N-i))/2;
end
X = cat(2,X,flip(X));
Y = cat(2,Y1,Y2);
for i = 1:2*N-2
    CPCoord(1,i) = X(i);
    CPCoord(2,i) = Y(i);
end


%Panel Length (s)
s = zeros(1,2*N-2);
for i = 1:2*N-3
    s(1,i) = ( (panel_end_x(i+1) - panel_end_x(i))^2 + (panel_end_y(i+1) - panel_end_y(i))^2 )^(0.5);
end
s(1,2*N-2) = ( (panel_end_x(1) - panel_end_x(2*N-2))^2 + (panel_end_y(1) - panel_end_y(2*N-2))^2 )^(0.5);


% The Distance Matrix
distance_matrix = zeros(2*N-2, 2*N-2);
for i = 1:2*N-2
    for j = i:2*N-2
        if i==j
            r=0;
        else
            r = ( (CPCoord(1,i)-CPCoord(1,j))^2 + (CPCoord(2,i)-CPCoord(2,j))^2 )^(0.5);
        end
        distance_matrix(i,j) = r;
        distance_matrix(j,i) = r;
    end
end


% The Alpha Matrix
alpha_matrix = zeros(2*N-2, 2*N-2);
for i = 1:2*N-2
    for j = 1:2*N-2
        alpha_matrix(i,j) = panel_inclination(i) - panel_inclination(j);
    end
end
alpha_matrix;

% Finding Theta1 and Theta2 for every pairing of panel
theta_matrix = zeros(2*N-2, 2*N-2, 2);
for i = 1:2*N-2
    for j = 1:2*N-2
        Normal_Pointij = transpose([tan(panel_inclination(j)), -1; -cot(panel_inclination(j)), -1]\[tan(panel_inclination(j))*CPCoord(1,j) - CPCoord(2,j);-CPCoord(1,i)*cot(panel_inclination(j)) - CPCoord(2,i)]);
        if j==2*N-2
            theta = t(Normal_Pointij, [panel_end_x(j), panel_end_y(j)], [panel_end_x(1), panel_end_y(1)], [CPCoord(1,i), CPCoord(2,i)]);
        else
            theta = t(Normal_Pointij, [panel_end_x(j), panel_end_y(j)], [panel_end_x(j+1), panel_end_y(j+1)], [CPCoord(1,i), CPCoord(2,i)]);
        end
        theta_matrix(i,j,1) = theta(1);
        theta_matrix(i,j,2) = theta(2);
    end
end
theta_matrix;

% Rnot matrix
Rnot = zeros(2*N-2, 2*N-2);
for i = 1:2*N-2
    for j = 1:2*N-2
        Normal_Pointij = transpose([tan(panel_inclination(j)), -1; -cot(panel_inclination(j)), -1]\[tan(panel_inclination(j))*CPCoord(1,j) - CPCoord(2,j);-CPCoord(1,i)*cot(panel_inclination(j)) - CPCoord(2,i)]);
        r = LeastDistance(Normal_Pointij, [CPCoord(1,i), CPCoord(2,i)]);
        Rnot(i,j) = r;
    end
end
Rnot;


% We are assuming each panel to have a uniform Gamma.
% Constructing the Matrix that calculates the normal velocity induced.
A = zeros(2*N-2, 2*N-2);
for i = 1:2*N-2
    for j = 1:2*N-2
        if i==j
            A(i,j) = 0;
        else
            A(i,j) = (s(j)/(4*pi*Rnot(i,j)))*(cos(theta_matrix(i,j,1) - alpha_matrix(i,j)) - cos(theta_matrix(i,j,2) - alpha_matrix(i,j))); % FIGURE OUT WHY THIS IS GIVING OUT OF BOUNDS INDEX ERRORRRRRR!!!!!!!!!!!!!!!!!!!
        end
    end
end
A;

B = zeros(2*N-2,1);
for i = 1:2*N-2
    B(i,1) = cos(alpha_matrix(i));
end
B;

%Think of a way to ensure KUTTA CONDITION
Solution = A\B










plot(x,yu)
hold
plot(x,yl)
plot(x,yc)
axis equal




function y = func(x,p,m,t)
if x <p
    yc = (m/(p^2))*(2*p*x - x^2);
    slope = (m/(p^2))*(2*p - 2*x);
end

if x >=p;
    yc = (m/(1-p)^2)*((1-2*p) + 2*p*x - x^2);
    slope = (m/(1-p)^2)*(2*p - 2*x);   
end

yt = (t/0.2)*(0.2969*x^(0.5) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4);

theta = atan(slope);

yu = yc + yt*cos(theta);
yl = yc - yt*cos(theta);

y = [yu,yl,yc,yt];

end


function theta = t(NP, p1, p2, p)
R1 = ((p(1)-p1(1))^2+(p(2)-p1(2))^2)^(0.5);
R2 = ((p(1)-p2(1))^2+(p(2)-p2(2))^2)^(0.5);
R  = ((p(1)-NP(1))^2+(p(2)-NP(2))^2)^(0.5);
r1 = ((NP(1)-p1(1))^2+(NP(2)-p1(2))^2)^(0.5);
r2 = ((NP(1)-p2(1))^2+(NP(2)-p2(2))^2)^(0.5);
if (NP(1)-p1(1))*(NP(1)-p2(1))<=0
    t1 = -acos(R/R1);
    t2 =  acos(R/R2);
elseif (NP(1)-p1(1))*(NP(1)-p2(1))>0
    if r1<r2
        t1 =  acos(R/R1);
        t2 =  acos(R/R2);
    else
        t1 = -acos(R/R1);
        t2 = -acos(R/R2);
    end
end
theta = [t1,t2];
end

function r = LeastDistance(NP, p)
r = ((p(1)-NP(1))^2+(p(2)-NP(2))^2)^(0.5);
end



