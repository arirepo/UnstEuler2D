function [A B f1] = S(colength,z,NBSCHTYPE,INTSCHTYPE)
%This function returns an integration operator for the column operand of
%length "colength" and depth "z", where "z" is the number of equations when 
%applied to a system of equations. That means this function constructs  
% matrices [[A]], [[B]], and column vector [f1] given in equation: 
%            [[A]] [f_tilde] = dx*([f1]*f(x1) + [[B]]*[f]) 
% when solved for unknown vector [f_tilde], it gives us the numerical
% integral of discretetized function [f] whose integral is ZERO at the
% start of the interval which is x1 and the value of the function at this
% point is f(x1) as given above.
%NBSCHTYPE determines the type of 
%near-boundary schemes. For example [4 4 6] says 
%that scheme at point 1 must be of type 4 and at point 2 must be of type 4 and at
%point 3 must be of type 6. Refer to the following database to see the
%nomenclature for type and other characteristics of schemes.
%For the interor scheme, INTSCHTYPE is a number that determines the type of
%interior scheme. Refer to the following database to see the
%nomenclature for type and other characteristics of interior schemes.
% ----------------- EXAMPLES ---------------------------------------------
%S(n,1,[3 3],2) : the arguments are:
%   n = length of the column
%   1 = number of variables in system of equations
%   [3 3] number and order of boundary scheme
%       length([3 3]) means number of near boundary points
%       3 means fourth-order
%   2 = type of interior scheme: 2 means fourth-order
%example: S(n,1,[1],1) return same thing with seond-order interior and near
%boundary schemes. Interestingly it is trapezodial or Crank - Nicholson
%scheme!!!

colength = colength + 1;

%COEFFICIENTS DATABASE (AQUIRED USING MAPLE CODE)
% Allocating memory for interior schemes.
Thetai = zeros(4);
 Lambdai = zeros(4,4); 
% Allocating memory for near boundary schemes.
Thetanbs = zeros(4,8);
 Lambdanbs = zeros(4,8,9); 
% The Interior Scheme.
% NOMENCLATURE------------------
%Thetai(TYPE) = returns the coefficient Theta of the interior scheme given
%the TYPE of scheme.
%Lambdai(TYPE, j) = returns the jth coefficient Lambda of the interior 
%scheme given the TYPE of scheme.
%NOTE : STARTING FROM SECOND-0RDER ACCURACY
Thetai(1)=-1/2;
Thetai(2)=-1/2;
Thetai(3)=-1/2;
Thetai(4)=-1/2;
Lambdai(1,1)=-1/4;
Lambdai(1,2)=0;
Lambdai(1,3)=0;
Lambdai(1,4)=0;
Lambdai(2,1)=-7/24;
Lambdai(2,2)=1/48;
Lambdai(2,3)=0;
Lambdai(2,4)=0;
Lambdai(3,1)=-179/576;
Lambdai(3,2)=13/360;
Lambdai(3,3)=-11/2880;
Lambdai(3,4)=0;
Lambdai(4,1)=-5561/17280;
Lambdai(4,2)=163/3456;
Lambdai(4,3)=-23/2688;
Lambdai(4,4)=191/241920;
% Near Boundary Schemes.
% NOMENCLATURE------------------
%Thetanbs(POINT,TYPE) = returns the coefficient Theta of the NBS given
%the TYPE of scheme at point = POINT.
%Lambdanbs(POINT,TYPE, j) = returns the jth coefficient Lambda of NBS 
%scheme given the TYPE of scheme at point = POINT.
%NOTE : STARTING FROM SECOND-0RDER ACCURACY
Thetanbs(1,1)=-1;
Thetanbs(1,2)=-1;
Thetanbs(1,3)=-1;
Thetanbs(1,4)=-1;
Thetanbs(1,5)=-1;
Thetanbs(1,6)=-1;
Thetanbs(1,7)=-1;
Thetanbs(1,8)=-1;
Thetanbs(2,1)=-1/2;
Thetanbs(2,2)=-1/2;
Thetanbs(2,3)=-1/2;
Thetanbs(2,4)=-1/2;
Thetanbs(2,5)=-1/2;
Thetanbs(2,6)=-1/2;
Thetanbs(2,7)=-1/2;
Thetanbs(2,8)=-1/2;
Thetanbs(3,1)=-1/2;
Thetanbs(3,2)=-1/2;
Thetanbs(3,3)=-1/2;
Thetanbs(3,4)=-1/2;
Thetanbs(3,5)=-1/2;
Thetanbs(3,6)=-1/2;
Thetanbs(3,7)=-1/2;
Thetanbs(3,8)=-1/2;
Thetanbs(4,1)=-1/2;
Thetanbs(4,2)=-1/2;
Thetanbs(4,3)=-1/2;
Thetanbs(4,4)=-1/2;
Thetanbs(4,5)=-1/2;
Thetanbs(4,6)=-1/2;
Thetanbs(4,7)=-1/2;
Thetanbs(4,8)=-1/2;
Lambdanbs(1,1,1)=-1/2;
Lambdanbs(1,1,2)=-1/2;
Lambdanbs(1,1,3)=0;
Lambdanbs(1,1,4)=0;
Lambdanbs(1,1,5)=0;
Lambdanbs(1,1,6)=0;
Lambdanbs(1,1,7)=0;
Lambdanbs(1,1,8)=0;
Lambdanbs(1,1,9)=0;
Lambdanbs(1,2,1)=-5/12;
Lambdanbs(1,2,2)=-2/3;
Lambdanbs(1,2,3)=1/12;
Lambdanbs(1,2,4)=0;
Lambdanbs(1,2,5)=0;
Lambdanbs(1,2,6)=0;
Lambdanbs(1,2,7)=0;
Lambdanbs(1,2,8)=0;
Lambdanbs(1,2,9)=0;
Lambdanbs(1,3,1)=-3/8;
Lambdanbs(1,3,2)=-19/24;
Lambdanbs(1,3,3)=5/24;
Lambdanbs(1,3,4)=-1/24;
Lambdanbs(1,3,5)=0;
Lambdanbs(1,3,6)=0;
Lambdanbs(1,3,7)=0;
Lambdanbs(1,3,8)=0;
Lambdanbs(1,3,9)=0;
Lambdanbs(1,4,1)=-251/720;
Lambdanbs(1,4,2)=-323/360;
Lambdanbs(1,4,3)=11/30;
Lambdanbs(1,4,4)=-53/360;
Lambdanbs(1,4,5)=19/720;
Lambdanbs(1,4,6)=0;
Lambdanbs(1,4,7)=0;
Lambdanbs(1,4,8)=0;
Lambdanbs(1,4,9)=0;
Lambdanbs(1,5,1)=-95/288;
Lambdanbs(1,5,2)=-1427/1440;
Lambdanbs(1,5,3)=133/240;
Lambdanbs(1,5,4)=-241/720;
Lambdanbs(1,5,5)=173/1440;
Lambdanbs(1,5,6)=-3/160;
Lambdanbs(1,5,7)=0;
Lambdanbs(1,5,8)=0;
Lambdanbs(1,5,9)=0;
Lambdanbs(1,6,1)=-19087/60480;
Lambdanbs(1,6,2)=-2713/2520;
Lambdanbs(1,6,3)=15487/20160;
Lambdanbs(1,6,4)=-586/945;
Lambdanbs(1,6,5)=6737/20160;
Lambdanbs(1,6,6)=-263/2520;
Lambdanbs(1,6,7)=863/60480;
Lambdanbs(1,6,8)=0;
Lambdanbs(1,6,9)=0;
Lambdanbs(1,7,1)=-5257/17280;
Lambdanbs(1,7,2)=-139849/120960;
Lambdanbs(1,7,3)=4511/4480;
Lambdanbs(1,7,4)=-123133/120960;
Lambdanbs(1,7,5)=88547/120960;
Lambdanbs(1,7,6)=-1537/4480;
Lambdanbs(1,7,7)=11351/120960;
Lambdanbs(1,7,8)=-275/24192;
Lambdanbs(1,7,9)=0;
Lambdanbs(1,8,1)=-1070017/3628800;
Lambdanbs(1,8,2)=-2233547/1814400;
Lambdanbs(1,8,3)=2302297/1814400;
Lambdanbs(1,8,4)=-2797679/1814400;
Lambdanbs(1,8,5)=31457/22680;
Lambdanbs(1,8,6)=-1573169/1814400;
Lambdanbs(1,8,7)=645607/1814400;
Lambdanbs(1,8,8)=-156437/1814400;
Lambdanbs(1,8,9)=33953/3628800;
Lambdanbs(2,1,1)=1/2;
Lambdanbs(2,1,2)=-1/2;
Lambdanbs(2,1,3)=0;
Lambdanbs(2,1,4)=0;
Lambdanbs(2,1,5)=0;
Lambdanbs(2,1,6)=0;
Lambdanbs(2,1,7)=0;
Lambdanbs(2,1,8)=0;
Lambdanbs(2,1,9)=0;
Lambdanbs(2,2,1)=1/4;
Lambdanbs(2,2,2)=0;
Lambdanbs(2,2,3)=-1/4;
Lambdanbs(2,2,4)=0;
Lambdanbs(2,2,5)=0;
Lambdanbs(2,2,6)=0;
Lambdanbs(2,2,7)=0;
Lambdanbs(2,2,8)=0;
Lambdanbs(2,2,9)=0;
Lambdanbs(2,3,1)=5/24;
Lambdanbs(2,3,2)=1/8;
Lambdanbs(2,3,3)=-3/8;
Lambdanbs(2,3,4)=1/24;
Lambdanbs(2,3,5)=0;
Lambdanbs(2,3,6)=0;
Lambdanbs(2,3,7)=0;
Lambdanbs(2,3,8)=0;
Lambdanbs(2,3,9)=0;
Lambdanbs(2,4,1)=3/16;
Lambdanbs(2,4,2)=5/24;
Lambdanbs(2,4,3)=-1/2;
Lambdanbs(2,4,4)=1/8;
Lambdanbs(2,4,5)=-1/48;
Lambdanbs(2,4,6)=0;
Lambdanbs(2,4,7)=0;
Lambdanbs(2,4,8)=0;
Lambdanbs(2,4,9)=0;
Lambdanbs(2,5,1)=251/1440;
Lambdanbs(2,5,2)=79/288;
Lambdanbs(2,5,3)=-91/144;
Lambdanbs(2,5,4)=37/144;
Lambdanbs(2,5,5)=-25/288;
Lambdanbs(2,5,6)=19/1440;
Lambdanbs(2,5,7)=0;
Lambdanbs(2,5,8)=0;
Lambdanbs(2,5,9)=0;
Lambdanbs(2,6,1)=95/576;
Lambdanbs(2,6,2)=119/360;
Lambdanbs(2,6,3)=-445/576;
Lambdanbs(2,6,4)=4/9;
Lambdanbs(2,6,5)=-131/576;
Lambdanbs(2,6,6)=5/72;
Lambdanbs(2,6,7)=-3/320;
Lambdanbs(2,6,8)=0;
Lambdanbs(2,6,9)=0;
Lambdanbs(2,7,1)=19087/120960;
Lambdanbs(2,7,2)=1315/3456;
Lambdanbs(2,7,3)=-1771/1920;
Lambdanbs(2,7,4)=2399/3456;
Lambdanbs(2,7,5)=-1649/3456;
Lambdanbs(2,7,6)=421/1920;
Lambdanbs(2,7,7)=-205/3456;
Lambdanbs(2,7,8)=863/120960;
Lambdanbs(2,7,9)=0;
Lambdanbs(2,8,1)=5257/34560;
Lambdanbs(2,8,2)=1145/2688;
Lambdanbs(2,8,3)=-18689/17280;
Lambdanbs(2,8,4)=3499/3456;
Lambdanbs(2,8,5)=-7/8;
Lambdanbs(2,8,6)=9289/17280;
Lambdanbs(2,8,7)=-755/3456;
Lambdanbs(2,8,8)=101/1920;
Lambdanbs(2,8,9)=-275/48384;
Lambdanbs(3,1,1)=1/2;
Lambdanbs(3,1,2)=-1/2;
Lambdanbs(3,1,3)=0;
Lambdanbs(3,1,4)=0;
Lambdanbs(3,1,5)=0;
Lambdanbs(3,1,6)=0;
Lambdanbs(3,1,7)=0;
Lambdanbs(3,1,8)=0;
Lambdanbs(3,1,9)=0;
Lambdanbs(3,2,1)=-1/4;
Lambdanbs(3,2,2)=1;
Lambdanbs(3,2,3)=-3/4;
Lambdanbs(3,2,4)=0;
Lambdanbs(3,2,5)=0;
Lambdanbs(3,2,6)=0;
Lambdanbs(3,2,7)=0;
Lambdanbs(3,2,8)=0;
Lambdanbs(3,2,9)=0;
Lambdanbs(3,3,1)=-1/24;
Lambdanbs(3,3,2)=3/8;
Lambdanbs(3,3,3)=-1/8;
Lambdanbs(3,3,4)=-5/24;
Lambdanbs(3,3,5)=0;
Lambdanbs(3,3,6)=0;
Lambdanbs(3,3,7)=0;
Lambdanbs(3,3,8)=0;
Lambdanbs(3,3,9)=0;
Lambdanbs(3,4,1)=-1/48;
Lambdanbs(3,4,2)=7/24;
Lambdanbs(3,4,3)=0;
Lambdanbs(3,4,4)=-7/24;
Lambdanbs(3,4,5)=1/48;
Lambdanbs(3,4,6)=0;
Lambdanbs(3,4,7)=0;
Lambdanbs(3,4,8)=0;
Lambdanbs(3,4,9)=0;
Lambdanbs(3,5,1)=-19/1440;
Lambdanbs(3,5,2)=73/288;
Lambdanbs(3,5,3)=11/144;
Lambdanbs(3,5,4)=-53/144;
Lambdanbs(3,5,5)=17/288;
Lambdanbs(3,5,6)=-11/1440;
Lambdanbs(3,5,7)=0;
Lambdanbs(3,5,8)=0;
Lambdanbs(3,5,9)=0;
Lambdanbs(3,6,1)=-3/320;
Lambdanbs(3,6,2)=83/360;
Lambdanbs(3,6,3)=77/576;
Lambdanbs(3,6,4)=-4/9;
Lambdanbs(3,6,5)=67/576;
Lambdanbs(3,6,6)=-11/360;
Lambdanbs(3,6,7)=11/2880;
Lambdanbs(3,6,8)=0;
Lambdanbs(3,6,9)=0;
Lambdanbs(3,7,1)=-863/120960;
Lambdanbs(3,7,2)=3713/17280;
Lambdanbs(3,7,3)=347/1920;
Lambdanbs(3,7,4)=-1807/3456;
Lambdanbs(3,7,5)=673/3456;
Lambdanbs(3,7,6)=-149/1920;
Lambdanbs(3,7,7)=337/17280;
Lambdanbs(3,7,8)=-271/120960;
Lambdanbs(3,7,9)=0;
Lambdanbs(3,8,1)=-275/48384;
Lambdanbs(3,8,2)=24587/120960;
Lambdanbs(3,8,3)=85/384;
Lambdanbs(3,8,4)=-10439/17280;
Lambdanbs(3,8,5)=8/27;
Lambdanbs(3,8,6)=-61/384;
Lambdanbs(3,8,7)=1039/17280;
Lambdanbs(3,8,8)=-335/24192;
Lambdanbs(3,8,9)=13/8960;
Lambdanbs(4,1,1)=1/2;
Lambdanbs(4,1,2)=-1/2;
Lambdanbs(4,1,3)=0;
Lambdanbs(4,1,4)=0;
Lambdanbs(4,1,5)=0;
Lambdanbs(4,1,6)=0;
Lambdanbs(4,1,7)=0;
Lambdanbs(4,1,8)=0;
Lambdanbs(4,1,9)=0;
Lambdanbs(4,2,1)=-3/4;
Lambdanbs(4,2,2)=2;
Lambdanbs(4,2,3)=-5/4;
Lambdanbs(4,2,4)=0;
Lambdanbs(4,2,5)=0;
Lambdanbs(4,2,6)=0;
Lambdanbs(4,2,7)=0;
Lambdanbs(4,2,8)=0;
Lambdanbs(4,2,9)=0;
Lambdanbs(4,3,1)=5/24;
Lambdanbs(4,3,2)=-7/8;
Lambdanbs(4,3,3)=13/8;
Lambdanbs(4,3,4)=-23/24;
Lambdanbs(4,3,5)=0;
Lambdanbs(4,3,6)=0;
Lambdanbs(4,3,7)=0;
Lambdanbs(4,3,8)=0;
Lambdanbs(4,3,9)=0;
Lambdanbs(4,4,1)=1/48;
Lambdanbs(4,4,2)=-1/8;
Lambdanbs(4,4,3)=1/2;
Lambdanbs(4,4,4)=-5/24;
Lambdanbs(4,4,5)=-3/16;
Lambdanbs(4,4,6)=0;
Lambdanbs(4,4,7)=0;
Lambdanbs(4,4,8)=0;
Lambdanbs(4,4,9)=0;
Lambdanbs(4,5,1)=11/1440;
Lambdanbs(4,5,2)=-17/288;
Lambdanbs(4,5,3)=53/144;
Lambdanbs(4,5,4)=-11/144;
Lambdanbs(4,5,5)=-73/288;
Lambdanbs(4,5,6)=19/1440;
Lambdanbs(4,5,7)=0;
Lambdanbs(4,5,8)=0;
Lambdanbs(4,5,9)=0;
Lambdanbs(4,6,1)=11/2880;
Lambdanbs(4,6,2)=-13/360;
Lambdanbs(4,6,3)=179/576;
Lambdanbs(4,6,4)=0;
Lambdanbs(4,6,5)=-179/576;
Lambdanbs(4,6,6)=13/360;
Lambdanbs(4,6,7)=-11/2880;
Lambdanbs(4,6,8)=0;
Lambdanbs(4,6,9)=0;
Lambdanbs(4,7,1)=271/120960;
Lambdanbs(4,7,2)=-433/17280;
Lambdanbs(4,7,3)=533/1920;
Lambdanbs(4,7,4)=191/3456;
Lambdanbs(4,7,5)=-1265/3456;
Lambdanbs(4,7,6)=133/1920;
Lambdanbs(4,7,7)=-257/17280;
Lambdanbs(4,7,8)=191/120960;
Lambdanbs(4,7,9)=0;
Lambdanbs(4,8,1)=13/8960;
Lambdanbs(4,8,2)=-2267/120960;
Lambdanbs(4,8,3)=883/3456;
Lambdanbs(4,8,4)=191/1920;
Lambdanbs(4,8,5)=-91/216;
Lambdanbs(4,8,6)=1961/17280;
Lambdanbs(4,8,7)=-71/1920;
Lambdanbs(4,8,8)=191/24192;
Lambdanbs(4,8,9)=-191/241920;


%Checking if the depth of the column operand is correctly specified. 
if(z<1)
    MESS_TEXT = ['The depth of discrete space (the second argument) '...
        'is invalid. Use z>=1.'];
    error(MESS_TEXT);
end

%Evaluating the length of the colomn operand.
if(colength <3)
    MESS_TEXT = ['The length of discrete space (the first arg) is not '...
        'valid. Use colength>=1'];
    error(MESS_TEXT);
end

%------------------------------------------------------------------------
%Generating the LHS tensor [[A]].
%------------------------------------------------------------------------

% Local variables.
A = cell(colength,colength);
A(:) = {zeros(z)}; 
B = cell(colength,colength);
B(:) = {zeros(z)};

%Making tensor [[A]]...
%Defining Near Boundary Schemes for [[A]] at point i=1.
%Theta1 = 2;
Theta1 = -1;
A(1,1) = {eye(z)};
A(1,2) = {Theta1*eye(z)};

%Defining Interior Schemes for [[A]].
Theta1 = -1/2;
for i = 2:(colength-1)
    A(i,i-1) = {Theta1*eye(z)};
    A(i,i) = {eye(z)};
    A(i,i+1) = {Theta1*eye(z)};
end
    
%Defining Near Boundary Schemes for [[A]] at point i="colength".
Theta1 = -1;
A(colength,colength) = {eye(z)};
A(colength,colength-1) = {Theta1*eye(z)};

%------------------------------------------------------------------------
%Generating the RHS tensor [[B]].
%------------------------------------------------------------------------

%Defining Interior Schemes for [[B]].

for i = (INTSCHTYPE+1):(colength-INTSCHTYPE)
    for j = 1:INTSCHTYPE
    B(i,i-j) = {-Lambdai(INTSCHTYPE,j)*eye(z)};
    B(i,i+j) = {Lambdai(INTSCHTYPE,j)*eye(z)};
    end
end

%Defining Near Boundary Schemes for [[B]] at point i=1.

for i = 1:length(NBSCHTYPE)
    for j = 1:(max(max(NBSCHTYPE))+1)
        B(i,j) = {Lambdanbs(i,NBSCHTYPE(i),j)*eye(z)};
    end
end
        
% Defining Near Boundary Schemes for [[B]] at points near i="colength".
for i = 1:length(NBSCHTYPE)
    for j = 1:(max(max(NBSCHTYPE))+1)
        B(colength-i+1,colength-j+1) = {-Lambdanbs(i,NBSCHTYPE(i),j)*eye(z)};
    end
end

%converting [[A]] and [[B]] from cell array to matrices 
A = cell2mat(A);
B = cell2mat(B);

%slicing the matrices
A(:,1:z) = []; %removing the first column of [[A]].
A((end-z+1):end,:) = []; %removing the last row of [[A]].

f1 = B(:,1:z);
f1((end-z+1):end,:) = []; %removing the last entity
B(:,1:z) = []; %removing the first column of [[B]].
B((end-z+1):end,:) = []; %removing the last row of [[B]].



























