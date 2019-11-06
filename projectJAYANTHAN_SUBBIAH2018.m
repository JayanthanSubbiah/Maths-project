fprintf('\n***************question1************\n');
syms s X t;

Y = dsolve('D2x + 2*Dx + 26*x = 10*cos(t)*heaviside(t-pi)','x(0) = 0.5, Dx(0) = 1'); 
%"dsolve" is used to solve the differential equations
%"heaviside" is a step function which returns value
pretty(Y);

%laplace Transform for heaviside function
F = 10*(cos(t))*heaviside(t-pi); 
FO = laplace(F, t, s);
%"t"is independent variable
%"s"is transformation variable
%initialising the value
X1 = s*X - 0.5;
X2 = s*X1 - 1; 
%Inverse laplace transformation
Laplace = solve(X2 + 2*X1 + 26*X - FO, X);
F = ilaplace(Laplace, s, t);

%graphical representation
%"flot" funtion to plot graph in the specified time intreval
%"axis" to scale the graph
figure(1)
fplot(F, [0, 20]);
axis([0 20 -0.8 0.8]);


fprintf('\n***************question2************\n');
clear;
syms Y0 Y1 Y2 X1;
%creating a table
T = table([1;2;3;4;5;6],[2;3;1;5;4;1],'VariableNames',{'X','Y'} );
%assignning the values to the variables
h=1; 
x=4;
%forward
%f'(x)=[f(x+h)-f(x)]/h
forward=(T.Y(x + h) - T.Y(x))/h;
fprintf('Forward=%f\n',forward); 

%backward
%f'(x)=[f(x)-f(x-h)]/h
backward=(T.Y(x) - T.Y(x - h))/h;
fprintf('Backward=%f\n',backward);  

%centered   
%f'(x)=[f(x+h)-f(x-h)]/2h
centered=(T.Y(x + h) - T.Y(x - h))/2*h;
fprintf('Centered=%f\n',centered);

%second oreder polynomial interpolation   
eqn1 = Y0 + Y1*T.X(x - 1) + Y2*((T.X(x - 1))^2) == T.Y(x - 1);          
eqn2 = Y0 + Y1*T.X(x) + Y2*((T.X(x))^2) == T.Y(x);                                          
eqn3 = Y0 + Y1*T.X(x + 1) + Y2*((T.X(x + 1))^2) == T.Y(x + 1);  

%converting equation to the matrix format
[A,B] = equationsToMatrix([eqn1, eqn2, eqn3], [Y0, Y1, Y2]);

%"linsolve"to solve the linear equation
fprintf('\nsolving linear equation \n');
C = linsolve(A,B);   
fprintf('\tC0=%f\tC1=%f\tC2=%f\n',C);
P(X1) = C(1,1) + X1*C(2,1) + (X1^2)*C(3,1);
PO(X1) = diff(P, X1);
%estimating f'(4)
fprintf("f'(4)=%f\n",PO(4));

%trapezoidal rule 
S=0;
A=1;
B=6;
N=6;
for i = 1:N-1
    S = S+(T.Y(i) + T.Y(i + 1));
end
fprintf('\nusing for loop= %f\n', (S * h/2));
fprintf('Trapezoidal %f\n', trapz(T.Y));
fprintf('The output of both for loop and build in function are same\n');


fprintf('\n***************question3************\n');
clear
syms NumIters p
format long
%%System Starting Point
a=[1;0];
 
%given function
y = @(x)[x(1)^3+x(1)*x(2)+x(2)^3-2; x(1)*3-x(1)^5*x(2)+x(2)^3-3];
%differentiation
dy=@(x)[3*x(1)^2+3*x(2)^2,2;3,3*x(1)^2+5*x(1)^4+3*x(2)^2];
%itterations
n=6;
for i=1:n-1
    dx=-dy(a)\y(a); %solve the icrement
    a=a+dx;         %new starting
    y(a)            %to check zero
end
%figure(2)
%plot(y(a),dy(a));

%using build in function
x1 = fsolve(y,a);
fprintf('Newtons method using buildin function=%f\n',x1);



fprintf('\n***************question4************\n');

clear
%initialization
alpha1 = 0.5 , alpha2 = 0.9 , alpha3 = 0.95 , alpha4 = 0.99;  
OA1 = sqrt(1 - alpha1*alpha1);
OA2 = sqrt(1 - alpha2*alpha2);
OA3 = sqrt(1 - alpha3*alpha3);
OA4 = sqrt(1 - alpha4*alpha4);
%normal distribution of random numbers
Z1(1) = randn(); 
Z2(1) = randn();
Z3(1) = randn();
Z4(1) = randn();
Z5(1) = Z1(1) +  Z2(1) + Z3(1) + Z4(1);
%recursively
for i=(2:500)
    Z1(i) = alpha1 * Z1(i-1) + OA1 * randn();
    Z2(i) = alpha2 * Z2(i-1) + OA2 * randn();
    Z3(i) = alpha3 * Z3(i-1) + OA3 * randn();
    Z4(i) = alpha4 * Z4(i-1) + OA4 * randn();
    Z5(i) = Z1(i) +  Z2(i) + Z3(i) +Z4(i);
end
figure(3)
subplot(2,3,1)  
plot((1:500), Z1);
title('0.5')
subplot(2,3,2)  
plot((1:500), Z2);
title('0.9')        
subplot(2,3,3) 
plot((1:500), Z3);
title('0.95')
subplot(2,3,4)
plot((1:500), Z4);
title('0.99')
subplot(2,3,5) 
plot((1:500), Z5);