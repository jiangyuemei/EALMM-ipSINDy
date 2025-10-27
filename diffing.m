function dx = diffing(t,x)
dx = [
x(2);
-0.2*x(2)-0.2*x(1)-x(1)^3;
];