function dx = nonlinear_pendulum(t,x)
dx = [
x(2);
-sin(x(1));
];