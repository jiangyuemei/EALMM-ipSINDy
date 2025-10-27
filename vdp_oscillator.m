function dx = vdp_oscillator(t,x)
dx = [x(2);
4*x(2)-x(1)-4*x(2)*x(1)^2;
];