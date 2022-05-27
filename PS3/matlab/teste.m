syms x y
fplot(chebyshevT(0:4,x))
axis([-1 1 -1 1])
grid on

title('Chebyshev polynomials of the first kind')