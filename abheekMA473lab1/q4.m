function q4
	close all; clear;
	h = 0.25;
	k = 0.05;
	m = 1/h;
	n = 50;

	X = (0:n)*k;
	Y = (0:m)*h;
	U = Crank(h, k, m, n, @fun, @f, @g1, @g2);
	% size(X), size(Y), size(U)
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); title('Crank-Nicolson');
	saveas(gcf, 'plots/q4_1.png');
	figure; surf(X, Y, U'); surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)'); title('Crank-Nicolson');
	saveas(gcf, 'plots/q4_2.png');
end

function [y] = fun(x, t)
	y = 0;
end

function [y] = f(x)
	y = x.*(1-x);
end

function [y] = g1(t)
	y = 0;
end

function [y] = g2(t)
	y = t.^2;
end

function [U] = Crank(h, k, m, n, fun, f, g1, g2)
	fprintf('\nRunning Crank Nicolson\n');
	lamda = k / h^2;
	U = zeros(n+1, m+1);

	U(1, 1:end) = f((0:m)*h);
	U(1:end, 1) = g1((0:n)*k);
	U(1:end, end) = g2((0:n)*k);

	for i = 2:n+1
		A = zeros(m+1, m+1);
		b = zeros(m+1, 1);

		A(1:m+2:end) = 1 + lamda;
		A(2:m+2:end) = -lamda/2;
		A(m+2:m+2:end) = -lamda/2;

		A(1,1) = 1;
		A(1,2) = 0;
		A(m+1,m+1) = 1/h;
		A(m+1,m) = -1/h;

		b(2:m) = U(i-1,1:m-1)*lamda/2 + (1-lamda)*U(i-1,2:m) + U(i-1,3:m+1)*lamda/2 + k*fun((1:m-1)*h,(i-1)*k);
		b(1) = U(i,1);
		b(end) = U(i,end);

		U(i,:) = (A\b)';
	end

	U;
end