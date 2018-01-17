function q1
	close all; clear;
	h = 0.05;
	k = 0.3 * h^2;
	m = 1/h;
	n = 50;

	X = (0:n)*k;
	Y = (0:m)*h;

	A = zeros(n+1, m+1);
	for i = 1:n+1
		for j = 1:m+1
			A(i,j) = exp(-pi^2*k*(i-1)) * sin(pi*(j-1)*h);
		end
	end
	surf(A)

	U = FTCS(h, k, m, n, @fun, @f, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)');
	hold on; plot(Y, A(end, :)); hold off;
	% saveas(gcf, 'plots/q1_1.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)');
	hold on; surf(X, Y, A'); hold off;
	% saveas(gcf, 'plots/q1_2.png');
	
	U = BTCS(h, k, m, n, @fun, @f, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)');
	hold on; plot(Y, A(end, :)); hold off;
	% saveas(gcf, 'plots/q1_3.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)');
	hold on; surf(X, Y, A'); hold off;
	% saveas(gcf, 'plots/q1_4.png');

	U = Crank(h, k, m, n, @fun, @f, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)');
	hold on; plot(Y, A(end, :)); hold off;
	% saveas(gcf, 'plots/q1_3.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)');
	hold on; surf(X, Y, A'); hold off;
	% saveas(gcf, 'plots/q1_4.png');
end

function [y] = fun(x, t)
	y = 0;
end

function [y] = f(x)
	y = sin(pi*x);
end

function [y] = g1(t)
	y = 0;
end

function [y] = g2(t)
	y = 0;
end

function [U] = FTCS(h, k, m, n, fun, f, g1, g2)
	fprintf('\nRunning FTCS\n');
	lamda = k / h^2;
	U = zeros(n+1, m+1);

	U(1, 1:end) = f((0:m)*h);
	U(1:end, 1) = g1((0:n)*k);
	U(1:end, end) = g2((0:n)*k);

	for i = 2:n+1
		for j = 2:m
			t = (i-1)*k;
			x = (j-1)*h;
			U(i, j) = lamda*U(i-1,j-1) + (1-2*lamda)*U(i-1,j) + lamda*U(i-1,j+1) + k*fun(x,t);
		end
	end

	U;
end

function [U] = BTCS(h, k, m, n, fun, f, g1, g2)
	fprintf('\nRunning BTCS\n');
	lamda = k / h^2;
	U = zeros(n+1, m+1);

	U(1, 1:end) = f((0:m)*h);
	U(1:end, 1) = g1((0:n)*k);
	U(1:end, end) = g2((0:n)*k);

	for i = 2:n+1
		A = zeros(m+1, m+1);
		b = zeros(m+1, 1);

		A(1:m+2:end) = 1 + 2*lamda;
		A(2:m+2:end) = -lamda;
		A(m+2:m+2:end) = -lamda;

		A(1,1) = 1;
		A(1,2) = 0;
		A(m+1,m+1) = 1;
		A(m+1,m) = 0;

		b(2:m) = U(i-1,2:m) + k*fun((1:m-1)*h,(i-1)*k);
		b(1) = U(i,1);
		b(end) = U(i,end);

		U(i,:) = (A\b)';
	end

	U;
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
		A(m+1,m+1) = 1;
		A(m+1,m) = 0;

		b(2:m) = U(i-1,1:m-1)*lamda/2 + (1-lamda)*U(i-1,2:m) + U(i-1,3:m+1)*lamda/2 + k*fun((1:m-1)*h,(i-1)*k);
		b(1) = U(i,1);
		b(end) = U(i,end);

		U(i,:) = (A\b)';
	end

	U;
end