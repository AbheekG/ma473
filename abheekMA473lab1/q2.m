function q2
	close all; clear;
	h = 0.025;
	k = h^2/2.1;
	m = 1/h;
	n = 3;

	X = (0:n)*k;
	Y = (0:m)*h;

	fprintf('f(x) = sin(pi*x)');
	U = FTCS(h, k, m, n, @fun, @f1, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); title('FTCS');
	saveas(gcf, 'plots/q2_1.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)'); title('FTCS');
	saveas(gcf, 'plots/q2_2.png');
	
	U = BTCS(h, k, m, n, @fun, @f1, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); title('BTCS');
	saveas(gcf, 'plots/q2_3.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)'); title('BTCS');
	saveas(gcf, 'plots/q2_4.png');

	U = CTCS2(h, k, m, n, @fun, @f1, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); title('CTCS');
	% U = CTCS1(h, k, m, n, @fun, @f1, @g1, @g2);
	% hold on; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); hold off;
	saveas(gcf, 'plots/q2_5.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)'); title('CTCS');
	saveas(gcf, 'plots/q2_6.png');

	fprintf('f(x) = x*(1-x)');
	U = FTCS(h, k, m, n, @fun, @f2, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); title('FTCS');
	saveas(gcf, 'plots/q2_7.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)'); title('FTCS'); title('FTCS');
	saveas(gcf, 'plots/q2_8.png');
	
	U = BTCS(h, k, m, n, @fun, @f2, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); title('BTCS');
	saveas(gcf, 'plots/q2_9.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)'); title('BTCS'); title('BTCS');
	saveas(gcf, 'plots/q2_10.png');

	U = CTCS2(h, k, m, n, @fun, @f2, @g1, @g2);
	figure; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); title('CTCS');
	% U = CTCS1(h, k, m, n, @fun, @f2, @g1, @g2);
	% hold on; plot(Y, U(end, :)); xlabel('x'); ylabel('u(x, T)'); hold off;
	saveas(gcf, 'plots/q2_11.png');
	figure; surf(X, Y, U'); xlabel('t'); ylabel('x'); zlabel('u(t,x)'); title('CTCS');
	saveas(gcf, 'plots/q2_12.png');
end

function [y] = fun(x, t)
	y = 0;
end

function [y] = f1(x)
	y = sin(pi*x);
end

function [y] = f2(x)
	y = x.*(1-x);
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

function [U] = CTCS1(h, k, m, n, fun, f, g1, g2)
	fprintf('\nRunning CTCS implementation 1\n');
	lamda = k / h^2;
	U = zeros(n+1, m+1);

	U(1, 1:end) = f((0:m)*h);
	U(1:end, 1) = g1((0:n)*k);
	U(1:end, end) = g2((0:n)*k);

	for i = 2:n+1
		if i == 2
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
		else
			for j = 2:m
				t = (i-2)*k;
				x = (j-1)*h;
				U(i, j) = U(i-2,j) + 2*lamda*U(i-1,j-1) - 4*lamda*U(i-1,j) + 2*lamda*U(i-1,j+1) + 2*k*fun(x,t);
			end
		end
	end

	U;
end

function [U] = CTCS2(h, k, m, n, fun, f, g1, g2)
	fprintf('\nRunning CTCS implementation 2\n');
	lamda = k / h^2;
	U = zeros(n+1, m+1);

	A = zeros((m+1)*(n+1), (m+1)*(n+1));
	b = zeros((m+1)*(n+1), 1);
	U(1, 1:end) = f((0:m)*h);
	U(1:end, 1) = g1((0:n)*k);
	U(1:end, end) = g2((0:n)*k);

	for i = 1:n+1
		for j = 1:m+1
			% (i-1)(m+1) + j
			t = (i-1)*k;
			x = (j-1)*h;
			if (i >= 2) && (j >= 2) && (i <= n) && (j <= m)
				A((i-1)*(m+1) + j, i*(m+1) + j) = 1;
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j-1) = -2*lamda;
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j) = 4*lamda;
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j+1) = -2*lamda;
				A((i-1)*(m+1) + j, (i-2)*(m+1) + j) = -1;

				b((i-1)*(m+1) + j) = 2*k*fun(x,t);
			elseif i == 1
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j) = 1;
				b((i-1)*(m+1) + j) = f((j-1)*h);
			elseif j == 1
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j) = 1;
				b((i-1)*(m+1) + j) = g1((i-1)*k);
			elseif j == m + 1
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j) = 1;
				b((i-1)*(m+1) + j) = g2((i-1)*k);
			else
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j) = 1 + 2*lamda;
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j-1) = -lamda;
				A((i-1)*(m+1) + j, (i-1)*(m+1) + j+1) = -lamda;
				A((i-1)*(m+1) + j, (i-2)*(m+1) + j) = -1;
				b((i-1)*(m+1) + j) = k*fun(x,t);
			end
		end			
		
		% A, b
		X = A\b;
		U = reshape(X, m+1, n+1)';
	end

	U;
end