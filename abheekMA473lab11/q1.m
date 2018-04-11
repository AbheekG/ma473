function q1
	close all; clear; im_num = 1;
	% Option Parameters
	T = 1;
	K = 10;
	r = 0.06;
	sig = 0.3;
	delta = 0;

	q = 2*r/sig^2;
	qd = 2*(r-delta)/sig^2;

	% Computational Parameters
	x_max = 1;
	x_min = -5;

	h = 0.05;
	k = h^2/2;
	m = (x_max - x_min)/h;
	n = ceil((T*sig^2/2)/k);

	X = x_min:h:x_max;
	Tau = 0:k:T*sig^2/2;

	S = K*exp(X);
	Time = T - 2*Tau/sig^2;

	Methods = ['Trapezoidal rule with piecewise linear functions'; 'Simpson’s  rule  with piecewise linear functions'];
	for meth = 1:2
		U = Crank(@fun, @f, @g1, @g2, T, K, r, sig, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, Methods(meth, :));
		figure; plot(S, U(end, :)); hold on; plot(S, U(1, :)); hold off;
		legend('Cost of option at t = 0', 'Cost of option at t = T'); xlabel('S'); ylabel('u(S, t)'); title(sprintf('Crank-Nicolson using %s method', Methods(meth, :)));
		saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
		figure; surf(S, Time, U); xlabel('S'); ylabel('t'); zlabel('u(S,t)'); title(sprintf('Crank-Nicolson using %s method', Methods(meth, :)));
		saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
	end
end

function [y] = fun(x, t)
	y = 0;
end

function [y] = f(x, qd)
	temp1 = zeros(size(x));
	temp2 = exp(x*(qd + 1)/2 ) - exp(x*(qd - 1)/2);
	y = max([temp1; temp2]);
end

function [y] = g1(x, t, qd)
	y = 0;
end

function [y] = g2(x, t, qd)
	y = exp(x.*(qd + 1)/2 + t.*(qd + 1)^2/4);
end

function [y] = transform(U, X, Tau, q, qd, K)
	y = zeros(size(U));
	if length(Tau) == 1
		for j = 1:length(X)
			y(j) = U(j) * K * exp(-0.5* (qd-1)*X(j) - (0.25*(qd-1)^2 + q)*Tau);
		end
	else
		for i = 1:length(Tau)
			for j = 1:length(X)
				y(i, j) = U(i, j) * K * exp(-0.5* (qd-1)*X(j) - (0.25*(qd-1)^2 + q)*Tau(i));
			end
		end
	end
end


function [U] = Crank(fun, f, g1, g2, T, K, r, sig, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, method)
	fprintf('\nRunning Crank Nicolson\n');
	fprintf('Using %s method\n\n', method);
	lamda = k / h^2;
	U = zeros(n+1, m+1);

	U(1:end, 1) = g1(x_min, Tau, qd);
	U(1:end, end) = g2(x_max, Tau, qd);
	U(1, 1:end) = f(X, qd);

	A = zeros(m-1, m-1);
	B = zeros(m-1, m-1);

	


	for i = 2:n+1
		b = zeros(m+1, 1);

		A(1:m+2:end) = 1 + lamda;
		A(2:m+2:end) = -lamda/2;
		A(m+2:m+2:end) = -lamda/2;

		A(1,1) = 1;
		A(1,2) = 0;
		A(m+1,m+1) = 1;
		A(m+1,m) = 0;

		b(2:m) = U(i-1,1:m-1)*lamda/2 + (1-lamda)*U(i-1,2:m) + U(i-1,3:m+1)*lamda/2;
		b(1) = U(i,1);
		b(end) = U(i,end);

		U(i,:) = (A\b)';
	end

	U = transform(U, X, Tau, q, qd, K);
end

% function [U] = Crank(fun, f, g1, g2, T, K, r, sig, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, method)
% 	fprintf('\nRunning Crank Nicolson\n');
% 	fprintf('Using %s method\n\n', method);
% 	lamda = k / h^2;
% 	U = zeros(n+1, m+1);

% 	U(1:end, 1) = g1(x_min, Tau, qd);
% 	U(1:end, end) = g2(x_max, Tau, qd);
% 	U(1, 1:end) = f(X, qd);

% 	for i = 2:n+1
% 		A = zeros(m+1, m+1);
% 		b = zeros(m+1, 1);

% 		A(1:m+2:end) = 1 + lamda;
% 		A(2:m+2:end) = -lamda/2;
% 		A(m+2:m+2:end) = -lamda/2;

% 		A(1,1) = 1;
% 		A(1,2) = 0;
% 		A(m+1,m+1) = 1;
% 		A(m+1,m) = 0;

% 		b(2:m) = U(i-1,1:m-1)*lamda/2 + (1-lamda)*U(i-1,2:m) + U(i-1,3:m+1)*lamda/2;
% 		b(1) = U(i,1);
% 		b(end) = U(i,end);

% 		U(i,:) = (A\b)';
% 	end

% 	U = transform(U, X, Tau, q, qd, K);
% end
