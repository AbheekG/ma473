function q1
	close all; clear; im_num = 1;
	% Terminal Condition flag
	isTerminal = true;

	T = 1;
	K = 10;
	r = 0.06;
	sig = 0.3;
	delta = 0;

	% Boundary
	S_min = 0;
	S_max = 30;

	h = 0.5;
	k = h^2/2.1;
	m = (S_max - S_min)/h;
	n = ceil(T/k);

	S = S_min:h:S_max;
	Time = 0:k:T;

	% Tau = T - t
	U = FTCS(T, K, r, sig, delta, S, Time, h, k, isTerminal);
	% length(X), length(U(end, :))
	U(1, :)
	figure; plot(S, U(1, :)); hold on; plot(S, U(1, :)); hold off;
	legend('Cost of option at t = 0', 'Cost of option at t = T'); xlabel('S'); ylabel('u(S, t)'); title('FTCS');
	saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
	figure; surf(S, Time, U); xlabel('S'); ylabel('t'); zlabel('u(x,t)'); title('FTCS');
	saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
	
	% Methods = ['Direct'; 'GaussS'; 'Jacobi'; 'SOR   '];
	% for meth = 1:4
	% 	U = BTCS(@fun, @f, @g1, @g2, T, K, r, sig, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, Methods(meth, :));
	% 	figure; plot(S, U(end, :)); hold on; plot(S, U(1, :)); hold off;
	% 	legend('Cost of option at t = 0', 'Cost of option at t = T'); xlabel('S'); ylabel('u(S, t)'); title(sprintf('BTCS using %s method', Methods(meth, :)));
	% 	saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
	% 	figure; surf(S, Time, U); xlabel('S'); ylabel('t'); zlabel('u(S,t)'); title(sprintf('BTCS using %s method', Methods(meth, :)));
	% 	saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
	% end

	% for meth = 1:4
	% 	U = Crank(@fun, @f, @g1, @g2, T, K, r, sig, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, Methods(meth, :));
	% 	figure; plot(S, U(end, :)); hold on; plot(S, U(1, :)); hold off;
	% 	legend('Cost of option at t = 0', 'Cost of option at t = T'); xlabel('S'); ylabel('u(S, t)'); title(sprintf('Crank-Nicolson using %s method', Methods(meth, :)));
	% 	saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
	% 	figure; surf(S, Time, U); xlabel('S'); ylabel('t'); zlabel('u(S,t)'); title(sprintf('Crank-Nicolson using %s method', Methods(meth, :)));
	% 	saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
	% end
end

function [y] = f(S, K)
	temp1 = zeros(size(S));
	temp2 = S - K;
	y = max([temp1; temp2]);
end

function [y] = g1(r, s, Tau, K)
	y = 0;
end

function [y] = g2(r, s, Tau, K)
	y = s - K*exp(-r.*Tau);
end

function [y] = fa(sig, s)
	y = sig.^2 .* s.^2 / 2;
end

function [y] = fb(r, delta, s);
	y = (r - delta) .* s;
end

function [y] = fc(r)
	y = -r;
end

function [U] = FTCS(T, K, r, sig, delta, S, Tau, h, k, isTerminal)
	fprintf('\nRunning FTCS\n');
	m = length(S);
	n = length(Time);
	U = zeros(n, m);

	if isTerminal
		k = -k;
	end

	U(1:end, 1) = g1(r, S(1), Time, K);
	U(1:end, end) = g2(r, S(end), Time, K);
	U(1, 1:end) = f(S, K);

	for i = 2:n+1
		for j = 2:m-1
			aa = fa(sig, S(j));
			bb = fb(r, delta, S(j));
			cc = fc(r);
			U(i, j) = (-aa*k/h^2 + 0.5*bb*k/h)*U(i-1,j-1) + (1 + 2*aa*k/h^2 - cc*k)*U(i-1,j) + (-aa*k/h^2 - 0.5*bb*k/h)*U(i-1,j+1);
		end
	end

	if isTerminal
		U = flipud(U);
	end
end

function [U] = BTCS(fun, f, g1, g2, T, K, r, sig, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, method)
	fprintf('\nRunning BTCS\n');
	fprintf('Using %s method\n', method);
	m = length(S);
	n = length(Time);
	U = zeros(n, m);

	if isTerminal
		k = -k;
		Time = flip(Time);
	end

	U(1:end, 1) = g1(r, s, t, T);
	U(1:end, end) = g2(r, s, t, T);
	U(1, 1:end) = f(S, K);

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

		b(2:m) = U(i-1,2:m);
		b(1) = U(i,1);
		b(end) = U(i,end);

		if method == 'Direct'
			U(i,:) = (A\b)';
		elseif method == 'GaussS'
			U(i,:) = gauss_seidel(A,b,1000,1e-5);
		elseif method == 'Jacobi'
			U(i,:) = jacobi(A,b,1000,1e-5);
		else
			U(i,:) = sor(A,b,1000,1e-5);
		end			
			
	end

	U;
end

function [U] = Crank(fun, f, g1, g2, T, K, r, sig, delta, q, qd, x_min, x_max, h, k, m, n, X, Tau, method)
	fprintf('\nRunning Crank Nicolson\n');
	fprintf('Using %s method\n', method);
	m = length(S);
	n = length(Time);
	U = zeros(n, m);

	if isTerminal
		k = -k;
		Time = flip(Time);
	end

	U(1:end, 1) = g1(r, s, t, T);
	U(1:end, end) = g2(r, s, t, T);
	U(1, 1:end) = f(S, K);

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

		b(2:m) = U(i-1,1:m-1)*lamda/2 + (1-lamda)*U(i-1,2:m) + U(i-1,3:m+1)*lamda/2;
		b(1) = U(i,1);
		b(end) = U(i,end);

		if method == 'Direct'
			U(i,:) = (A\b)';
		elseif method == 'GaussS'
			U(i,:) = gauss_seidel(A,b,1000,1e-5);
		elseif method == 'Jacobi'
			U(i,:) = jacobi(A,b,1000,1e-5);
		else
			U(i,:) = sor(A,b,1000,1e-5);
		end	
	end

	U;
end