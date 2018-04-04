function q1
	close all; clear; im_num = 1;
	% Terminal Condition flag
	isTerminal = true;

	T = 0.2;
	K = 100;
	r = 0.05;
	sig = 0.05;

	% Boundary
	R_min = 0;
	R_max = 1;

	h = 0.01;
	k = 0.01;
	m = (R_max - R_min)/h;
	n = ceil(T/k);

	R = R_min:h:R_max;
	Time = 0:k:T;

	H = Crank(T, K, r, sig, R, Time, h, k, isTerminal);
	figure; plot(R, H(1, :)); hold on; plot(R, H(end, :)); hold off;
	legend('H(R,0) at t = 0', 'H(R,T) at t = T'); xlabel('R'); ylabel('H(R, t)'); title('H(R,t) using Crank-Nicolson');
	saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
	figure; surf(R, Time, H); xlabel('R'); ylabel('t'); zlabel('H(R,t)'); title('H(R,t) using Crank-Nicolson');
	saveas(gcf, sprintf('plots/q1_%d.png', im_num)); im_num = im_num + 1;
end

function [y] = f(R, T)
	temp1 = zeros(size(R));
	temp2 = 1 - R./T;
	y = max([temp1; temp2]);
end

function [y] = fa(sig, R)
	y = sig.^2 .* R.^2 / 2;
end

function [y] = fb(r, R);
	y = 1 - r*R;
end

function [y] = fc(r)
	y = 0;
end

function [U] = Crank(T, K, r, sig, R, Tau, h, k, isTerminal)
	fprintf('\nRunning Crank Nicolson\n');
	m = length(R);
	n = length(Tau);
	U = zeros(n, m);

	if isTerminal
		k = -k;
	end

	% Initial Condition
	U(1, 1:end) = f(R, T);

	for i = 2:n
		A = zeros(m, m);
		b = zeros(m, 1);

		aa = fa(sig, R);
		bb = fb(r, R);
		cc = fc(r);

		A(1:m+1:end) = 2 - 2*aa*k/h^2 + cc*k;
		A(2:m+1:end) = aa(2:m)*k/h^2 - bb(2:m)*k/(2*h);
		A(m+1:m+1:end) = aa(1:m-1)*k/h^2 + bb(1:m-1)*k/(2*h);

		A(1,1) = 1 - (3*k)/(2*h);
		A(1,2) = (4*k)/(2*h);
		A(1,3) = (-k)/(2*h);
		% A(1,1) = 1 - k/h;
		% A(1,2) = k/h;
		A(m,m) = 1;
		A(m,m-1) = 0;

		b(2:m-1) = (-aa(2:m-1)*k/h^2 + bb(2:m-1)*k/(2*h)) .* U(i-1,1:m-2) ...
					+ (2 + 2*aa(2:m-1)*k/h^2 - cc*k) .* U(i-1,2:m-1) ...
					+ (-aa(2:m-1)*k/h^2 - bb(2:m-1)*k/(2*h)) .* U(i-1,3:m);
		b(1) = U(i-1,1);
		b(end) = 0;

		U(i,:) = (A\b)';
	end

	if isTerminal
		U = flipud(U);
	end
end