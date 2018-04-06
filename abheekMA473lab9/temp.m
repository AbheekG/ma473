function temp
	close all; clear; im_num = 1;
	x_min = 0;
	x_max = 1;
	N = 10;
	h = (x_max - x_min)/N;
	X = x_min:h:x_max;


	deg = 2;
	i = 6;
	X =  x_min:0.01:x_max;
	Y = [];
	for j = 1:length(X)
		Y = [Y, dphi2(x_min, x_max, N, i, X(j), 0)];
	end
	plot(X, Y)

	% dphi1(x_min, x_max, N, 5, 0.5, 0)
	% dphi1(x_min, x_max, N, 5, 0.5, 1)
end