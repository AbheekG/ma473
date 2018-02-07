A = [29, 29];
a = '''\includegraphics{"q%d_%d"}
\pagebreak

''';

for i in range(len(A)):
	for j in range(A[i]):
		print(a % (i+1, j+1));
