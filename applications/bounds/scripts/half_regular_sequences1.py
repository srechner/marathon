import sys

def half_regular(n):
	print(n-1, end=",")
	print(n-2, end=",")
	print(3, end=";")
	for i in range(n-1):
		print(2, end=",")
	print(2)

def half_regular2(n):
	print(n-1, end=",")
	print(n-2, end=",")
	print('2,1', end=";")
	for i in range(n-1):
		print(2, end=",")
	print(2)

for i in range(4, int(sys.argv[1])+1):
	half_regular(i)
