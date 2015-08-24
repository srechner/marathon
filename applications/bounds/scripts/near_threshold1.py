import sys

def near_threshold1(n):
	for i in range(n-1):
		print(n-i, end=",")
	print(2, end=";")
	for i in range(n-1):
		print(n-i, end=",")
	print(2)


for i in range(4, int(sys.argv[1])+1):
	near_threshold1(i)
