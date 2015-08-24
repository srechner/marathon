import sys

def half_regular3(n):
	print(n-2, end=",")
	print(n-3, end=",")
	print(5, end=";")
	for i in range(n-1):
		print(2, end=",")
	print(2)


for i in range(4, int(sys.argv[1])+1):
	half_regular3(i)
