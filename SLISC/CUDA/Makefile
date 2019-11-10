goal:
	nvcc -c --std=c++11 cuBLAS_wrapper.cpp
	g++ -c --std=c++17 cuBLAS_test.cpp
	nvcc -o main.x cuBLAS_test.o cuBLAS_wrapper.o -lcublas
clean:
	*.x *.o
