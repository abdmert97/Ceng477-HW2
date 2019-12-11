src = *.cpp

rasterizer_cpp:
	g++ $(src) -std=c++11 -o rasterizer

all:
	g++ $(src) -std=c++11 -o rasterizer
