all:
	g++ -std=gnu++2a -fopenmp -Wall -Werror -pedantic-errors -Ofast -march=native -shared -fPIC -I/usr/include -I/usr/include/pybind11 `python3 -m pybind11 --includes` pymodule.cpp -o smic`python3-config --extension-suffix`
