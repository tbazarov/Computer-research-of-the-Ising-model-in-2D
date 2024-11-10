g++ -std=c++17 randomc.cpp -c -W -Wall -O3 -o randomc.o
g++ -std=c++17 main.cpp -c -W -Wall -O3 -o main.o
g++ *.o -o start
./start
