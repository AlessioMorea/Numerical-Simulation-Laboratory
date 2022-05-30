CC = g++
CFLAGS = -Wall -O3 --std=c++11
main1.exe : main1.o random.o
	$(CC) random.o main1.o -o main1.exe

main1.o : main1.cpp
	$(CC) -c main1.cpp -o main1.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o main1.exe
