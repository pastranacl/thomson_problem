# Cesar L. Pastrana, 2022

CC := g++
FLAGS := -std=c++11 -lm --fast-math -fopenmp -O3
LIBS := -I /usr/include/eigen3
BIN_PATH := ./bin
SRC_PATH := ./src

SRCS := $(wildcard ${SRC_PATH}/*.cpp)


comp:
	@clear
	@echo -n "Compiling... "
	@${CC} ${SRCS} ${FLAGS} ${LIBS} -o ${BIN_PATH}/thomson
	@echo "Done!"


