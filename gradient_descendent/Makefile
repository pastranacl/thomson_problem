# Cesar L. Pastrana, 2021

CC := gcc
FLAGS := -std=c99 -lm --fast-math -fopenmp -O3

BIN_PATH := ./bin
SRC_PATH := ./src

SRCS := ${SRC_PATH}/thomson.c ${SRC_PATH}/rands.c


comp:
	@clear
	@echo -n "Compiling... "
	@${CC} ${SRCS} ${FLAGS} -o ${BIN_PATH}/mcsphere
	@echo "Done!"


