# Cesar L. Pastrana, 2021

CC := gcc
FLAGS := -std=c99 -lm --fast-math -fopenmp -O3

BIN_PATH := ./bin
SRC_PATH := ./src

SRCS := $(wildcard ${SRC_PATH}/*.c)


comp:
	@clear
	@echo -n "Compiling... "
	@${CC} ${SRCS} ${FLAGS} -o ${BIN_PATH}/mcsphere
	@echo "Done!"


