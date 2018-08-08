#
CC = g++
GUROBI_HOME=/home/vanduc/Downloads/gurobi801/linux64
CFLAGS = -std=c++11 -m64 -g
INC      = $(GUROBI_HOME)/include/
CPPLIB   = -L${GUROBI_HOME}/lib -lgurobi_c++ -lgurobi80
TARGET = DecisionEngine_test
SRC_DIR = src
# Sources files
SRCS = ${SRC_DIR}/DecisionEngine_test.cpp ${SRC_DIR}/DecisionEngine.cpp ${SRC_DIR}/Metadata.cpp ${SRC_DIR}/Projection.cpp ${SRC_DIR}/common.cpp ${SRC_DIR}/GaussianFilter.cpp
# Object files
OBJS = $(SRCS:.c=.o)
# Build
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) -I$(INC) $(CPPLIB) -lm

