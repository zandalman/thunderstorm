# compiler and flags
CXX := mpic++
CXXFLAGS := -std=c++17 -Wall -Wextra -O3 -march=native -finline-functions

# directories
SRC_DIR := ../src
INCLUDE_DIR := ../include
BUILD_DIR := .

# files
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
INC := -I$(INCLUDE_DIR)
TARGET := thunderstorm

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INC) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

clean:
	rm -f $(BUILD_DIR)/*.o $(TARGET)

.PHONY: all clean
