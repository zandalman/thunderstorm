# compiler and flags
CXX := mpic++
CXXFLAGS_RELEASE := -std=c++17 -Wall -Wextra -O2 -march=native -finline-functions
CXXFLAGS_DEBUG1 := -std=c++17 -Wall -Wextra -O0 -fno-tree-vectorize -g

# Use DEBUG=1 to enable debug flags
ifeq ($(DEBUG),1)
  CXXFLAGS := $(CXXFLAGS_DEBUG1)
else
  CXXFLAGS := $(CXXFLAGS_RELEASE)
endif

# directories
SRC_DIR := ../src
INCLUDE_DIR := ../include
BUILD_DIR := .

# files
SRC := $(wildcard $(SRC_DIR)/*.cpp)
OBJ := $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
INC := -I$(INCLUDE_DIR)
TARGET := thunder

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INC) $^ -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

clean:
	rm -f $(BUILD_DIR)/*.o $(TARGET)

.PHONY: all clean
