# Variables
CXX = g++
#CXXFLAGS = `root-config --cflags` -Wall -O2
CXXFLAGS = -g3 -O0 -Wall -Wextra -DDEBUG
LDFLAGS = `root-config --libs`
SRC_DIR = src
BUILD_DIR = build
INCLUDE_DIR = include
TARGET = FulldatCoinCheck

# Fichiers source et objets
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Règle principale
all: $(TARGET)

# Compilation du programme principal
$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compilation des fichiers objets
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) -I$(INCLUDE_DIR) -c $< -o $@ $(CXXFLAGS)

# Nettoyage
clean:
	rm -rf $(BUILD_DIR) $(TARGET)
