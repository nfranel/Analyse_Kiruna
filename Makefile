# Variables
CXX = g++
CXXFLAGS = -g3 -O0 -Wall -Wextra -DDEBUG
LDFLAGS = `root-config --libs`
SRC_DIR = src
BUILD_DIR = build
INCLUDE_DIR = include

TARGET1 = FulldatCoinCheck
TARGET2 = thesis_analysis

SRCS1 = $(SRC_DIR)/FulldatCoinCheck.cpp $(SRC_DIR)/CorrectionTimes.cpp
SRCS2 = $(SRC_DIR)/thesis_analysis.cpp $(SRC_DIR)/CorrectionTimes.cpp

OBJS1 = $(SRCS1:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
OBJS2 = $(SRCS2:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Règle principale
all: $(TARGET1) $(TARGET2)

# Compilation des exécutables
$(TARGET1): $(OBJS1)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(TARGET2): $(OBJS2)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Compilation des fichiers objets
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(BUILD_DIR)
	$(CXX) -I$(INCLUDE_DIR) -c $< -o $@ $(CXXFLAGS)

# Nettoyage
clean:
	rm -rf $(BUILD_DIR) $(TARGET1) $(TARGET2)
