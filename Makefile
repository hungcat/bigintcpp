### settings for compiling ###
CC      = g++
CFLAGS  = -std=c++11 -g -MMD -MP -O2 -Wall -Wextra -Winit-self -Wno-unused-parameter -Wfloat-equal  
LDFLAGS = 
LIBS    = 
INCLUDE = 

### settings for target files ###
TARGET  = ./bin/$(shell basename `readlink -f .`)
SRC_DIR = ./src
OBJ_DIR = ./obj
SOURCES = $(shell ls $(SRC_DIR)/*.cpp) 
OBJS    = $(subst $(SRC_DIR),$(OBJ_DIR), $(SOURCES:.cpp=.o))
DEPENDS = $(OBJS:.o=.d)


all: $(TARGET)

$(TARGET): $(OBJS) 
	$(CC) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp 
	@if [ ! -d $(OBJ_DIR) ]; then \
		echo "mkdir -p $(OBJ_DIR)"; mkdir -p $(OBJ_DIR); \
	fi
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $< 

clean:
	$(RM) $(OBJS) $(TARGET) $(DEPENDS)

### make run ###
# If the first argument is "run"...
ifeq (run, $(firstword $(MAKECMDGOALS)))
  # use the rest as arguments for "run"
  RUN_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  # ...and turn them into do-nothing targets
  $(eval $(RUN_ARGS):;@:)
endif

run: 
	./$(TARGET) $(RUN_ARGS)

-include $(DEPENDS)

# phony
.PHONY: all clean run $(RUN_ARGS)
