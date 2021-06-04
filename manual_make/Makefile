# To build the library, from the directory where this Makefile is, call:
#    make all
#
# NOTE: You may need to supply executables for $(AR) and $(CC)
# to override or supplement the environment variables.
#
# To clean the previous build, call:
#    make clean
#
# To check the actual value of a variable $(VAR):
#    make print-VAR
print-%  : ; @echo $* = $($*)

SRC_DIR := ../utils/
ARRAYS_DIR := $(SRC_DIR)arrays/
BIOCHEMICAL_DIR := $(SRC_DIR)biochemical/
COMMON_DIR := $(SRC_DIR)common/

SRCS_BIOCHEMICAL := pdb.cpp
SRCS_COMMON := threads.cpp msgr.cpp exceptions.cpp

SRCS_BIOCHEMICAL_ALL := $(patsubst %,$(BIOCHEMICAL_DIR)%,$(SRCS_BIOCHEMICAL))
SRCS_COMMON_ALL := $(patsubst %,$(COMMON_DIR)%,$(SRCS_COMMON))

SRCS_ALL += $(SRCS_BIOCHEMICAL_ALL)
SRCS_ALL += $(SRCS_COMMON_ALL)

.PHONY: all

OBJS_DIR := ./obj/
OBJS += $(SRCS_BIOCHEMICAL:.cpp=.o)
OBJS += $(SRCS_COMMON:.cpp=.o)
OBJS_ALL := $(patsubst %,$(OBJS_DIR)%,$(OBJS))

all : $(OBJS_ALL)
$(OBJS_ALL) : | $(OBJS_DIR)
$(OBJS_DIR) :
	mkdir -p $(OBJS_DIR)

TARGET_DIR := ../lib/
TARGET_NAME := libutils.a
TARGET := $(TARGET_DIR)$(TARGET_NAME)
all : $(TARGET)

CC = g++
CFLAGS := -c -Wall -O3 -msse4 -march=native -std=c++20 -fPIC

$(TARGET) : $(OBJS_ALL)
	mkdir -p $(TARGET_DIR)
	$(AR) -rcs $@ $(OBJS_ALL) && \
    ranlib $(TARGET)

$(OBJS_DIR)threads.o : $(COMMON_DIR)threads.cpp \
	$(COMMON_DIR)constants.h \
	$(COMMON_DIR)misc.h \
	$(COMMON_DIR)msgr.h \
	$(COMMON_DIR)threads.h
	$(CC) $(CFLAGS) $< -o $@
	
$(OBJS_DIR)msgr.o : $(COMMON_DIR)msgr.cpp \
	$(COMMON_DIR)constants.h \
	$(COMMON_DIR)misc.h \
	$(COMMON_DIR)msgr.h
	$(CC) $(CFLAGS) $< -o $@

$(OBJS_DIR)pdb.o : $(BIOCHEMICAL_DIR)pdb.cpp \
	$(ARRAYS_DIR)array3.h \
	$(COMMON_DIR)msgr.h
	$(CC) $(CFLAGS) $< -o $@

$(OBJS_DIR)exceptions.o : $(COMMON_DIR)exceptions.cpp \
	$(COMMON_DIR)constants.h \
	$(COMMON_DIR)misc.h \
	$(COMMON_DIR)msgr.h \
	$(COMMON_DIR)exceptions.h
	$(CC) $(CFLAGS) $< -o $@

# Cleaning:
.PHONY: clean
clean:
	rm -rf $(OBJS_DIR)*.o $(TARGET)
