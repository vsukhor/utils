# NOTE: You may need to supply executables for $(AR), $(AS) and $(CC)
# to override or supplement the environment variables
# e.g., the call from the directory where the Makefile is, may look like:
# make all AS=/usr/bin/as AR=/usr/bin/ar CC=/usr/local/gcc-9.1.0/bin/gcc-9.1.0
# Alternatively, $(AR), $(AS) and $(CC) may be hardcoded directly in the local copy of this file.
#
# To check the actual value of a variable $(VAR):
# make print-VAR

print-%  : ; @echo $* = $($*)

SRC_DIR := ../utils/
ARRAYS_DIR := $(SRC_DIR)arrays/
BIOCHEMICAL_DIR := $(SRC_DIR)biochemical/
COMMON_DIR := $(SRC_DIR)common/
RANDOM_DIR := $(SRC_DIR)random/

SRCS_BIOCHEMICAL := pdb.cpp
SRCS_COMMON := threads.cpp misc.cpp msgr.cpp exceptions.cpp

SRCS_BIOCHEMICAL_ALL := $(patsubst %,$(BIOCHEMICAL_DIR)%,$(SRCS_BIOCHEMICAL))
SRCS_COMMON_ALL := $(patsubst %,$(COMMON_DIR)%,$(SRCS_COMMON))

SRCS_ALL += $(SRCS_BIOCHEMICAL_ALL)
SRCS_ALL += $(SRCS_COMMON_ALL)

ASMS_DIR := ./asm/
ASMS += $(SRCS_BIOCHEMICAL:.cpp=.s)
ASMS += $(SRCS_COMMON:.cpp=.s)
ASMS_ALL := $(patsubst %,$(ASMS_DIR)%,$(ASMS))

.PHONY: all

all : $(ASMS_ALL)
$(ASMS_ALL) : | $(ASMS_DIR)
$(ASMS_DIR) :
	mkdir -p $(ASMS_DIR)

OBJS_DIR := ./obj/
OBJS += $(SRCS_BIOCHEMICAL:.cpp=.o)
OBJS += $(SRCS_COMMON:.cpp=.o)
OBJS_ALL := $(patsubst %,$(OBJS_DIR)%,$(OBJS))

all : $(OBJS_ALL)
$(OBJS_ALL) : | $(OBJ_DIR)
$(OBJ_DIR) :
	mkdir -p $(OBJ_DIR)

TARGET_DIR := ../lib/
TARGET_NAME := libutils.a
TARGET := $(TARGET_DIR)$(TARGET_NAME)
all : $(TARGET)

INCDIR := -I$(BOOST_DIR)

CC = gcc
CFLAGS := -c -Wall -O3 -msse4 -march=native -std=c++17 -fPIC

$(TARGET) : $(OBJS_ALL)
	mkdir -p $(TARGET_DIR)
	$(AR) -rcs $@ $(OBJS_ALL)

$(OBJS_DIR)%.o : $(ASMS_DIR)%.s
	$(AS) $< -o $@

$(ASMS_DIR)misc.s : $(COMMON_DIR)misc.cpp \
	$(ARRAYS_DIR)array2.h \
	$(ARRAYS_DIR)array3.h \
	$(ARRAYS_DIR)array4.h \
	$(COMMON_DIR)misc.h
	$(CC) $(CFLAGS) -S $< -o $@

$(ASMS_DIR)threads.s : $(COMMON_DIR)threads.cpp \
	$(COMMON_DIR)misc.h \
	$(COMMON_DIR)threads.h
	$(CC) $(CFLAGS) -S $< -o $@
	
$(ASMS_DIR)msgr.s : $(COMMON_DIR)msgr.cpp \
	$(COMMON_DIR)misc.h \
	$(COMMON_DIR)msgr.h
	$(CC) $(CFLAGS) -S $< -o $@

$(ASMS_DIR)pdb.s : $(BIOCHEMICAL_DIR)pdb.cpp \
	$(ARRAYS_DIR)array3.h \
	$(COMMON_DIR)msgr.h
	$(CC) $(CFLAGS) -S $< -o $@

$(ASMS_DIR)exceptions.s : $(COMMON_DIR)exceptions.cpp \
	$(COMMON_DIR)msgr.h \
	$(COMMON_DIR)misc.h \
	$(COMMON_DIR)exceptions.h
	$(CC) $(CFLAGS) -S $< -o $@

#clean:

.PHONY: clean

clean:
	rm -rf $(OBJS_DIR)*.o $(ASMS_DIR)*.s $(TARGET)
