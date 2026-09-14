CC ?= cc
CFLAGS ?= -O2 -std=c99 -Wall -Wextra -pedantic

TARGET := kinetic_fokker_planck_solver
SOURCE := kinetic_fokker_planck_solver.c

UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
OPENBLAS_PREFIX := $(shell brew --prefix openblas 2>/dev/null)
ifeq ($(OPENBLAS_PREFIX),)
$(error OpenBLAS was not found. Install it with: brew install openblas)
endif
CPPFLAGS += -I$(OPENBLAS_PREFIX)/include
LDFLAGS += -L$(OPENBLAS_PREFIX)/lib
LDLIBS += -lopenblas
else
LDLIBS += -llapacke -llapack -lblas -lm
endif

.PHONY: all run clean

all: $(TARGET)

$(TARGET): $(SOURCE)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
