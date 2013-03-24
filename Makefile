LDFLAGS += -Wl,--no-as-needed $(shell root-config --libs) -ggdb3
CXXFLAGS += $(shell root-config --cflags) -ggdb3 -I.
CXXFLAGS += -Wextra

DEP_OBJ=CastorTreeVariables.o CastorTreeVariablesDict.o ParticleInfo.o
EXE=csana2

all: $(EXE)

.PHONY: all clean

CastorTreeVariablesDict.cc: CastorTreeVariables.h CastorTreeVariablesLinkDef.h
	rootcint -f $@  -c  $^

$(EXE): $(DEP_OBJ)

clean:
	rm -f *~ *.o *.so $(EXE) *.d
