LDFLAGS += -Wl,--no-as-needed $(shell root-config --libs) -ggdb3
CXXFLAGS += $(shell root-config --cflags) -ggdb3 -I.
CXXFLAGS += -Wextra
CXXFLAGS += -std=c++0x

DEP_OBJ=CastorTreeVariables.o CastorTreeVariablesDict.o ParticleInfo.o
EXE=treeana

all: $(EXE)

.PHONY: all clean

CastorTreeVariablesDict.cc: CastorTreeVariables.h CastorTreeVariablesLinkDef.h
	rootcint -f $@  -c  $^

$(EXE): $(DEP_OBJ)

clean:
	rm -f *~ *.o *.so $(EXE) *.d
