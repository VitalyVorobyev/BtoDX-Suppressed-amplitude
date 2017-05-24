SourceDir=src/
ObjectDir=obj/
BinDir=bin/
Sources=$(notdir $(wildcard $(SourceDir)*.cpp))
Executable=cpvfit
CFlags=-c -Wall -g -Iinc -I. `root-config --cflags` -std=c++14
LDFlags= -I. -Wl,--no-as-needed `root-config --glibs` -lm -lstdc++ -ltatami -lfaddeeva -lMinuit2

CC=g++ -O2
RM=rm

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
Objects=$(Sources:.cpp=.o)
CSources=$(addprefix $(SourceDir),$(Sources))
CObjects=$(addprefix $(ObjectDir),$(Objects))
CExecutable=$(addprefix $(BinDir),$(Executable))

all: $(CSources) $(CExecutable)

$(CExecutable): $(CObjects)
	$(CC) $(LDFlags) $(CObjects) -o $@

$(ObjectDir)%.o: $(SourceDir)%.cpp
	$(CC) $(CFlags) $< -o $@

clean:
	$(RM) $(CObjects)
