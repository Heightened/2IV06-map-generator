CXX = $(shell wx-config --cxx)

OBJDIR = obj
SRCDIR = src
BINDIR = bin

PROGRAM = $(BINDIR)/generator

OBJ_NAMES = Window
OBJECTS = $(addsuffix .o, $(addprefix $(OBJDIR)/, $(OBJ_NAMES)))

.SUFFIXES: .o .cpp

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)
	$(CXX) -c `wx-config --cxxflags` -o $@ $<

all: $(PROGRAM)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(PROGRAM): $(OBJECTS) $(BINDIR)
	$(CXX) -o $(PROGRAM) $(OBJECTS) `wx-config --libs`

clean:
	rm -rf $(BINDIR)
	rm -rf $(OBJDIR)
