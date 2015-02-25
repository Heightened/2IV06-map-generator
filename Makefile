CXX = $(shell wx-config --cxx)

OBJDIR = obj
OBJVENDORDIR = obj/vendor
SRCDIR = src
SRCVENDORDIR = src/vendor
BINDIR = bin

WX_FLAGS = `wx-config --cxxflags`

WX_LIBS = `wx-config --libs --gl-libs`
LIBS = $(WX_LIBS) -lGL -lGLU -lglew

PROGRAM = $(BINDIR)/generator

OBJ_NAMES = Canvas Window Generator
OBJECTS = $(addsuffix .o, $(addprefix $(OBJDIR)/, $(OBJ_NAMES)))

OBJ_VENDOR = VoronoiDiagramGenerator
OBJECTSVENDOR = $(addsuffix .o, $(addprefix $(OBJVENDORDIR)/, $(OBJ_VENDOR)))

.SUFFIXES: .o .cpp

all: $(PROGRAM)

run: $(PROGRAM)
	./$(PROGRAM)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJVENDORDIR):
	mkdir -p $(OBJVENDORDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)
	$(CXX) -c $(WX_FLAGS) -o $@ $<

$(OBJVENDORDIR)/%.o: $(SRCVENDORDIR)/%.cpp $(OBJVENDORDIR)
	$(CXX) -c $(WX_FLAGS) -o $@ $<

$(PROGRAM): $(OBJECTS) $(OBJECTSVENDOR) $(BINDIR)
	$(CXX) -o $(PROGRAM) $(OBJECTS) $(OBJECTSVENDOR) $(LIBS)

clean:
	rm -rf $(BINDIR)
	rm -rf $(OBJDIR)
