CXX = $(shell wx-config --cxx)

OBJDIR = obj
SRCDIR = src
BINDIR = bin

WX_FLAGS = `wx-config --cxxflags`

WX_LIBS = `wx-config --libs --gl-libs`
LIBS = $(WX_LIBS) -lGL -lGLU -lglew

PROGRAM = $(BINDIR)/generator

OBJ_NAMES = Canvas Window
OBJECTS = $(addsuffix .o, $(addprefix $(OBJDIR)/, $(OBJ_NAMES)))

.SUFFIXES: .o .cpp

all: $(PROGRAM)

run: $(PROGRAM)
	./$(PROGRAM)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)
	$(CXX) -c $(WX_FLAGS) -o $@ $<

$(PROGRAM): $(OBJECTS) $(BINDIR)
	$(CXX) -o $(PROGRAM) $(OBJECTS) $(LIBS)

clean:
	rm -rf $(BINDIR)
	rm -rf $(OBJDIR)
