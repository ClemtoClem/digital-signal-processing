PROJECT_NAME := signaux
CC        := g++
SRCDIR    := src
HEADERDIR := src
BUILDDIR  := build
BINDIR    := bin
TARGET    := $(BINDIR)/$(PROJECT_NAME)
SOURCES   := $(wildcard $(SRCDIR)/*.c $(SRCDIR)/*.cpp)
HEDEARS   := $(wildcard $(HEADERDIR)/*.h $(HEADERDIR)/*.hpp)
OBJECTS   := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(addsuffix .o,$(basename $(SOURCES))))
DEPS      := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(addsuffix .d,$(basename $(SOURCES))))
CFLAGS    := -Wall -g
LIB       := -lm
INC       := -I $(SRCDIR) # -I include

# List of all files in the directory
ALL_FILES := $(wildcard $(BINDIR)/*) $(wildcard $(BUILDDIR)/*)
# List of .dll files
DLL_FILES := $(wildcard $(BINDIR)/*.dll)
# List of files to be deleted (excluding .dll files)
DELETE_FILES := $(filter-out $(DLL_FILES), $(ALL_FILES))

GREEN=`tput setaf 2`
RESET=`tput sgr0`

define print_green
	@echo "$(GREEN)$(1)$(RESET)"
endef

all: $(TARGET)

clean:
	rm -r $(DELETE_FILES)

$(TARGET): $(BINDIR) $(BUILDDIR) $(OBJECTS)
	$(call print_green,"Linking object files...")
	@$(CC) $(OBJECTS) -o $(TARGET) $(LIB)
	$(call print_green,"$(TARGET) has been created!")

$(BUILDDIR) :
	mkdir $(BUILDDIR)

$(BINDIR):
	mkdir $(BINDIR)
	
$(BUILDDIR)/%.o: $(SRCDIR)/%.c*
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INC) -M $< -MT $@ > $(@:.o=.td)
	@cp $(@:.o=.td) $(@:.o=.d); 
	@sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	-e '/^$$/ d' -e 's/$$/ :/' < $(@:.o=.td) >> $(@:.o=.d); 
	@rm -f $(@:.o=.td)


-include $(DEPS)

.PHONY: clean all
