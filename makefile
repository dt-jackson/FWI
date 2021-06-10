TARGET := xfwi

CC := gcc
CFLAGS := -Wall -Wextra -Wpedantic -O1 -g

# path to include files (do not include -I)
INCDIR := src
# path to shared libraries (do not include -L)
LIBDIR := src


# no need to include the "-l" or "lib" at the beginning
LIBS := fftw3 qhull_r crypto lapacke cblas m


SRCDIR := src
BUILDDIR := obj
TARGETDIR := bin


SRCEXT := c
OBJEXT := o
DEPEXT := d

#================================================================
SOURCES := $(shell find $(SRCDIR) -maxdepth 1 -type f -name "*.$(SRCEXT)")
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

LIBS :=   $(addprefix -l, $(LIBS))
LIBDIR := $(addprefix -L, $(LIBDIR))
INCDIR := $(addprefix -I, $(INCDIR))



$(TARGET) : $(TARGETDIR) $(BUILDDIR) $(TARGETDIR)/$(TARGET)

#link
$(TARGETDIR)/$(TARGET) : $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $(LIBDIR) $(OBJECTS) $(LIBS)

# include dependencies
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT));

# compile
$(BUILDDIR)/%.$(OBJEXT) : $(SRCDIR)/%.$(SRCEXT)
	$(CC) $(CFLAGS) $(INCDIR) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDIR) -MM -MT $(BUILDDIR)/$*.$(OBJEXT) $< > $(BUILDDIR)/$*.$(DEPEXT)


$(TARGETDIR) : 
	mkdir $(TARGETDIR)

$(BUILDDIR) :
	mkdir $(BUILDDIR)

clean :
	@rm -f $(TARGETDIR)/$(TARGET)

cleaner : clean
	@rm -f $(OBJECTS) $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

.PHONY : clean cleaner $(TARGET)

