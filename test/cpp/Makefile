# pbdagcon source code
SRCDIR := ../../src/cpp

# contains some PBI-specific build environment stuff
include $(SRCDIR)/pbi.mk

# Flags passed to the preprocessor.
CPPFLAGS += $(INCDIRS) -I$(SRCDIR) -I$(GTEST)/include

# Flags passed to the C++ compiler.
CXXFLAGS += -g -std=c++0x -fpermissive

# All tests produced by this Makefile.
TESTS = test_alngraph \
		test_alignment \
		test_simple_aligner

# Source files needed for tests
SRCS := $(SRCDIR)/AlnGraphBoost.cpp \
		AlnGraphBoostTest.cpp \
		$(SRCDIR)/Alignment.cpp \
		AlignmentTest.cpp \
		$(SRCDIR)/SimpleAligner.cpp \
		SimpleAlignerTest.cpp

OBJS := $(SRCS:.cpp=.o)
LIBS := -lblasr -llog4cpp -lpthread 

all : check

check : $(TESTS)
	./$<

clean :
	rm -f $(TESTS) gtest.a gtest_main.a $(OBJS) 

%.o : %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c "$<" -o "$@"

test_% : $(OBJS) gtest_main.a
	$(CXX) $(LIBDIRS) -static $^ -o $@ $(LIBS)

# The rest is building gtest
GTEST_SRCS_ = $(GTEST)/src/*.cc $(GTEST)/src/*.h 

gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(GTEST)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^
