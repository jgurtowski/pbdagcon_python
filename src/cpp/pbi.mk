# Defines variables specific to the PBI build env using relative paths

THIRDPARTY ?= ../../../../../ThirdParty
BOOST := $(THIRDPARTY)/boost/boost_1_47_0
LOG4CPP := $(THIRDPARTY)/log4cpp/log4cpp-1.0
BLASR ?= ../../../../lib/cpp/alignment
INCDIRS := -I $(BLASR) -I $(BOOST) -I $(LOG4CPP)/include
LIBDIRS := -L $(BLASR) -L $(BOOST)/stage/linux_gcc/lib -L $(LOG4CPP)/src/.libs
