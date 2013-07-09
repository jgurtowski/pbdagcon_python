THIRDPARTY := ../../../../../ThirdParty
BOOST := $(THIRDPARTY)/boost/boost_1_47_0
BOOST_THREAD_ARCH := $(BOOST)/stage/linux_gcc/lib/libboost_thread.a
BOOST_OPTS_ARCH := $(BOOST)/stage/linux_gcc/lib/libboost_program_options.a
LOG4CPP := $(THIRDPARTY)/log4cpp/log4cpp-1.0
BLASR_COMMON := ../../../../../private/blasr/src
BLASR_ARCH := ../../../../../private/blasr/src/libblasr.a
CXX_OPTS := -I $(BLASR_COMMON) -I $(LOG4CPP)/include -I $(BOOST) -std=c++0x -g -Wall -Wuninitialized -Wno-div-by-zero -fpermissive -pedantic -c -fmessage-length=0 -MMD -MP
#CXX_OPTS := -I $(LOG4CPP)/include -I $(BOOST) -std=c++0x -O2 -c -fmessage-length=0 -MMD -MP

all: pbdagcon 

pbdagcon: main.o AlnGraphBoost.o Alignment.o BlasrM5AlnProvider.o SimpleAligner.o $(BLASR_ARCH) $(LOG4CPP)/src/.libs/liblog4cpp.a $(BOOST_THREAD_ARCH) $(BOOST_OPTS_ARCH)
	$(CXX) -o $@ $^ -lpthread

%.o: %.cpp | $(LOG4CPP)/include/log4cpp/config.h
	$(CXX) $(CXX_OPTS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"

$(LOG4CPP)/src/.libs/liblog4cpp.a:
	$(MAKE) -C $(LOG4CPP)/src

$(LOG4CPP)/include/log4cpp/config.h:
	cd $(LOG4CPP) && ./configure

clean:
	rm -f *.d
	rm -f *.o
	rm -f pbdagcon

.PHONY: all clean
