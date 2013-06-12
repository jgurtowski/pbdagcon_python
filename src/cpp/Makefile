THIRDPARTY := ../../../../../ThirdParty
BOOST := $(THIRDPARTY)/boost/boost_1_47_0
LOG4CPP := $(THIRDPARTY)/log4cpp/log4cpp-1.0
CXX_OPTS := -I $(LOG4CPP)/include -I $(BOOST) -std=c++0x -g -Wall -Wuninitialized -pedantic -c -fmessage-length=0 -MMD -MP
#CXX_OPTS := $(BOOST) -std=c++0x -O2 -c -fmessage-length=0 -MMD -MP

all: gcon 

gcon: main.o AlnGraphBoost.o Alignment.o $(LOG4CPP)/src/.libs/liblog4cpp.a
	$(CXX) -o gcon $^ -lpthread

%.o: %.cpp
	$(CXX) $(CXX_OPTS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"

$(LOG4CPP)/src/.libs/liblog4cpp.a: $(LOG4CPP)/config.log
	$(MAKE) -C $(LOG4CPP)

$(LOG4CPP)/config.log:
	cd $(LOG4CPP) && ./configure --with-pthreads

clean:
	rm -f *.d
	rm -f *.o
	rm -f gcon

.PHONY: all clean
