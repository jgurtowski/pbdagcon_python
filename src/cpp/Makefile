BOOST := -I ../../../../../ThirdParty/boost/boost_1_47_0
CXX_OPTS := $(BOOST) -std=c++0x -g -Wall -Wuninitialized -pedantic -c -fmessage-length=0 -MMD -MP
#CXX_OPTS := $(BOOST) -std=c++0x -O2 -c -fmessage-length=0 -MMD -MP

all: gcon 

gcon: main.o AlnGraphBoost.o Alignment.o
	$(CXX) -o gcon $^

%.o: %.cpp
	$(CXX) $(CXX_OPTS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"

clean:
	rm -f *.d
	rm -f *.o
	rm -f gcon

.PHONY: all clean
