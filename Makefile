SHELL = /bin/bash -e
LIBS =  pbdata/libpbdata.a alignment/libblasr.a

all : cpp
clean : cpp-clean
test : cpp-check

# C++ project build directives
cpp : $(LIBS)
	$(MAKE) -C src/cpp

cpp-check : cpp
	$(MAKE) -C test/cpp

cpp-clean :
	$(MAKE) -C src/cpp clean
	$(MAKE) -C test/cpp clean
	$(MAKE) -C pbdata clean
	$(MAKE) -C alignment clean

$(LIBS) : 
	$(MAKE) -C $(@D) nohdf=1

