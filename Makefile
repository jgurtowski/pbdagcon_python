.PHONY: all clean install dev-install test
SHELL = /bin/bash -e

all: install

install:
	python setup.py install

dev-install:
	python setup.py develop

clean:
	rm -rf build/;\
	find . -name "*.egg-info" | xargs rm -rf;\
	rm -rf dist/

test:
	nosetests src/tests/*
	cram examples/*.t
	
# C++ project build directives
cpp:
	$(MAKE) -C src/cpp

cpp-check: cpp
	$(MAKE) -C test/cpp

cpp-clean:
	$(MAKE) -C src/cpp clean
	$(MAKE) -C test/cpp clean
