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
	
