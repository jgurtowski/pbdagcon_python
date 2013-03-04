from setuptools import setup, Extension, find_packages
#from distutils.core import setup
from distutils.extension import Extension
from distutils.extension import Extension

# setuptools DWIM monkey-patch madness
# http://mail.python.org/pipermail/distutils-sig/2007-September/thread.html#8204
import sys
if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__


setup(
    setup_requires=['setuptools_cython'],
    name = 'pbtools.pbdagcon',
    version='0.2.0',
    author='pbiDevNet',
    author_email='pbiDevNet@pacificbiosciences.com',
    license='LICENSE.txt',
    packages = find_packages('src'), 
    namespace_packages = ['pbtools'],
    scripts = ["src/gcon.py", "src/hcon.py", "src/q-sense.py"],
    ext_modules = [ Extension("pbtools.pbdagcon.c_aligngraph", ["./src/pbtools/pbdagcon/c_aligngraph.pyx", "./src/pbtools/pbdagcon/c_aligngraph.pxd"]),
                    Extension("pbtools.pbdagcon.c_utils", ["./src/pbtools/pbdagcon/c_utils.pyx"]) ],
    package_dir = {'':'src'},
    zip_safe = False,
    install_requires=[
           "pbcore >= 0.1"
    ]
    )
