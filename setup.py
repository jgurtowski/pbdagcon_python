from setuptools import setup, Extension, find_packages

setup(
    name = 'pbtools.pbdagcon',
    version='0.1.0',
    author='pbiDevNet',
    author_email='pbiDevNet@pacificbiosciences.com',
    license='LICENSE.txt',
    packages = find_packages('src'), 
    namespace_packages = ['pbtools'],
    scripts = ["src/gcon.py", "src/hcon.py"],
    package_dir = {'':'src'},
    zip_safe = False,
    install_requires=[
           "pbcore >= 0.1"
    ]
    )
