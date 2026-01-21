from distutils.core import setup, Extension

ctitest_module = Extension('ctitest',sources = ['ctitest.cpp'])

setup(name = 'CtiTest',version='1.0',description = 'This is a cti test',ext_modules = [ctitest_module])
