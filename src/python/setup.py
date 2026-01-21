from distutils.core import setup, Extension

#cti_module = Extension('cti',sources = ['ctimodule.cpp'])
cti_module = Extension('cti',sources = [])

setup(name = 'CtiInterface',version='1.0',description = 'Provides access to cti classes and functions',ext_modules = [cti_module])
