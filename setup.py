import distutils
from distutils.core import setup, Extension

SWIG_WRAPPER = Extension('dspace.SWIG._dspace_interface',
                         define_macros = [('MAJOR_VERSION', '0'),
                                          ('MINOR_VERSION', '9')],
                         include_dirs = ['/usr/local/include/', '/usr/local/include/designspace/'],
                         libraries = ['designspace'],
                         library_dirs = ['/usr/local/lib'],
                         runtime_library_dirs= ['/usr/local/lib'],
                         extra_compile_args=['-w'],
                         sources = ['dspace/SWIG/designspacetoolbox_wrap.c'])

setup(name = 'dspace',
      version = '1.0',
      description = 'This is a python interface for the design space toolbox',
      author = 'Jason G. Lomnitz',
      author_email = 'jlomn@ucdavis.edu',
      url = 'http://www.bme.ucdavis.edu/savageaulab/',
      long_description = '''
      This is a python interface for the C version of the design space toolbox.
      ''',
      ext_modules = [SWIG_WRAPPER],
      packages=['dspace', 
                'dspace.SWIG', 
                'dspace.models', 
                'dspace.plotutils',
                'dspace.graphs',
                'dspace.display', 
                'dspace.display.UI',
                'dspace.examples'])
