from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy, os

#bourbon_dir = os.path.abspath(os.path.dirname(bourbon.__file__))
#zen_dir = os.path.abspath(os.path.dirname(zen.__file__))
np_dir = os.path.abspath(os.path.dirname(__file__))
include_dirs = [numpy.get_include(),".",np_dir,"../src"] #, os.path.abspath(os.path.join(bourbon_dir, '..')), os.path.abspath(os.path.join(zen_dir, '..'))]

def gppg(path):
	return os.path.join('../src', path)
	
modules = [
		Extension('gppg.model.sequence', 
					['gppg/model/sequence.pyx', gppg('Model/Sequence/Data.cpp'), gppg('Model/Sequence/Operation.cpp')], 
					include_dirs=include_dirs, language='c++'),
#		Extension('netpop.util.nobject', ['netpop/util/nobject.pyx'], include_dirs=include_dirs),
]

setup(
	name = 'GPPG Library',
	version = '0.1',
	author = 'Troy Ruths',
	author_email= 'troyruths@rice.edu',
	ext_modules = modules,
	packages = ['gppg', 'gppg.model'],
	package_data = {	'gppg' : ['*.pxd'],
						'gppg.model' : ['*.pxd']
					},
	cmdclass = {'build_ext': build_ext}
)


