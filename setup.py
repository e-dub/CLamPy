# -*- coding: utf-8 -*-

import os
def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

if __name__ == '__main__':
    from distutils.core import setup
    extra_files = package_files('LamPy')
    setup(name='LamPy',
          version='1.0',
          description='classical LAMinate theory for the lightweight design of structures and systems in PYthon',
          author='E. M. Gioia and E. J. Wehrle',
          package_data={'': extra_files},
          license='GNU Lesser General Public License',
          packages=['LamPy'])
