from os.path import abspath, dirname, join
from setuptools import setup
from setuptools.command import develop, build_py
import subprocess

here = abspath(dirname(__file__))

package = 'rrn'
cmd_build = 'make'
shared_libs = [join(here, package, 'visibility', 'visibilitylib.so'),
               join(here, package, 'corrsift', 'libcorrsift.so')]

about = {}
with open(join(here, package, '__about__.py'), 'r', encoding='utf-8') as f:
    exec(f.read(), about)

with open(join(here, 'requirements.txt'), 'r', encoding='utf-8') as f:
    requirements = [r.rstrip() for r in f]


def readme():
    with open(join(here, 'README.md'), 'r', encoding='utf-8') as f:
        return f.read()


def run_with_verbose_failure(cmd):
    try:
        p = subprocess.run(cmd_build, shell=True, check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           encoding='utf-8')
    except subprocess.CalledProcessError as e:
        print('----------------------------------------\n'
              f'\nERROR running command "{cmd}". (exit code: {e.returncode})\n'
              f'\nstdout:\n\n{e.stdout}\nstderr:\n\n{e.stderr}')
        raise


class CustomDevelop(develop.develop, object):
    """
    Class needed for "pip install -e ."
    """
    def run(self):
        run_with_verbose_failure(cmd_build)
        super(CustomDevelop, self).run()


class CustomBuildPy(build_py.build_py, object):
    """
    Class needed for "pip install <package>"
    """
    def run(self):
        libdir = join(self.build_lib, package)
        self.mkpath(libdir)
        run_with_verbose_failure(cmd_build)
        run_with_verbose_failure(f'cp -r {" ".join(shared_libs)} {libdir}')
        super(CustomBuildPy, self).run()


setup(name=about['__title__'],
      version=about['__version__'],
      description=about['__description__'],
      long_description=readme(),
      long_description_content_type='text/markdown',
      url=about['__url__'],
      author=about['__author__'],
      author_email=about['__author_email__'],
      packages=[package],
      install_requires=requirements,
      python_requires='>=3.5',
      cmdclass={'develop': CustomDevelop,
                'build_py': CustomBuildPy},
      entry_points={'console_scripts':['rrn = rrn.normalization:cli']})

