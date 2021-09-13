from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

GENIE3_extension = Extension(
    name="GENIE3",
    sources=["genie3_wrap.pyx"],
    libraries=["genie3"],
    library_dirs=["lib"],
    include_dirs=["lib"]
)
setup(
    name="GENIE3",
    ext_modules = cythonize("genie3_wrap.pyx", language_level=3),
)