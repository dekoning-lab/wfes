from distutils.core import setup, Extension
from Cython.Build import cythonize
from glob import glob
from platform import system

extra_libs = ["-lmkl_intel_ilp64", "-lmkl_core", "-lmkl_intel_thread",
              "-lpthread", "-lmkl_avx", "-lmkl_vml_avx", "-liomp5",
              "-lm", "-ldl"]

extra_flags = ["-std=c99", "-DMKL_ILP64", "-m64", "-Wno-unused-variable"]

if system() == "Linux":
    extra_flags += ["-fopenmp"]

# This assumes that the compiler is CLANG
elif system() == "Darwin":
    extra_flags += ["-Wno-unused-function",
                    "-Wno-unneeded-internal-declaration"]


extensions = [Extension("wfes",
                        sources=["wfes.pyx"] + glob("src/*.c"),
                        include_dirs=["include/"],
                        extra_compile_args=extra_flags,
                        extra_link_args=extra_libs)]

setup(
    ext_modules=cythonize(extensions)
)
