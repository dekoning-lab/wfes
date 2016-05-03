from distutils.core import setup, Extension
from Cython.Build import cythonize
from glob import glob

extra_libs = ["-lmkl_intel_ilp64", "-lmkl_core", "-lmkl_intel_thread",
              "-lpthread", "-lmkl_avx", "-lmkl_vml_avx", "-lm", "-liomp5"]

extensions = [Extension("wfes",
                        sources=["wfes.pyx"] + glob("src/*.c"),
                        include_dirs=["include/"],
                        extra_compile_args=["-Wno-unused"],
                        extra_link_args=extra_libs)]

setup(
    ext_modules=cythonize(extensions)
)
