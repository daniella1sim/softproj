from setuptools import setup, Extension

module = Extension("symnmf", sources=["symnmf.c", "symnmfmodule.c"])

setup(
    name="symnmf",
    version="1.0",
    description="Python interface for the symnmf C extension",
    ext_modules=[module],
)
