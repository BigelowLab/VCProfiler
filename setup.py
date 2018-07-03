import io
from setuptools import setup

setup(
    name='vcprofiler',
    version=0.0,
    url='',
    license='',
    author='Julia Brown',
    author_email='julia@bigelow.org',
    description='',
    long_description=__doc__,
    py_modules=['vcprofiler'],
    packages=['vcprofiler'],
    install_requires=[],
    entry_points='''
    [console_scripts]
    vcprofiler=vcprofiler.__main__:cli
    '''
)
