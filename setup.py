from setuptools import setup
import sys

if sys.version_info < (3, 9):
    sys.exit('Sorry, Python < 3.9 is not supported')

# todo: requirements
setup(
    name='pytomebio',
    version='0.0.1',
    author='Fulcrum Genomics',
    author_email='no-reply@fulcrumgenomics.com',
    maintainer='Fulcrum Genomics',
    maintainer_email='no-reply@fulcrumgenomics.com',
    description='Tools and APIs for off-target detection assays',
    url='https://github.com/tomebio/tb-bfx',
    packages=['pytomebio'],
    package_dir={'': 'src/python'},
    entry_points={
        'console_scripts': ['tomebio-tools=pytomebio.tools.__main__:main']
    },
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        'Environment :: Console',
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
    ]
)
