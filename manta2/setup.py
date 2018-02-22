from setuptools import setup

setup(
    name='manta2',
    packages=['manta2'],
    include_package_data=True,
    install_requires=[
        'flask',
        'pymongo'
    ],
)
