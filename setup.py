from setuptools import setup, find_packages

setup(
    name='petrosim',
    version='0.0.2',
    description='A collection of petrological simulation tools.',
    packages=[
        'petrosim',
        'petrosim.models',
        'petrosim.models.ecafc',
        'petrosim.models.ecafc.utils'
    ],
    package_data={
        '': [
            'LICENSE',
            'README.md',
        ],
        'petrosim': [
            'models/ecafc/example.in',
            'models/ecafc/test/*'
        ]
    },
    include_package_data=True,
    install_requires = [
        "setuptools>=61.0",
        "ruamel.yaml>=0.15",
        "voluptuous>=0.13",
        "numpy>=1.19",
        "pytest>=6.2"
    ]
)