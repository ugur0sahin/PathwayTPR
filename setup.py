from setuptools import setup, find_packages

setup(
    name='PathwayTPR',
    version='1.0.0',
    description='A description of your tool',
    author='Your name',
    author_email='your_email@example.com',
    packages=find_packages(),
    install_requires=[
        'plotly',
        'pandas'
    ],
    entry_points={
        'console_scripts': [
            'PathwayTPR=PathwayTPR.main:main'
        ]
    }
)
