from setuptools import setup


setup(name='seqalign',
      version='0.0.1',
      description='',
      url='https://github.com/uzh-dqbm-cmi/SeqAlign',
      packages=['seqalign'],
      python_requires='>=3.6.0',
      install_requires=[
            'numpy',
            'scipy',
            'matplotlib'
      ],
      zip_safe=False)