from setuptools import setup

setup(name='guacamol',
      version='0.1.4',
      author='BenevolentAI',
      author_email='guacamol@benevolent.ai',
      description='Guacamol: benchmarks for de novo molecular design',
      url='https://github.com/BenevolentAI/guacamol',
      packages=['guacamol', 'guacamol.data', 'guacamol.utils'],
      license='MIT',
      install_requires=[
          'joblib==0.12.5',
          'numpy==1.15.2',
          'tqdm==4.26.0',
          'fcd==1.0'
      ],
      dependency_links=['https://github.com/avaucher/FCD/tarball/master#egg=fcd-1.0'],
      python_requires='>=3.6',
      extras_require={
          'rdkit': ['rdkit>=2018.09.1.0'],
      },
      include_package_data=True,
      zip_safe=False,
      )
