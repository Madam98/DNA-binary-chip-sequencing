import distutils.core
import Cython.Build

distutils.core.setup(
    name='main-cython',
    #compiler_directives={'language_level' : "3"},
    ext_modules = Cython.Build.cythonize(
        ['main.pyx',
        'createdClasses/colors.pyx',
        'createdClasses/element.pyx',
        'createdClasses/genethicClass.pyx',
        'createdClasses/queue.pyx',
        'createdClasses/text.pyx',
        'createdFunctions/firstPopulation.pyx',
        'createdFunctions/log.pyx',
        'createdFunctions/read.pyx',
        'createdFunctions/selectionModel.pyx',
        'createdFunctions/stopCondition.pyx',
        'createdFunctions/stringOperations.pyx',
        'createdThread/threadPopulation.pyx']
    ),
    #ext_modules=Cython.Build.cythonize("*.pyx"),
    )

if __name__ == '__main__':
    import main
